#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <complex>    // std::complex
#include <cstddef>    // size_t
#include <stdexcept>  // std::runtime_error
#include <string>     // std::to_string
#include <unordered_map>

typedef std::complex<double> complex128;

typedef cudaStream_t gpuStream_t;
typedef cudaEvent_t gpuEvent_t;
typedef cudaError_t gpuError_t;

struct Context {
  int num_streams;
  int num_events;
  gpuStream_t *streams;
  gpuEvent_t *events;
  gpuError_t lasterror;
  Context(int nstreams, int nevents)
      : num_streams(nstreams), num_events(nevents), lasterror((gpuError_t)0) {
    streams = new gpuStream_t[nstreams];
    events = new gpuEvent_t[nevents];
  }
  ~Context() {
    delete[] streams;
    delete[] events;
  }
};

Context *gpu_context = nullptr;

extern "C" {

int init_cuda() {
    int count;

    // Check that we are able to run cuda code
    if (cudaGetDeviceCount(&count) != cudaSuccess)
    {
        printf("ERROR: GPU drivers are not configured or cuda-capable device "
               "not found\n");
        return 1;
    }
    if (count == 0)
    {
        printf("ERROR: No cuda-capable devices found\n");
        return 2;
    }

    // Initialize cuda before we run the application
    float *dev_X;
    cudaMalloc((void **) &dev_X, 1);
    cudaFree(dev_X);

    

    gpu_context = new Context(15, 16);

    // Create cuda streams and events
    for(int i = 0; i < 15; ++i) {
        cudaStreamCreateWithFlags(&gpu_context->streams[i], cudaStreamNonBlocking);
    }
    for(int i = 0; i < 16; ++i) {
        cudaEventCreateWithFlags(&gpu_context->events[i], cudaEventDisableTiming);
    }

    

    return 0;
}

void exit_cuda() {
    

    // Destroy cuda streams and events
    for(int i = 0; i < 15; ++i) {
        cudaStreamDestroy(gpu_context->streams[i]);
    }
    for(int i = 0; i < 16; ++i) {
        cudaEventDestroy(gpu_context->events[i]);
    }

    delete gpu_context;
}

bool gpu_set_stream(int streamid, gpuStream_t stream)
{
    if (streamid < 0 || streamid >= 15)
        return false;

    gpu_context->streams[streamid] = stream;

    return true;
}

void gpu_set_all_streams(gpuStream_t stream)
{
    for (int i = 0; i < 15; ++i)
        gpu_context->streams[i] = stream;
}

static void CheckCusolverDnError(cusolverStatus_t const& status) {
  if (status != CUSOLVER_STATUS_SUCCESS) {
    throw std::runtime_error("cuSOLVER failed with error code: " +
                             std::to_string(status));
  }
}

static cusolverDnHandle_t CreateCusolverDnHandle(int device) {
  if (cudaSetDevice(device) != cudaSuccess) {
    throw std::runtime_error("Failed to set CUDA device.");
  }
  cusolverDnHandle_t handle;
  CheckCusolverDnError(cusolverDnCreate(&handle));
  return handle;
}

/**
 * CUSOLVERDN wrapper class for DaCe. Once constructed, the class can be used to
 * get or create a CUSOLVERDN library handle (cusolverDnHandle_t) for a given
 * GPU ID. The class is constructed when the CUSOLVERDN DaCe library is used.
 **/
class CusolverDnHandle {
 public:
  CusolverDnHandle() = default;
  CusolverDnHandle(CusolverDnHandle const&) = delete;

  cusolverDnHandle_t& Get(int device) {
    auto f = handles_.find(device);
    if (f == handles_.end()) {
      // Lazily construct new cuSolverDn handle if the specified key does not
      // yet exist
      auto handle = CreateCusolverDnHandle(device);
      f = handles_.emplace(device, handle).first;
    }
    return f->second;
  }

  ~CusolverDnHandle() {
    for (auto& h : handles_) {
      CheckCusolverDnError(cusolverDnDestroy(h.second));
    }
  }

  CusolverDnHandle& operator=(CusolverDnHandle const&) = delete;

  std::unordered_map<int, cusolverDnHandle_t> handles_;
};

CusolverDnHandle cusolverDn_handle;

void invert_matrix(complex128 * __restrict__ _ain, int N) {

    {
        complex128 * _ainout;
        cudaMalloc((void**)&_ainout, N * N * sizeof(complex128));
        int * _pivots;
        cudaMalloc((void**)&_pivots, N * sizeof(int));
        int * _info;
        cudaMalloc((void**)&_info, 1 * sizeof(int));

        complex128 * _aout = new complex128[N * N];
        complex128 * _dev_aout;
        cudaMalloc((void**)&_dev_aout, N * N * sizeof(complex128));

        std::fill(_aout, _aout + N * N, std::complex<double>{0.0, 0.0});
        for (auto i = 0; i < N; ++i) {
            _aout[i * N + i] = 1;
        }

        cudaMemcpyAsync(_dev_aout, _aout, N * N * sizeof(complex128), cudaMemcpyHostToDevice, gpu_context->streams[0]);
        cudaMemcpyAsync(_ainout, _ain, N * N * sizeof(complex128), cudaMemcpyHostToDevice, gpu_context->streams[1]);
        {
            complex128 * _xin = &_ainout[0];
            int* _ipiv = _pivots;
            int* _res = _info;
            complex128* _xout = _ainout;

            ///////////////////
            int __dace_current_stream_id = 1;
            cudaStream_t __dace_current_stream = gpu_context->streams[__dace_current_stream_id];
            const int __dace_cuda_device = 0;
            cusolverDnHandle_t &__dace_cusolverDn_handle = cusolverDn_handle.Get(__dace_cuda_device);
            cusolverDnSetStream(__dace_cusolverDn_handle, __dace_current_stream);

            int __dace_workspace_size = 0;
            cuDoubleComplex* __dace_workspace;
            cusolverDnZgetrf_bufferSize(
            __dace_cusolverDn_handle, N, N, (cuDoubleComplex*)_xin,
            1000, &__dace_workspace_size);
            cudaMalloc<cuDoubleComplex>(
            &__dace_workspace,
            sizeof(cuDoubleComplex) * __dace_workspace_size);
            cusolverDnZgetrf(
            __dace_cusolverDn_handle, N, N, (cuDoubleComplex*)_xin,
            1000, __dace_workspace, _ipiv, _res);
            cudaFree(__dace_workspace);

            cudaEventRecord(gpu_context->events[3], gpu_context->streams[1]);
            cudaStreamWaitEvent(gpu_context->streams[3], gpu_context->events[3], 0);
            cudaEventRecord(gpu_context->events[5], gpu_context->streams[1]);
            cudaStreamWaitEvent(gpu_context->streams[4], gpu_context->events[5], 0);
            ///////////////////

        }

        cudaEventRecord(gpu_context->events[2], gpu_context->streams[0]);
        cudaStreamWaitEvent(gpu_context->streams[1], gpu_context->events[2], 0);
        {
            complex128 * _a = &_ainout[0];
            int * _ipiv = &_pivots[0];
            complex128 * _rhs_in = &_aout[0];
            int* _res = _info;
            complex128* _rhs_out = _aout;

            ///////////////////
            int __dace_current_stream_id = 1;
            cudaStream_t __dace_current_stream = gpu_context->streams[__dace_current_stream_id];
            const int __dace_cuda_device = 0;
            cusolverDnHandle_t &__dace_cusolverDn_handle = cusolverDn_handle.Get(__dace_cuda_device);
            cusolverDnSetStream(__dace_cusolverDn_handle, __dace_current_stream);

            cusolverDnZgetrs(
            __dace_cusolverDn_handle, CUBLAS_OP_N, N, N,
            (cuDoubleComplex*)_a, N, _ipiv, (cuDoubleComplex*)_rhs_in, N, _res);

            cudaEventRecord(gpu_context->events[4], gpu_context->streams[1]);
            cudaStreamWaitEvent(gpu_context->streams[2], gpu_context->events[4], 0);
            ///////////////////

        }
        cudaStreamSynchronize(gpu_context->streams[1]);
        cudaStreamSynchronize(gpu_context->streams[2]);
        cudaStreamSynchronize(gpu_context->streams[3]);

        cudaMemcpy(_aout, _dev_aout, N * N * sizeof(complex128), cudaMemcpyDeviceToHost);

        cudaFree(_ainout);
        cudaFree(_pivots);
        cudaFree(_info);
        cudaFree(_dev_aout);
        delete[] _aout;

    }
    
}

}