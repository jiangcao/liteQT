module invert_cpp

    use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_double_complex

    public :: init_cuda, end_cuda, invert_matrix

    interface

        integer(c_int) function init_cuda () bind (C, name = 'init_cuda')
            import c_int
        end function init_cuda

        subroutine end_cuda () bind (C, name = 'end_cuda')
        end subroutine end_cuda

        subroutine invert_matrix (A, N) bind (C, name = 'invert_matrix')
            import c_int, c_double_complex
            complex (kind=c_double_complex), intent(inout), dimension(*) :: A
            integer (kind=c_int) :: N
        end subroutine invert_matrix
  
    end  interface 

 end  module invert_cpp
