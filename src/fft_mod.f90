! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE. 
!
module fft_mod

implicit none 

private

public :: corr1d, conv1d, do_mkl_dfti_conv
public :: corr1d2, conv1d2

complex(8), parameter :: cone = dcmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = dcmplx(0.0d0,0.0d0)
complex(8), parameter :: c1i  = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: tpi = 3.14159265359d0*2.0d0

CONTAINS

! Z(nop) = sum_ie X(ie) Y(ie-nop)
function corr1d(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n)
complex(8)::Z(2*n-1)
complex(8),allocatable,dimension(:) :: X_, Y_
integer::i,ie
select case (trim(method))
  case default ! explicit index
    Z=czero
    do i=-n+1,n-1
      do ie = max(i+1,1),min(n,n+i)
        Z(i+n)=Z(i+n) + X(ie)*Y(ie-i)
      enddo
    enddo
  case('simple')
    do i=-n+1,n-1
      Z(i+n)=sum(X(max(1,1+i):min(n+i,n))*Y(max(1-i,1):min(n,n-i)))
    enddo
  case('fft')  
    allocate(X_(n*2-1))
    allocate(Y_(n*2-1))
    X_=czero ! pad by zero
    Y_=czero
    X_(1:n) = X
    Y_(1:n) = Y(n:1:-1)
    
    call do_mkl_dfti_conv(n*2-1,X_,Y_,Z)
    
    deallocate(X_,Y_)
end select
end function corr1d


! Z(nop) = sum_ie X(ie) Y(ie-nop)
function corr1d2(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n)
complex(8)::Z(n)
complex(8),allocatable,dimension(:) :: X_, Y_, Z_
integer::i,ie
select case (trim(method))
  case default ! explicit index
    Z=czero
    do i=-n/2+1,n/2-1
      do ie = max(i+1,1),min(n,n+i)
        Z(i+n/2)=Z(i+n/2) + X(ie)*Y(ie-i)
      enddo
    enddo
  case('simple')
    do i=-n/2+1,n/2-1
      Z(i+n/2)=sum(X(max(1,1+i):min(n+i,n))*Y(max(1-i,1):min(n,n-i)))
    enddo
  case('fft')  
    allocate(X_(n*2-1))
    allocate(Y_(n*2-1))
    allocate(Z_(n*2-1))
    X_=czero ! pad by zero
    Y_=czero
    X_(1:n) = X
    Y_(1:n) = Y(n:1:-1)
    
    call do_mkl_dfti_conv(n*2-1,X_,Y_,Z_)
    
    Z=Z_(n-n/2:n+n/2-1)
    deallocate(X_,Y_,Z_)
end select
end function corr1d2

! Z(ie) = sum_nop X(ie-nop) Y(nop)
function conv1d(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n*2-1)
complex(8)::Z(n)
complex(8),allocatable,dimension(:)::x_in, y_in, z_in
integer::i,ie
select case (trim(method))
  case default ! explicit index
    Z=czero
    do ie=1,n
      do i= -n+1,n-1
        if ((ie .gt. max(i,1)).and.(ie .lt. min(n,(n+i)))) then
          Z(ie)=Z(ie) + X(ie-i)*Y(i+n)
        endif
      enddo
    enddo
  case('simple')
    do i=1,n
      Z(i)=sum(X(n:1:-1)*Y(i:i+n-1))
    enddo
  case('fft')
    allocate(X_in(n*2-1))
    allocate(Y_in(n*2-1))
    allocate(Z_in(n*2-1))
    X_in=czero
    X_in(1:n)=X    
    Y_in=cshift(Y,-n)    
    call do_mkl_dfti_conv(n*2-1,Y_in,X_in,Z_in)
    Z=Z_in(1:n)
    deallocate(X_in,Y_in,Z_in)
end select
end function conv1d


! Z(ie) = sum_nop X(ie-nop) Y(nop)
function conv1d2(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n)
complex(8)::Z(n)
complex(8),allocatable,dimension(:)::x_in, y_in, z_in,tmp
integer::i,ie,iop
select case (trim(method))
  case default ! explicit index
    Z = czero  
    do i= -n/2+1,n/2-1  
      iop = i+n/2
      Z((max(i,1)+1):min(n,(n+i)))=Z(max(i,1)+1:min(n,(n+i))) + X((max(i,1)+1-i):(min(n,(n+i))-i))*Y(iop)
!      do ie=max(i,1)+1,min(n,(n+i))      
!        !if ((ie .gt. max(i,1)).and.(ie .lt. (n+i))) then
!          Z(ie)=Z(ie) + X(ie-i)*Y(iop)
!        !endif
!      enddo
    enddo    
  case('simple')
    allocate(Y_in(n*2-1))    
    Y_in=czero
    Y_in(n/2+1:n/2+n-1)=Y(1:n-1)    
    do i=1,n
      Z(i) = sum( Y_in(i:i+n-1)*X(n:1:-1) )
    enddo
    deallocate(Y_in)
  case('fft')
    allocate(X_in(n*2-1))
    allocate(Y_in(n*2-1))
    allocate(Z_in(n*2-1))    
    X_in=czero
    Y_in=czero
    X_in(1:n)=X    
    Y_in(1:n/2)=Y(n/2:n-1)
    Y_in(n*2-n/2+1:n*2-1)=Y(1:n/2-1)
     
    call do_mkl_dfti_conv(n*2-1,Y_in,X_in,Z_in)
    
    Z=Z_in(1:n)
    deallocate(X_in,Y_in,Z_in)
end select
end function conv1d2

subroutine do_mkl_dfti_conv(n,X_in,Y_in,Z_out)
! 1D complex to complex
Use MKL_DFTI
integer :: n
Complex(8) :: X_in(n),Y_in(n),Z_out(n)
Complex(8) :: X_out(n),Y_out(n),Z_in(n)
type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
Integer :: Status
! Perform a complex to complex transform
Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n )
Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
Status = DftiCommitDescriptor( My_Desc1_Handle )
Status = DftiComputeForward( My_Desc1_Handle, X_in, X_out )
Status = DftiComputeForward( My_Desc1_Handle, Y_in, Y_out )
!
Z_in(:) = X_out(:) * Y_out(:)
!
Status = DftiComputeBackward( My_Desc1_Handle, Z_in, Z_out )
Status = DftiFreeDescriptor(My_Desc1_Handle)
Z_out(:) = Z_out(:) / dble(n)
end subroutine do_mkl_dfti_conv

end module fft_mod
