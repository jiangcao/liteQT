!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 11/2022

module green

implicit none 

private

public :: green_calc_g, green_calc_polarization

CONTAINS

! P<>(hw) = G<>(E) * G><(E+hw)
subroutine green_calc_polarization(ne,nopmin,nopmax,E,nm_dev,G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater)
implicit none
integer, intent(in) :: ne
integer, intent(in) :: nopmin,nopmax
integer, intent(in) :: nm_dev
real(8), intent(in) :: E(ne)  ! energy vector
complex(8), intent(in) :: G_retarded(nm_dev,nm_dev,ne) ! Green's functions
complex(8), intent(in) :: G_lesser(nm_dev,nm_dev,ne)
complex(8), intent(in) :: G_greater(nm_dev,nm_dev,ne)
complex(8), intent(inout) :: P_retarded(nm_dev,nm_dev,ne) ! polarization functions
complex(8), intent(inout) :: P_lesser(nm_dev,nm_dev,ne)
complex(8), intent(inout) :: P_greater(nm_dev,nm_dev,ne)
integer :: i,j,nm,ie,nop
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
complex(8), allocatable, dimension(:,:) :: B,C
complex(8) :: dE
allocate(B(nm_dev,nm_dev))
do nop=nopmin,nopmax
    P_lesser(:,:,nop) = dcmplx(0.0d0,0.0d0)
    P_greater(:,:,nop) = dcmplx(0.0d0,0.0d0)
    dE = ( E(ie) - E(ie-1) ) / 2.0d0 / pi * dcmplx(0.0d0 , -1.0d0)
    do ie = nop+1,ne 
        call zgemm('n','t',nm_dev,nm_dev,nm_dev,cone,G_lesser(:,:,ie),nm_dev,G_greater(:,:,ie-nop),nm_dev,czero,B,nm_dev) 
        P_lesser(:,:,nop) = P_lesser(:,:,nop) + B * dE
        call zgemm('n','t',nm_dev,nm_dev,nm_dev,cone,G_greater(:,:,ie),nm_dev,G_lesser(:,:,ie-nop),nm_dev,czero,B,nm_dev) 
        P_greater(:,:,nop) = P_greater(:,:,nop) + B * dE
    end do
end do
end subroutine green_calc_polarization


subroutine green_calc_g(ne,E,num_lead,nm_dev,nm_lead,max_nm_lead,Ham,H00,H10,T,G_retarded,G_lesser,G_greater,mu,temp)
implicit none
integer, intent(in) :: num_lead ! number of leads/contacts
integer, intent(in) :: nm_dev   ! size of device Hamiltonian
integer, intent(in) :: nm_lead(num_lead) ! size of lead Hamiltonians
integer, intent(in) :: max_nm_lead ! max size of lead Hamiltonians
real(8), intent(in) :: E(ne)  ! energy vector
integer, intent(in) :: ne
complex(8), intent(in) :: Ham(nm_dev,nm_dev)
complex(8), intent(in) :: H00(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian diagonal blocks
complex(8), intent(in) :: H10(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian off-diagonal blocks
complex(8), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
complex(8), intent(inout) :: G_retarded(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_lesser(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_greater(nm_dev,nm_dev,ne)
real(8), intent(in), optional :: mu(num_lead), temp
integer :: i,j,nm,ie
complex(8), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
real(8) :: fd
allocate(sig(nm_dev,nm_dev))  
do ie = 1, ne
  G_retarded(:,:,ie) = - Ham(:,:)  
  if ((present(G_lesser)).or.(present(G_greater))) then    
    if (ie .eq. 1) then
      allocate(sig_lesser(nm_dev,nm_dev))
      allocate(B(nm_dev,nm_dev))
      allocate(C(nm_dev,nm_dev))
    end if
    sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
  end if    
  ! compute and add contact self-energies    
  do i = 1,num_lead
    NM = nm_lead(i)    
    allocate(S00(nm,nm))
    allocate(G00(nm,nm))
    allocate(GBB(nm,nm))
    allocate(A(nm_dev,nm))    
    call identity(S00,nm)
    call sancho(NM,E(ie),S00,H00(1:nm,1:nm,i),H10(1:nm,1:nm,i),G00,GBB)
    call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
    G_retarded(:,:,ie) = G_retarded(:,:,ie) - sig(:,:)
    if ((present(G_lesser)).or.(present(G_greater))) then      
      fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP))		
      B(:,:) = conjg(sig(:,:))
      C(:,:) = transpose(B(:,:))
      B(:,:) = sig(:,:) - C(:,:)
      sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd	        
    end if
    deallocate(S00,G00,GBB,A)
  end do  
  do i = 1,nm_dev
    G_retarded(i,i,ie) = G_retarded(i,i,ie) + dcmplx(E(ie),0.0d0)
  end do
  !
  call invert(G_retarded(:,:,ie),nm_dev) 
  !
  if ((present(G_lesser)).or.(present(G_greater))) then    
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev) 
    if (present(G_lesser)) then
      G_lesser(:,:,ie) = C
    end if
    if (present(G_greater)) then
      B(:,:) = conjg(G_retarded(:,:,ie))
      G_greater(:,:,ie) = transpose(B(:,:))
      B(:,:) = G_retarded(:,:,ie)-G_greater(:,:,ie)
      G_greater(:,:,ie) = C + B
    end if      
  end if 
end do  
deallocate(sig)
if ((present(G_lesser)).or.(present(G_greater))) then      
  deallocate(B,C,sig_lesser)
end if
end subroutine green_calc_g


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)
implicit none
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  integer i,j,k,nm,ny,nz,nmax
  COMPLEX(8) :: z
  real(8) :: E,error
  REAL(8) :: TOL=1.0D-100  ! [eV]
  COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
  COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)

  COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:)
  COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)

  COMPLEX(8), ALLOCATABLE :: WORK(:)
  COMPLEX(8), EXTERNAL :: ZLANGE


  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( Id(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )

  Allocate( tmp(nm,nm) )

  nmax=100

  z = cmplx(E,1.0d-3)

  Id=0.0d0
  tmp=0.0d0
  do i=1,nm
     Id(i,i)=1.0d0
     tmp(i,i)=cmplx(0.0d0,1.0d0)
  enddo

  
  H_BB = H00
  H_10 = H10

  H_01 = TRANSPOSE( CONJG( H_10 ) )
  H_SS = H00

  do i = 1, nmax

	A = z*S00 - H_BB

	call invert(A,nm)

	call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm) 
	call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
	H_SS = H_SS + C
	H_BB = H_BB + C

	call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm) 
	call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm) 
	call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)  
	H_10 = C    
	H_BB = H_BB + A

	call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 

	H_01 = C 

	! NORM --> inspect the diagonal of A
	error=0.0d0
	DO k=1,nm
	 DO j=1,nm
		error=error+sqrt(aimag(C(k,j))**2+Dble(C(k,j))**2)
	 END DO
	END DO
	!write(90,*)E,i,error
	tmp=H_SS
	IF ( abs(error) < TOL ) THEN
		!rite(90,*) 'SR: Exited, abs(error)=',i,abs(error)
		EXIT
	ELSE
	END IF

	IF (i .EQ. nmax) THEN
		write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error
    call abort
	END IF
    
  enddo
  


  G00 = z*S00 - H_SS
  call invert(G00,nm)
 
  GBB = z*S00 - H_BB
  call invert(GBB,nm)



!!$  Deallocate( WORK)
!!$
  Deallocate( tmp )

  Deallocate( A )
  Deallocate( B )
  Deallocate( C )

  Deallocate( H_BB )
  Deallocate( H_SS )
  Deallocate( H_01 )
  Deallocate( H_10 )
  Deallocate( Id )

end subroutine sancho



subroutine identity(A,n)
  implicit none      
  integer, intent(in) :: n        
  complex(8), dimension(n,n), intent(inout) :: A
  integer :: i
  A = dcmplx(0.0d0,0.0d0)
  do i = 1,n
    A(i,i) = dcmplx(1.0d0,0.0d0)
  end do
end subroutine identity


subroutine invert(A,nn)
  implicit none      
  integer :: info,lda,lwork,nn      
  integer, dimension(:), allocatable :: ipiv
  complex(8), dimension(nn,nn),intent(inout) :: A
  complex(8), dimension(:), allocatable :: work
  allocate(work(nn*nn))
  allocate(ipiv(nn))
  call zgetrf(nn,nn,A,nn,ipiv,info)
  call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
  deallocate(work)
  deallocate(ipiv)
end subroutine invert

Function ferm(a)
	Real (8) a,ferm
	ferm=1.0d0/(1.0d0+Exp(a))
End Function ferm

end module green
