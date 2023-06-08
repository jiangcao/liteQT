!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 06/2023

MODULE wannierHam3d

IMPLICIT NONE 

private 
COMPLEX(8), ALLOCATABLE :: Hr(:,:,:,:,:)
COMPLEX(8), ALLOCATABLE :: rmn(:,:,:,:,:,:) ! position operator
COMPLEX(8), ALLOCATABLE :: pmn(:,:,:,:,:,:) ! momentum operator
COMPLEX(8), ALLOCATABLE :: Vmn(:,:,:,:,:) ! Coulomb operator
real(8), allocatable :: wannier_center(:,:)
REAL(8), DIMENSION(3) :: alpha,beta,gamm,xhat,yhat,zhat,b1,b2
REAL(8) :: Lx, Ly,Lz,cell(3,3) ! in Ang
REAL(8) :: CBM, VBM, Eg ! CB minimum and VB maximum
REAL(8) :: kt_CBM, kt_VBM ! transverse kt value of CBM and VBM location
INTEGER :: ymin,ymax,xmin,xmax,zmin,zmax,nb,nx,ny,nz,nvb, spin_deg
! neighbor index range 
! nb is the number of Wannier orbitals in a unit cell
COMPLEX(8), PARAMETER :: zzero = dcmplx(0.0d0,0.0d0)
COMPLEX(8), PARAMETER :: z1j = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: m0=5.6856D-16 ! eV s2 / cm2
REAL(8), PARAMETER :: hbar=6.58211899D-16 ! eV s
REAL(8), PARAMETER :: c0=2.998d8 ! m/s

public :: w90_free_memory, w90_load_from_file,w90_MAT_DEF
public :: NB, Lx, Ly,Lz, Nvb, VBM, CBM, kt_CBM, kt_VBM, Eg, spin_deg
public :: eig,cross,eigv,b1,b2,norm,wannier_center,alpha,beta, invert
public :: w90_MAT_DEF_full_device, w90_bare_coulomb_full_device

CONTAINS

SUBROUTINE w90_free_memory()
    deallocate(Hr)
    if (allocated(rmn)) deallocate(rmn)
    if (allocated(pmn)) deallocate(pmn)
    if (allocated(Vmn)) deallocate(Vmn)
    deallocate(wannier_center)
END SUBROUTINE w90_free_memory

SUBROUTINE w90_load_from_file(fid,lreorder_axis,axis)
implicit none
integer, intent(in) :: fid
logical, intent(in), optional :: lreorder_axis
integer, intent(in), optional :: axis(3)
integer :: n,i,nkx,nky,nkz
character(len=40) :: line, comment
REAL(8) :: dky,dkz,aux1(3),aux2(3,3)
integer, allocatable :: ind(:)
REAL(8), allocatable :: ham(:,:), energ(:,:), aux3(:,:)
    read(fid,*) nvb, spin_deg ! number of VBs, spin-degeneracy
    read(fid,*) comment
    read(fid,*) alpha
    read(fid,*) beta
    read(fid,*) gamm
    cell(:,1) = alpha
    cell(:,2) = beta
    cell(:,3) = gamm
    read(fid,*) comment
    if (present(lreorder_axis).and.(lreorder_axis)) then
        aux2(:,:) = cell(:,axis)
        cell = aux2
        alpha = cell(:,1)
        beta  = cell(:,2)
        gamm  = cell(:,3)   
    end if
    read(fid,*) n
    allocate(ham(n,7))    
    do i=1,n
        read(fid,*) ham(i,:)
    end do    
    if (present(lreorder_axis).and.(lreorder_axis)) then
        allocate(aux3(n,7))
        aux3 = ham
        aux3(:,1:3) = ham(:,axis)
        ham = aux3        
        deallocate(aux3)
    end if
    xmin=minval(ham(:,1))
    xmax=maxval(ham(:,1))
    ymin=minval(ham(:,2))
    ymax=maxval(ham(:,2))
    zmin=minval(ham(:,3))
    zmax=maxval(ham(:,3))
    nb = maxval(ham(:,4))
    nx = xmax-xmin+1
    ny = ymax-ymin+1
    nz = zmax-zmin+1
    xhat = alpha/norm(alpha)
    yhat = - cross(xhat,gamm)
    yhat = yhat/norm(yhat)
    zhat = gamm/norm(gamm)
    Ly=abs(dot_product(beta,yhat)); ! L is in unit of A
    Lx=abs(dot_product(alpha,xhat));
    Lz=norm(gamm)
    print '(A40)', 'reading Wannier H from file, info:'
    print '(8a5)', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'nb', 'nvb'
    print '(8i5)', xmin, xmax, ymin, ymax,zmin,zmax, nb, nvb
    print '(3a5)', 'Lx', 'Ly', 'Lz'
    print '(3f5.1)', Lx, Ly, Lz   
    b1=cross(beta,gamm)/dot_product(alpha,cross(beta,gamm))
    b2 = cross(gamm,alpha)/dot_product(beta,cross(gamm,alpha))    
    allocate(Hr(nb,nb,nx,ny,nz))
    Hr(:,:,:,:,:) = zzero
    do i = 1,n 
        Hr(floor(ham(i,4)),floor(ham(i,5)),floor(ham(i,1))-xmin+1,&
        &floor(ham(i,2))-ymin+1,floor(ham(i,3))-zmin+1) = ham(i,6) + z1j*ham(i,7);
    end do    
    print *, 'Find CBM and VBM'
    nkx=15
    nky=15    
    nkz=15
    if (ny .eq. 1) then
        nky = 1
    end if
    if (nx .eq. 1) then
        nkx = 1
    end if
    if (nz .eq. 1) then
        nkz = 1
    end if
    if (nky > 1) then 
        dky = 1.0d0 / dble(nky-1) * 2 * pi / Ly
    else 
        dky = 2 * pi / Ly
    end if
    if (nkz > 1) then 
        dkz = 1.0d0 / dble(nkz-1) * 2 * pi / Lz
    else 
        dkz = 2 * pi / Lz
    end if
    allocate(energ(nb,nkx*nky*nkz))
    allocate(ind(nb))
    call w90_PLOT_BZ(nkx,nky,nkz,energ)
        
    ind = maxloc(energ,2)
    kt_vbm = mod(ind(nvb),nkx)*dky - pi / Ly ! k transverse corresponding to VBMax      
    ind = minloc(energ,2)
    kt_cbm = mod(ind(nvb+1),nkx)*dky - pi / Ly ! k transverse corresponding to CBM
    
    if (ny .eq. 1) then
        kt_cbm=0.0d0
        kt_vbm=0.0d0
    end if
    ind = maxloc(energ,2)
    VBM = energ(nvb,ind(nvb))
    ind = minloc(energ,2)
    CBM = energ(nvb+1,ind(nvb+1))
    Eg = CBM-VBM
    print '(3A8)','CBM','VBM','Eg'
    print '(3f8.3)', CBM, VBM, CBM-VBM
    print '(3A8)','kt_CBM','kt_VBM','2pi/Ly'
    print '(2f8.3)', kt_cbm/(2.0*pi/Ly), kt_vbm/(2.0*pi/Ly)
    print '(A40)', 'reading Wannier centers from file'
    read(fid,*) comment    
    allocate(wannier_center(3,nb))
    do i=1,nb
        read(fid,*) wannier_center(:,i)
    end do
    if (present(lreorder_axis).and.(lreorder_axis)) then
        allocate(aux3(3,nb))
        aux3 = wannier_center
        aux3(1:3,:) = wannier_center(axis,:)
        wannier_center = aux3        
        deallocate(aux3)
    end if
    ! bring wannier_center into the first unit-cell, only works for orth.
    ! to-do: make this general
    do i=1,NB
        if (wannier_center(1,i)<0) wannier_center(:,i)=wannier_center(:,i)+alpha
        if (wannier_center(1,i)>Lx) wannier_center(:,i)=wannier_center(:,i)-alpha
        if (wannier_center(2,i)<0) wannier_center(:,i)=wannier_center(:,i)+beta
        if (wannier_center(2,i)>Ly) wannier_center(:,i)=wannier_center(:,i)-beta
        if (wannier_center(3,i)<0) wannier_center(:,i)=wannier_center(:,i)+gamm
        if (wannier_center(3,i)>Lz) wannier_center(:,i)=wannier_center(:,i)-gamm
    end do
    ! center z around 0
    wannier_center(3,:)=wannier_center(3,:)-(sum(wannier_center(3,:)/dble(NB)))
    deallocate(ham)
    deallocate(energ)
    deallocate(ind)
END SUBROUTINE w90_load_from_file


SUBROUTINE w90_PLOT_BZ(nkx,nky,nkz,EN)
implicit none
integer, intent(in) :: nkx,nky,nkz
real(8), dimension(nb,nkx*nky*nkz), intent(out) :: EN
integer :: i,j,l
real(8) :: dkx, dky, dkz,kx, ky, kz
complex(8), dimension(NB,NB) :: Hii    
    if (nkx > 1) then
        dkx = 1.0d0 / dble(nkx-1) * 2 * pi / Lx
    else 
        dkx = pi / Lx
    end if
    if (nky > 1) then
        dky = 1.0d0 / dble(nky-1) * 2 * pi / Ly
    else
        dky = pi / Ly
    end if
    if (nkz > 1) then
        dkz = 1.0d0 / dble(nkz-1) * 2 * pi / Lz
    else
        dkz = pi / Lz
    end if
    do i = 1,nkx
      do j = 1,nky
        do l = 1,nkz
          kx = i*dkx - pi / Lx
          ky = j*dky - pi / Ly
          kz = l*dky - pi / Lz
        
          call w90_MAT_DEF_3D(Hii, kx,ky,kz)    
        
          EN(1:NB,i+(j-1)*nkx+(l-1)*nkx*nky) = eig(NB,Hii)
        enddo
      enddo
    enddo
END SUBROUTINE w90_PLOT_BZ

!!! construct the diagonal and off-diagonal blocks H(I,I), H(I+1,I)
SUBROUTINE w90_MAT_DEF(Hii,H1i,kx, ky,kz,ns)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns) :: Hii, H1i
real(8), intent(in) :: ky,kx,kz
integer :: i,j,k,l
real(8), dimension(3) :: kv, r
Hii(:,:) = zzero
H1i(:,:) = zzero
do i = 1,ns
    do k = 1,ns    
        do j = ymin,ymax
          do l = zmin,zmax
            kv = kx*xhat + ky*yhat + kz*zhat
            r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
            if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) )           
            end if                 
            if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                   
                H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                & Hr(:,:,i-k-xmin+ns+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) )            
            end if
          enddo
        end do                 
    end do
end do  
END SUBROUTINE w90_MAT_DEF

!!! construct the fully periodic Hamiltonian matrix
SUBROUTINE w90_MAT_DEF_3D(Hii,kx,ky,kz)
! kx and ky is in unit of [2pi/Ang]
implicit none
REAL(8), INTENT(IN) :: kx, ky, kz
COMPLEX(8), INTENT(OUT), DIMENSION(NB,NB) :: Hii
real(8), dimension(3) :: kv, r
integer :: i,j,l
    Hii(:,:) = zzero    
    kv = kx*xhat + ky*yhat + kz*zhat
    do i = xmin,xmax           
      do j = ymin,ymax
        do l = zmin,zmax
            r = dble(i)*alpha + dble(j)*beta + dble(l)*gamm                         
            Hii(:,:) = Hii(:,:) + Hr(:,:,i-xmin+1,j-ymin+1,l-zmin+1) &
              &* exp(-z1j * dot_product(r,kv))             
        enddo
      enddo      
    enddo  
    !!! Hii = (Hii + transpose(conjg(Hii))) / 2.0d0     
END SUBROUTINE w90_MAT_DEF_3D


!!! construct the full-device Hamiltonian Matrix
SUBROUTINE w90_MAT_DEF_full_device(Ham,ky,kz,length,NS)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky,kz
complex(8), intent(inout), dimension(NB*length,NB*length) :: Ham
integer :: i,j, k,l
real(8), dimension(3) :: kv, r
complex(8) :: phi
Ham = dcmplx(0.0d0,0.0d0)
do i = 1, length
  do k = 1, length
    do j = ymin,ymax
      do l = zmin,zmax
        kv = ky*yhat + kz*zhat
        r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
        phi = dcmplx( 0.0d0, - dot_product(r,kv) )
        if (present(NS)) then
          if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
          endif                 
        else
          if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
          endif                 
        endif
      enddo
    enddo
  enddo
enddo
END SUBROUTINE w90_MAT_DEF_full_device



!!! construct the diagonal and off-diagonal blocks V(I,I), V(I+1,I)
SUBROUTINE w90_bare_coulomb_blocks(Hii,H1i,kx, ky,kz,eps,r0,ns,ldiag)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns) :: Hii, H1i
real(8), intent(in) :: ky,kx,kz
real(8), intent(in) :: eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
integer :: i,j,k,l
real(8), dimension(3) :: kv, r
Hii(:,:) = zzero
H1i(:,:) = zzero
do i = 1,ns
    do k = 1,ns    
        do j = ymin,ymax
          do l = zmin,zmax
            kv = kx*xhat + ky*yhat + kz*zhat
            r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
            if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )           
            end if                 
            if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                   
                H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )            
            end if
          enddo
        end do        
    end do
end do  
END SUBROUTINE w90_bare_coulomb_blocks

!!! construct the bare Coulomb Matrix for the full-device
SUBROUTINE w90_bare_coulomb_full_device(V,ky,kz,length,eps,r0,ldiag,NS,method)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky,kz, eps ! dielectric constant / to reduce V
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
character(len=*),intent(in),optional :: method
character(len=20) :: cmethod
complex(8), intent(out), dimension(NB*length,NB*length) :: V
integer :: i,j, k,l
real(8), dimension(3) :: kv, r
V = dcmplx(0.0d0,0.0d0)
if (present(method)) then
    cmethod=method
else
    cmethod='pointlike'
endif
select case(trim(cmethod))
  case('pointlike')
    do i = 1, length
        do k = 1, length
            do j = ymin,ymax
              do l = zmin,zmax
                kv = ky*yhat + kz*zhat
                r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                if (present(NS)) then
                    if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )           
                    end if                 
                else
                    V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                       & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )                           
                end if              
              enddo
            end do
        end do
    end do
  case('fromfile')
    call w90_load_coulomb_blocks
    do i = 1, length
        do k = 1, length
            do j = ymin,ymax
              do l = zmin,zmax
                kv = ky*yhat + kz*zhat
                r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                if (present(NS)) then
                    if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & Vmn(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) ) / eps          
                    endif                 
                else
                    if (abs(i-k)<4)then
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & Vmn(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) ) / eps                           
                    else
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )    
                    endif
                endif
              enddo
            enddo
        enddo
    enddo
end select
END SUBROUTINE w90_bare_coulomb_full_device
! function to calculate the bare coulomb potential for wannier orbitals between the (0,0) and (a1,a2) cells
FUNCTION bare_coulomb(a1,a2,a3,eps,r0,ldiag)
implicit none
integer, intent(in) :: a1, a2, a3
real(8), dimension(NB,NB) :: bare_coulomb
real(8), intent(in) :: eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
real(8), parameter :: pi=3.14159265359d0
real(8), parameter :: e=1.6d-19            ! charge of an electron (C)
real(8), parameter :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)
real(8) :: r(3),normr
real(8) :: maxV
integer :: i,j  
do i=1,NB
    do j=1,NB
        r = dble(a1)*alpha + dble(a2)*beta + dble(a3)*gamm + wannier_center(:,i) - wannier_center(:,j)
        normr = norm(r)
        if (normr >0.0d0) then
          bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*normr*1.0d-10) * tanh(normr/r0)  ! in eV
        else
          if (ldiag) then
            bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*1.0d-10) * (1.0d0/r0) ! self-interaction 
          else
            bare_coulomb(i,j) = 0.0d0
          endif
        endif
    end do
end do 
END FUNCTION bare_coulomb





FUNCTION norm(vector)
implicit none
REAL(8) :: vector(3),norm
norm = sqrt(dot_product(vector,vector))
END FUNCTION

FUNCTION cross(a, b)
implicit none
REAL(8), DIMENSION(3) :: cross
REAL(8), DIMENSION(3), INTENT(IN) :: a, b
cross(1) = a(2) * b(3) - a(3) * b(2)
cross(2) = a(3) * b(1) - a(1) * b(3)
cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

FUNCTION eig(NN, A)
implicit none
INTEGER, INTENT(IN) :: NN
COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
REAL(8) :: eig(NN)
real(8) :: W(1:NN)
integer :: INFO,LWORK,liwork, lrwork
complex(8), allocatable :: work(:)
real(8), allocatable :: RWORK(:)
!integer, allocatable :: iwork(:) 
lwork= max(1,2*NN-1)
lrwork= max(1,3*NN-2)
allocate(work(lwork))
allocate(rwork(lrwork))

CALL zheev( 'N','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )

deallocate(work,rwork)
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
   call abort
endif
eig(:)=W(:)
END FUNCTION eig

FUNCTION eigv(NN, A)
implicit none
INTEGER, INTENT(IN) :: NN
COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
REAL(8) :: eigv(NN)
real(8) :: W(1:NN)
integer :: INFO,LWORK,liwork, lrwork
complex(8), allocatable :: work(:)
real(8), allocatable :: RWORK(:)
!integer, allocatable :: iwork(:) 
lwork= max(1,2*NN-1)
lrwork= max(1,3*NN-2)
allocate(work(lwork))
allocate(rwork(lrwork))

CALL zheev( 'V','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )

deallocate(work,rwork)
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
   call abort
endif
eigv(:)=W(:)
END FUNCTION eigv

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

subroutine w90_load_coulomb_blocks()
implicit none
integer::i
real(8)::reV(NB)
real(8),parameter::ry2ev=13.6d0
real(8)::tmp(NB,NB)
allocate(Vmn(NB,NB,nx,ny,nz))
Vmn=dcmplx(0.0d0,0.0d0)

open(unit=10,file='V_CNT_0_0_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+1,1,1)=dcmplx(reV,0.0d0)*ry2ev
    Vmn(i,i,-xmin+1,1,1)=dcmplx(max(reV(i)*ry2ev,0.01),0.0d0)
enddo
close(10)
!! symmetrize V
tmp=transpose(conjg(Vmn(:,:,-xmin+1,1,1)))
Vmn(:,:,-xmin+1,1,1)=(Vmn(:,:,-xmin+1,1,1)+tmp(:,:))/2.0d0

open(unit=10,file='V_CNT_0_1_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+2,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

open(unit=10,file='V_CNT_0_2_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+3,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

open(unit=10,file='V_CNT_0_3_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+4,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

do i=1,3
    Vmn(:,:,-xmin+1-i,1,1) = transpose(conjg(Vmn(:,:,-xmin+1+i,1,1)))
enddo
end subroutine w90_load_coulomb_blocks

END MODULE wannierHam3d
