!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 09/2021

MODULE wannierHam

IMPLICIT NONE 

private 
COMPLEX(8), ALLOCATABLE :: Hr(:,:,:,:)
COMPLEX(8), ALLOCATABLE :: rmn(:,:,:,:,:) ! position operator
COMPLEX(8), ALLOCATABLE :: pmn(:,:,:,:,:) ! momentum operator
real(8), allocatable :: wannier_center(:,:)
REAL(8), DIMENSION(3) :: alpha,beta,gamm,xhat,yhat,b1,b2
REAL(8) :: Lx, Ly,cell(3,3) ! in Ang
REAL(8) :: CBM, VBM, Eg ! CB minimum and VB maximum
REAL(8) :: kt_CBM, kt_VBM ! transverse kt value of CBM and VBM location
INTEGER :: ymin,ymax,xmin,xmax,nb,nx,ny,nvb, spin_deg
! neighbor index range 
! nb is the number of Wannier orbitals in a unit cell
COMPLEX(8), PARAMETER :: zzero = dcmplx(0.0d0,0.0d0)
COMPLEX(8), PARAMETER :: z1j = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: m0=5.6856D-16 ! eV s2 / cm2
REAL(8), PARAMETER :: hbar=6.58211899D-16 ! eV s
REAL(8), PARAMETER :: c0=2.998d8 ! m/s

public :: w90_free_memory, w90_load_from_file, w90_MAT_DEF, w90_MAT_DEF_2D,w90_MAT_DEF_2D_kv, w90_plot_bz,w90_MAT_DEF_ribbon_simple
public :: w90_plot_x, NB, Lx, Ly, Nvb, VBM, CBM, kt_CBM, kt_VBM, Eg, spin_deg
public :: eig,cross,eigv,b1,b2,norm,wannier_center,alpha,beta, invert
public :: w90_ribbon_add_peierls, w90_MAT_DEF_full_device, w90_MAT_DEF_dot,w90_dot_add_peierls
public :: w90_bare_coulomb_full_device
public :: w90_momentum_full_device

CONTAINS

SUBROUTINE w90_free_memory()
    deallocate(Hr)
    if (allocated(rmn)) deallocate(rmn)
    if (allocated(pmn)) deallocate(pmn)
    deallocate(wannier_center)
END SUBROUTINE w90_free_memory

SUBROUTINE w90_load_from_file(fid,lreorder_axis,axis)
implicit none
integer, intent(in) :: fid
logical, intent(in), optional :: lreorder_axis
integer, intent(in), optional :: axis(3)
integer :: n,i,nkx,nky
character(len=40) :: line, comment
REAL(8) :: dky,aux1(3),aux2(3,3)
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
    nb = maxval(ham(:,4))
    nx = xmax-xmin+1
    ny = ymax-ymin+1
    xhat = alpha/norm(alpha)
    yhat = - cross(xhat,gamm)
    yhat = yhat/norm(yhat)
    Ly=abs(dot_product(beta,yhat)); ! L is in unit of A
    Lx=abs(dot_product(alpha,xhat));
    print '(A40)', 'reading Wannier H from file, info:'
    print '(6a5)', 'xmin', 'xmax', 'ymin', 'ymax', 'nb', 'nvb'
    print '(6i5)', xmin, xmax, ymin, ymax, nb, nvb
    print '(2f5.1)', Lx, Ly
    print *, 'a1, a2 ='
    print '(3f5.1)', alpha
    print '(3f5.1)', beta
    b1=cross(beta,gamm)/dot_product(alpha,cross(beta,gamm))
    b2 = cross(gamm,alpha)/dot_product(beta,cross(gamm,alpha))
    print *, 'b1, b2 ='
    print '(3f5.1)', b1
    print '(3f5.1)', b2
    allocate(Hr(nb,nb,nx,ny))
    Hr(:,:,:,:) = zzero
    do i = 1,n 
        Hr(floor(ham(i,4)),floor(ham(i,5)),floor(ham(i,1))-xmin+1,&
        &floor(ham(i,2))-ymin+1) = ham(i,6) + z1j*ham(i,7);
    end do    
    print *, 'Find CBM and VBM'
    nkx=25
    nky=25    
    if (ny .eq. 1) then
        nky = 1
    end if
    if (nx .eq. 1) then
        nkx = 1
    end if
    if (nky > 1) then 
        dky = 1.0d0 / dble(nky-1) * 2 * pi / Ly
    else 
        dky = 2 * pi / Ly
    end if
    allocate(energ(nb,nkx*nky))
    allocate(ind(nb))
    call w90_PLOT_BZ(nkx,nky,energ)
    ind = maxloc(energ,2)
    kt_vbm = mod(ind(nvb),nkx)*dky - pi / Ly ! k transverse corresponding to VBM
    VBM = energ(nvb,ind(nvb))
    ind = minloc(energ,2)
    kt_cbm = mod(ind(nvb+1),nkx)*dky - pi / Ly ! k transverse corresponding to CBM
    if (ny .eq. 1) then
        kt_cbm=0.0d0
        kt_vbm=0.0d0
    end if
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
    deallocate(ham)
    deallocate(energ)
    deallocate(ind)
END SUBROUTINE w90_load_from_file

SUBROUTINE w90_PLOT_LINE(kstart,kend,nk,EN)
! k in cartesian coordinate, in [2pi/L] unit
implicit none
integer, intent(in) :: nk
real(8), dimension(3), intent(in) :: kstart, kend
real(8), dimension(nb,nk), intent(out) :: EN
integer :: i
real(8) :: dkx, dky, kx, ky
complex(8), dimension(NB,NB) :: Hii
    dkx = (kend(1) - kstart(1)) / dble(nk-1) * 2 * pi / Lx
    dky = (kend(2) - kstart(2)) / dble(nk-1) * 2 * pi / Ly
    do i = 1,nk
        kx = kstart(1) + dkx*(i-1)
        ky = kstart(2) + dky*(i-1)    
        
        call w90_MAT_DEF_2D(Hii, kx,ky)    
        
        EN(1:NB,i) = eig(NB,Hii)
    end do
END SUBROUTINE w90_PLOT_LINE

SUBROUTINE w90_PLOT_BZ(nkx,nky,EN)
implicit none
integer, intent(in) :: nkx,nky
real(8), dimension(nb,nkx*nky), intent(out) :: EN
integer :: i,j
real(8) :: dkx, dky, kx, ky
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
    do i = 1,nkx
    do j = 1,nky
        kx = i*dkx - pi / Lx
        ky = j*dky - pi / Ly
        
        call w90_MAT_DEF_2D(Hii, kx,ky)    
        
        EN(1:NB,i+(j-1)*nkx) = eig(NB,Hii)
    end do
    end do
END SUBROUTINE w90_PLOT_BZ


SUBROUTINE w90_PLOT_X(ns,nk,ky,EN)
! ky is in [2pi/Ang]
implicit none
integer, intent(in) :: nk,ns    
real(8), dimension(nb*ns,nk), intent(out) :: EN
complex(8), dimension(NB*ns,NB*ns) :: Hii,H1i,Ht
real(8) :: dkx, dky, kx, ky
integer :: i
    dkx = 2.d0 * pi / dble(nk-1) /Lx/Ns
    do i = 1,nk
        kx = - pi/Lx/Ns + dkx*(i-1)                    
        call w90_MAT_DEF(Hii, H1i,kx, ky,ns)           
        Ht(:,:) = Hii(:,:) + H1i(:,:) * exp(-z1j*kx*ns*Lx) + &
                  transpose(conjg(H1i(:,:)))*exp(+z1j*kx*ns*Lx)
        EN(:,i) = eig(NB*ns,Ht);
    end do
END SUBROUTINE w90_PLOT_X

!!! construct the fully periodic Hamiltonian matrix
SUBROUTINE w90_MAT_DEF_2D(Hii,kx,ky)
! kx and ky is in unit of [2pi/Ang]
implicit none
REAL(8), INTENT(IN) :: kx, ky
COMPLEX(8), INTENT(OUT), DIMENSION(NB,NB) :: Hii
real(8), dimension(3) :: kv, r
integer :: i,j
    Hii(:,:) = zzero    
    kv = kx*xhat + ky*yhat
    do i = xmin,xmax           
        do j = ymin,ymax
            r = i*alpha + j*beta                          
            Hii(:,:) = Hii(:,:) + Hr(:,:,i-xmin+1,j-ymin+1) &
            &* exp(-z1j * dot_product(r,kv))                                  
        end do
        !!! Hii = (Hii + transpose(conjg(Hii))) / 2.0d0    
    end do   
END SUBROUTINE w90_MAT_DEF_2D

!!! construct the fully periodic Hamiltonian matrix
SUBROUTINE w90_MAT_DEF_2D_kv(Hii,kv)
! kv is in original cartesien coordinate of W90 [1/ang]
implicit none
REAL(8), INTENT(IN) :: kv(1:3)
COMPLEX(8), INTENT(OUT), DIMENSION(NB,NB) :: Hii
real(8), dimension(3) :: r
integer :: i,j
    Hii(:,:) = zzero        
    do i = xmin,xmax           
        do j = ymin,ymax
            r = i*alpha + j*beta                                      
            Hii(:,:) = Hii(:,:) + Hr(:,:,i-xmin+1,j-ymin+1) &
            &* exp(-z1j * dot_product(r,kv))                                  
        end do
        !!! Hii = (Hii + transpose(conjg(Hii))) / 2.0d0    
    end do   
END SUBROUTINE w90_MAT_DEF_2D_kv

!!! construct the diagonal and off-diagonal blocks H(I,I), H(I,I+1)
SUBROUTINE w90_MAT_DEF(Hii,H1i,kx, ky,ns)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns) :: Hii, H1i
real(8), intent(in) :: ky,kx
integer :: i,j,k
real(8) :: phiy, phix
real(8), dimension(3) :: kv, r
Hii(:,:) = zzero
H1i(:,:) = zzero
do i = 1,ns
    do k = 1,ns    
        do j = ymin,ymax
            kv = kx*xhat + ky*yhat
            r =  (i-k)*alpha + j*beta                    
            if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1) * exp(-z1j* dot_product(r,kv) )           
            end if                 
            if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                   
                H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                & Hr(:,:,i-k-xmin+ns+1,j-ymin+1) * exp(-z1j* dot_product(r,kv) )            
            end if
        end do        
    end do
end do  
END SUBROUTINE w90_MAT_DEF

!!! construct the full-device Hamiltonian Matrix
SUBROUTINE w90_MAT_DEF_full_device(Ham,ky,length,NS)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky
complex(8), intent(inout), dimension(NB*length,NB*length) :: Ham
integer :: i,j, k
real(8), dimension(3) :: kv, r
complex(8) :: phi
Ham = dcmplx(0.0d0,0.0d0)
do i = 1, length
    do k = 1, length
        do j = ymin,ymax
            kv = ky*yhat
            r =  dble(i-k)*alpha + dble(j)*beta                    
            phi = dcmplx( 0.0d0, - dot_product(r,kv) )
            if (present(NS)) then
                if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
                    Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                    & Hr(:,:,i-k-xmin+1,j-ymin+1) * exp( phi )           
                end if                 
            else
                if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                    Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                    & Hr(:,:,i-k-xmin+1,j-ymin+1) * exp( phi )           
                end if                 
            end if
        end do
    end do
end do
END SUBROUTINE w90_MAT_DEF_full_device

!!! construct the bare Coulomb Matrix for the full-device
SUBROUTINE w90_bare_coulomb_full_device(V,ky,length,eps,r0,NS)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky, eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
complex(8), intent(out), dimension(NB*length,NB*length) :: V
integer :: i,j, k
real(8), dimension(3) :: kv, r
V = dcmplx(0.0d0,0.0d0)
do i = 1, length
    do k = 1, length
        do j = ymin,ymax
            kv = ky*yhat
            r =  dble(i-k)*alpha + dble(j)*beta                    
            if (present(NS)) then
                if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                    V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                    & bare_coulomb(i-k,j,eps,r0) * exp(-z1j* dot_product(r,kv) )           
                end if                 
            else
                V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & bare_coulomb(i-k,j,eps,r0) * exp(-z1j* dot_product(r,kv) )                           
            end if
        end do
    end do
end do
END SUBROUTINE w90_bare_coulomb_full_device
! function to calculate the bare coulomb potential for wannier orbitals between the (0,0) and (a1,a2) cells
FUNCTION bare_coulomb(a1,a2,eps,r0)
implicit none
integer, intent(in) :: a1, a2
real(8), dimension(NB,NB) :: bare_coulomb
real(8), intent(in) :: eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
real(8), parameter :: pi=3.14159265359d0
real(8), parameter :: e=1.6d-19            ! charge of an electron (C)
real(8), parameter :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)
real(8) :: r(3),normr
real(8) :: maxV
integer :: i,j  
do i=1,NB
    do j=1,NB
        r = dble(a1)*alpha + dble(a2)*beta + wannier_center(:,i) - wannier_center(:,j)
        normr = norm(r)
        if (normr >0.0d0) then
          bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*normr*1.0d-10) * tanh(normr/r0)  ! in eV
        else
          bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*1.0d-10) * (1.0d0/r0) ! self-interaction 
        endif
!!        if (norm(r) .lt. r0) then
!!            bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*(norm(r)+r0)*1.0d-10)*2.0d0
!!        else
!!            bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*norm(r)*1.0d-10);  ! in eV
!!        end if
    end do
end do 
END FUNCTION bare_coulomb


!!! construct the Ribbon structure Hamiltonian's diagonal and off-diagonal blocks
SUBROUTINE w90_MAT_DEF_ribbon_simple(Hii,H1i,nn,width,axis)
implicit none
integer, intent(in) :: nn,width,axis
complex(8), intent(inout):: Hii(NB*width,NB*width),H1i(NB*width,NB*width,nn)
integer :: i,j,k
Hii=dcmplx(0.0d0,0.0d0)
H1i=dcmplx(0.0d0,0.0d0)
if (axis .eq. 1) then
do i = 1,width
    do j = 1,width
        if ((i-j <= ymax) .and. (i-j >= ymin)) then
            Hii(((i-1)*NB+1):i*NB,((j-1)*NB+1):j*NB) = Hr(:,:,-xmin+1,j-i-ymin+1)
            do k = 1,nn
                if (k <= xmax) then
                    H1i(((i-1)*NB+1):i*NB,((j-1)*NB+1):j*NB,k) = Hr(:,:,k-xmin+1,j-i-ymin+1)
                end if
            end do
        end if
    end do
end do
else 
do i = 1,width
    do j = 1,width
        if ((i-j <= xmax) .and. (i-j >= xmin)) then
            Hii(((i-1)*NB+1):i*NB,((j-1)*NB+1):j*NB) = Hr(:,:,j-i-xmin+1,-ymin+1)
            do k = 1,nn
                if (k <= ymax) then
                    H1i(((i-1)*NB+1):i*NB,((j-1)*NB+1):j*NB,k) = Hr(:,:,j-i-xmin+1,k-ymin+1)
                end if
            end do
        end if
    end do
end do
end if
END SUBROUTINE w90_MAT_DEF_ribbon_simple


!!! construct the quantum dot structure Hamiltonian
SUBROUTINE w90_MAT_DEF_dot(Ham,dot_shape,dot_size,cell_index)
implicit none
integer, intent(in) :: dot_size(:)
character(len=*), intent(in) :: dot_shape
complex(8), intent(out), allocatable:: Ham(:,:)
integer, intent(out), allocatable, optional :: cell_index(:,:)
integer, allocatable :: include_index(:,:)
integer :: i,j,k,nm,ix,iy,jx,jy,nn
real(8) :: rr,ri(3)
select case (dot_shape)
    case ('circle')
        ! find the cells inside circle of radius dot_size(1)*Lx
        rr = dot_size(1)*Lx
        allocate(include_index(2,(dot_size(1)*8)**2))
        i=0
        do ix = -dot_size(1)*4,dot_size(1)*4
            do iy = -dot_size(1)*4,dot_size(1)*4
                ri = alpha(:)*dble(ix) + beta(:)*dble(iy)
                if (norm(ri) <= rr) then
                    i=i+1                    
                    include_index(:,i) = (/ix,iy/)                    
                end if
            end do
        end do
        if (present(cell_index)) then
            allocate(cell_index(2,i))
            cell_index(:,:) = include_index(:,1:i)
        end if
        nm = i
        allocate(Ham(nb*nm,nb*nm))
        Ham = dcmplx(0.0d0,0.0d0)
        ! build H
        do i=1,nm
            ix=include_index(1,i)
            iy=include_index(2,i)
            do j=1,nm
                jx=include_index(1,j)
                jy=include_index(2,j)
                if ((ix-jx <= xmax ) .and. (ix-jx >= xmin )) then  
                    if ((iy-jy <= ymax ) .and. (iy-jy >= ymin )) then 
                        Ham((i-1)*nb+1:i*nb,(j-1)*nb+1:j*nb) = Hr(:,:,jx-ix-xmin+1,jy-iy-ymin+1)
                    end if
                end if
            end do
        end do
        deallocate(include_index)
    case default ! simple dot
        nm = (dot_size(1)*2+1)**2
        allocate(Ham(nb*nm,nb*nm))
        Ham = dcmplx(0.0d0,0.0d0)
        
        if (present(cell_index)) then
            allocate(cell_index(2,nm))
            cell_index=0
        end if
        
        i=0
        do ix = -dot_size(1),dot_size(1)
            do iy = -dot_size(1),dot_size(1)
                i = i+1
                if (present(cell_index)) then
                    cell_index(:,i) = (/ix,iy/)
                end if
                j=0
                do jx = -dot_size(1),dot_size(1)
                    do jy = -dot_size(1),dot_size(1)
                        j = j+1
                        if ((ix-jx <= xmax ) .and. (ix-jx >= xmin )) then  
                            if ((iy-jy <= ymax ) .and. (iy-jy >= ymin )) then 
                                Ham((i-1)*nb+1:i*nb,(j-1)*nb+1:j*nb) = Hr(:,:,jx-ix-xmin+1,jy-iy-ymin+1)
                            end if
                        end if
                    end do
                end do
            end do
        end do
end select
END SUBROUTINE w90_MAT_DEF_dot

!!! add the Peierls phase factor onto the Hamiltonian
!!! gauge (-By,Bx,0)
SUBROUTINE w90_dot_add_peierls(B,Ham,cell_index)
implicit none
real(8),intent(in) :: B
complex(8),intent(inout) :: Ham(:,:)
integer, intent(in)::cell_index(:,:)
integer::nm,i,j,ix,iy,jx,jy,io,jo
real(8):: r1(3),r2(3)
complex(8) :: phiB
real(8), parameter :: e = 1.60217663e-19 ! coulombs
real(8), parameter :: hbar = 1.05457182e-34 !m2*kg*/s
nm=size(cell_index,2)
do i=1,nm
    ix=cell_index(1,i)
    iy=cell_index(2,i)
    do io=1,NB
        r1(:) = alpha(:)*dble(ix) + beta(:)*dble(iy) + wannier_center(:,io)
        do j=1,nm
            jx=cell_index(1,j)
            jy=cell_index(2,j)
            do jo=1,NB
                r2(:) = alpha(:)*dble(jx) + beta(:)*dble(jy) + wannier_center(:,jo)
                phiB = e / hbar * B * (((-r1(1)+r2(1))*0.5d-10 * (r1(2)+r2(2)))*1d-10 + ((r1(1)+r2(1))*0.5d-10 * (-r1(2)+r2(2)))*1d-10)
                Ham((i-1)*NB+io,(j-1)*NB+jo) = Ham((i-1)*NB+io,(j-1)*NB+jo) * exp(dcmplx(0.0d0,1.0d0)*phiB)
            end do
        end do
    end do
end do
END SUBROUTINE w90_dot_add_peierls

!!! add the Peierls phase factor onto the Hamiltonian
!!! Landau gauge (By,0,0)
SUBROUTINE w90_ribbon_add_peierls(B,Hii,H1i,width,nn,axis,orb_pos)
implicit none
integer, intent(in) :: nn,width
integer, intent(in), optional :: axis
real(8), intent(in) :: B ! B field strength in T
real(8), intent(in), optional :: orb_pos(2,width*nb,nn)
complex(8), intent(inout) :: Hii(NB*width,NB*width), H1i(NB*width,NB*width,nn)
integer :: i,j,k,io,jo,ko
real(8) :: r1(3),r2(3),y0,x0
complex(8) :: phiB
real(8), parameter :: e = 1.60217663e-19 ! coulombs
real(8), parameter :: hbar = 1.05457182e-34 !m2*kg*/s
if ((present(axis)).and.(axis.eq.2)) then
    x0 = dot_product( width*alpha , xhat)
    do i = 1,width
        do io = 1,NB
            r1(:) = alpha(:)*dble(i) + wannier_center(:,io)
            do j = 1,width
                do jo = 1,NB
                    r2(:) = alpha(:)*dble(j) + wannier_center(:,jo)
                    phiB = e / hbar * B * (r1(2)-r2(2))*0.5d-10 * (r1(1)+r2(1)-x0)*1d-10
                    Hii((i-1)*NB+io,(j-1)*NB+jo) = Hii((i-1)*NB+io,(j-1)*NB+jo) * exp(dcmplx(0.0d0,1.0d0)*phiB)
                    do k = 1,NN
                        do ko = 1,NB
                            r2(:) = alpha(:)*dble(j) + beta(:)*dble(k) + wannier_center(:,ko)
                            phiB = e / hbar * B * (r1(2)-r2(2))*0.5d-10 * (r1(1)+r2(1)-x0)*1d-10
                            H1i((i-1)*NB+io,(j-1)*NB+jo,k) = H1i((i-1)*NB+io,(j-1)*NB+jo,k) * exp(dcmplx(0.0d0,1.0d0)*phiB)
                        end do
                    end do
                end do
            end do
        end do
    end do      
else 
    y0 = dot_product( width*beta , yhat)/2.0d0
    do i = 1,width
        do io = 1,NB
            r1(:) = beta(:)*dble(i) + wannier_center(:,io)
            do j = 1,width
                do jo = 1,NB
                    r2(:) = beta(:)*dble(j) + wannier_center(:,jo)
                    phiB = e / hbar * B * (r1(1)-r2(1))*0.5d-10 * (r1(2)+r2(2)-y0)*1d-10
                    Hii((i-1)*NB+io,(j-1)*NB+jo) = Hii((i-1)*NB+io,(j-1)*NB+jo) * exp(dcmplx(0.0d0,1.0d0)*phiB)
                    do k = 1,NN
                        do ko = 1,NB
                            r2(:) = beta(:)*dble(j) + alpha(:)*dble(k) + wannier_center(:,ko)
                            phiB = e / hbar * B * (r1(1)-r2(1))*0.5d-10 * (r1(2)+r2(2)-y0)*1d-10
                            H1i((i-1)*NB+io,(j-1)*NB+jo,k) = H1i((i-1)*NB+io,(j-1)*NB+jo,k) * exp(dcmplx(0.0d0,1.0d0)*phiB)
                        end do
                    end do
                end do
            end do
        end do
    end do      
end if
END SUBROUTINE w90_ribbon_add_peierls


SUBROUTINE w90_momentum_full_device(Ham,ky,length,NS,method)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky
complex(8), intent(inout), dimension(NB*length,NB*length,3) :: Ham ! momentum matrix [eV]
character(len=*),intent(in)::method
integer :: i,j, k,v
real(8), dimension(3) :: kv, r
complex(8) :: phi
Ham = dcmplx(0.0d0,0.0d0)
call calc_momentum_operator(method)
do v=1,3 ! cart direction
  do i = 1, length
    do k = 1, length
      do j = ymin,ymax
        kv = ky*yhat
        r =  dble(i-k)*alpha + dble(j)*beta                    
        phi = dcmplx( 0.0d0, - dot_product(r,kv) )
        if (present(NS)) then
            if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
                Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + &
                & pmn(v,:,:,i-k-xmin+1,j-ymin+1) * exp( phi )           
            end if                 
        else
            if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + &
                & pmn(v,:,:,i-k-xmin+1,j-ymin+1) * exp( phi )           
            end if                 
        end if
      end do
    end do
  end do
enddo
END SUBROUTINE w90_momentum_full_device
!
SUBROUTINE calc_momentum_operator(method)
implicit none
character(len=*),intent(in)::method
integer::io,jo,ix,iy,mx,my,mo
real(8)::r(3)
complex(8)::pre_fact
if (not(allocated(pmn))) allocate(pmn(3,NB,NB,nx,ny)) ! [eV]
pre_fact=1.0d-8 * dcmplx(0.0d0,1.0d0) *m0/hbar*1.0d2 * c0  ! multiply light-speed so c0*pmn in energy eV 
select case (method)
  case ('approx')
  ! use wannier centers, point-like orbitals
    pmn = 0.0d0
    do io=1,NB
      do jo=1,NB
        do ix=xmin,xmax
          do iy=ymin,ymax
            r = dble(ix)*alpha + dble(iy)*beta + wannier_center(:,io) - wannier_center(:,jo)
            pmn(:,io,jo,ix-xmin+1,iy-ymin+1)=Hr(io,jo,ix-xmin+1,iy-ymin+1)*r
          enddo
        enddo
      enddo
    enddo
    pmn(:,:,:,:,:)=pmn(:,:,:,:,:)*pre_fact
  case ('exact')
  ! use position operator : im_0/hbar sum_{R'l} H_{nl}(R-R') r_{lm}(R') - r_{nl}(R-R') H_{lm}(R')
    pmn = 0.0d0
    do io=1,NB
      do jo=1,NB
        do ix=xmin,xmax
          do iy=ymin,ymax
            do mo=1,NB
              do mx=xmin,xmax
                do my=ymin,ymax
                  if (((ix-mx)>=xmin).and.((ix-mx)<=xmax).and.((iy-my)>=ymin).and.((iy-my)<=ymax)) then
                  pmn(:,io,jo,ix-xmin+1,iy-ymin+1)=pmn(:,io,jo,ix-xmin+1,iy-ymin+1)+&
                    &Hr(io,mo,ix-mx-xmin+1,iy-my-ymin+1)*rmn(:,mo,jo,mx-xmin+1,my-ymin+1)-rmn(:,io,mo,ix-mx-xmin+1,iy-my-ymin+1)*Hr(mo,jo,mx-xmin+1,my-ymin+1)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    pmn(:,:,:,:,:)=pmn(:,:,:,:,:)*pre_fact
  case default
    print *, 'Unknown method!!'
    call abort
end select
END SUBROUTINE



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

END MODULE wannierHam
