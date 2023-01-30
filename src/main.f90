PROGRAM main
USE wannierHam, only : NB, w90_load_from_file, w90_free_memory,Ly, w90_MAT_DEF, CBM,VBM,eig,w90_MAT_DEF_ribbon_simple, w90_ribbon_add_peierls, w90_MAT_DEF_full_device, invert, Lx,w90_MAT_DEF_dot,w90_dot_add_peierls, w90_bare_coulomb_full_device,kt_CBM,w90_inverse_bare_coulomb_full_device
use green, only : green_calc_g, green_subspace_invert,green_solve_gw_1D
use mod_string, only : string
implicit none
real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,i,j,k,axis,num_B,ib,ncpu,xyz(3),length
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),Hii(:,:),H1i(:,:,:),BHii(:,:),BH1i(:,:,:),Ham(:,:),BHam(:,:),H00ld(:,:,:),H10ld(:,:,:),T(:,:,:),V(:,:),invV(:,:),V2(:,:),invV2(:,:)
complex(8), allocatable :: G_retarded(:,:,:),G_lesser(:,:,:),G_greater(:,:,:)
complex(8), allocatable :: P_retarded(:,:,:),P_lesser(:,:,:),P_greater(:,:,:)
complex(8), allocatable :: W_retarded(:,:,:),W_lesser(:,:,:),W_greater(:,:,:)
complex(8), allocatable :: Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:)
complex(8), allocatable :: Sig_retarded_new(:,:,:),Sig_lesser_new(:,:,:),Sig_greater_new(:,:,:)
complex(8)::ldos,pdos,ndos,prT,plT,pgT,wrT,wlT,wgT,srT,slT,sgT
logical :: reorder_axis, ltrans, lreadpot, lqdot
integer :: nen
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)
integer :: nm_dev, iter, niter
real(8) :: eps_screen, mud,mus,temps,tempd, alpha_mix

open(unit=10,file='input',status='unknown')
read(10,*) ns
read(10,*) width
read(10,*) nkx
read(10,*) num_B
allocate(B(num_B))
read(10,*) B
read(10,*) axis
read(10,*) ncpu
read(10,*) reorder_axis
if (reorder_axis) then
    read(10,*) xyz
end if    
read(10,*) lqdot
read(10,*) ltrans
if (ltrans) then
    read(10,*) length
    read(10,*) emin, emax, nen
    read(10,*) lreadpot
    read(10,*) niter
    read(10,*) eps_screen
    read(10,*) mus,mud
    read(10,*) temps, tempd
    read(10,*) alpha_mix
end if
close(10)

open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10, reorder_axis,xyz)
close(10)

call omp_set_num_threads(ncpu)

allocate(Hii(nb*width,nb*width))
allocate(H1i(nb*width,nb*width,ns))
allocate(BHii(nb*width,nb*width))
allocate(BH1i(nb*width,nb*width,ns))
allocate(H00(nb*width,nb*width))
allocate(ek(nb*width,nkx))

print *, 'define ribbon H'
call w90_MAT_DEF_ribbon_simple(BHii,BH1i, NS, width,axis)
open(unit=11,file='ek.dat',status='unknown')
allocate(phix(nkx))

do ib=1,num_B
    print *, 'B point=', B(ib)
    Hii = BHii
    H1i = BH1i
    call w90_ribbon_add_peierls(B(ib),Hii,H1i,width,NS,axis=axis)

    phix=(/(i, i=1,nkx, 1)/) / dble(nkx) * pi * 2.0d0 - pi
    
    !$omp parallel default(none) private(i,j,H00) shared(nkx,ns,Hii,H1i,phix,ek,NB,width)
    !$omp do
    do i=1,nkx
        H00(:,:) = Hii(:,:)
        do j=1,ns
            H00(:,:) = H00(:,:) + exp(+dcmplx(0.0d0,1.0d0)*phix(i)*dble(j))*H1i(:,:,j)+ exp(-dcmplx(0.0d0,1.0d0)*phix(i)*dble(j))*conjg(transpose(H1i(:,:,j)))        
        end do
        ek(:,i) = eig(NB*width,H00)
    end do
    !$omp end do
    !$omp end parallel

    do j=1,nb*width
        do i=1,nkx
            write(11,'(3F15.4)') phix(i), ek(j,i) , B(ib)   
        end do
        write(11,*)
    end do    

end do
close(11)
deallocate(Hii)
deallocate(H1i)
deallocate(H00)
deallocate(ek)

if (lqdot) then
    print *, 'Build the quantum dot H'
    print *, 'size=',width
    call w90_MAT_DEF_dot(Ham,'circle',(/width/2/),cell_index)
    
    allocate(BHam(size(Ham,1),size(Ham,2)))
    allocate(ek(size(Ham,1),num_B))
    
    do ib=1,num_B
        print *, 'B point=', B(ib)
        BHam(:,:) = Ham(:,:)
        call w90_dot_add_peierls(B(ib),BHam,cell_index)    
        ek(:,ib) = eig(size(BHam,1),BHam)        
    end do
    
    open(unit=11,file='dot_en.dat',status='unknown')
    do j=1, size(Ham,1)
        do ib=1,num_B
            write(11,'(2F15.4)') B(ib), ek(j,ib) 
        end do
        write(11,*)
    end do
    close(11)
    deallocate(Ham,ek,cell_index)
    deallocate(BHam)
end if

if (ltrans) then
    print *, 'Build the full device H'
    print *, 'length=',length
    allocate(Ham(nb*length,nb*length))
    allocate(   V(nb*length*3,nb*length*3))
    allocate(invV(nb*length*3,nb*length*3))
    allocate(V2(nb*length,nb*length))
    allocate(invV2(nb*length,nb*length))
    allocate(pot(length))
    allocate(H00ld(nb*NS,nb*NS,2))
    allocate(H10ld(nb*NS,nb*NS,2))
    allocate(T(nb*ns,nb*length,2))
    ! contact Ham blocks
    call w90_MAT_DEF(H00ld(:,:,1),H10ld(:,:,1),0.0d0, kt_CBM,NS)
    !
    H10ld(:,:,2) = H10ld(:,:,1)
    H10ld(:,:,1) = transpose(conjg(H10ld(:,:,1)))
    H00ld(:,:,2) = H00ld(:,:,1)
    
    T = dcmplx(0.0d0,0.0d0)
    T(1:nb*ns,1:nb*ns,1) = H10ld(:,:,1)
    T(1:nb*ns,nb*(length-ns)+1:nb*length,2) = H10ld(:,:,2)
    
    pot(:) = 0.0d0
    if (lreadpot) then
        open(unit=10,file='pot_dat',status='unknown')
        do i = 1,length
            read(10,*) pot(i)
        end do
        close(10)
    end if    

    allocate(en(nen))
    allocate(G_retarded(nb*length,nb*length,nen))
    allocate(G_lesser(nb*length,nb*length,nen))
    allocate(G_greater(nb*length,nb*length,nen))
    
    allocate(P_retarded(nb*length,nb*length,nen))
    allocate(P_lesser(nb*length,nb*length,nen))
    allocate(P_greater(nb*length,nb*length,nen))
    
    allocate(W_retarded(nb*length,nb*length,nen))
    allocate(W_lesser(nb*length,nb*length,nen))
    allocate(W_greater(nb*length,nb*length,nen))
    
    allocate(Sig_retarded(nb*length,nb*length,nen))
    allocate(Sig_lesser(nb*length,nb*length,nen))
    allocate(Sig_greater(nb*length,nb*length,nen))
    
    allocate(Sig_retarded_new(nb*length,nb*length,nen))
    allocate(Sig_lesser_new(nb*length,nb*length,nen))
    allocate(Sig_greater_new(nb*length,nb*length,nen))
    
    Sig_retarded = dcmplx(0.0d0,0.0d0)
    Sig_lesser = dcmplx(0.0d0,0.0d0)
    Sig_greater = dcmplx(0.0d0,0.0d0)
    
    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
    
    ! device Ham matrix
    call w90_MAT_DEF_full_device(Ham,kt_CBM,length)
    
    ! Coulomb operator
    call w90_bare_coulomb_full_device(V,0.0d0,length*3,eps_screen)
    
    open(unit=11,file='V.dat',status='unknown')
    do i=1, size(V,1)
        do j=1, size(V,2)
            write(11,'(2I6,2F15.4)') i,j, dble(V(i,j)), aimag(V(i,j))
        end do
        write(11,*)
    end do
    close(11)
    !
    !open(unit=11,file='Ham.dat',status='unknown')
    !do i=1, size(Ham,1)
    !    do j=1, size(Ham,2)
    !        write(11,'(2I6,2F15.4)') i,j, dble(Ham(i,j)), aimag(Ham(i,j))
    !    end do
    !    write(11,*)
    !end do
    !close(11)
    ! add on potential    
    do j = 1,length
        do ib = 1,nb
            Ham((j-1)*nb+ib,(j-1)*nb+ib)=Ham((j-1)*nb+ib,(j-1)*nb+ib)+pot(j)
        end do
    end do
    
    do ib = 1,nb*NS
        H00ld(ib,ib,1)=H00ld(ib,ib,1)+pot(1)
        H00ld(ib,ib,2)=H00ld(ib,ib,2)+pot(length)
    end do
    nm_dev=nb*length
    
!    invV=V
!    call green_subspace_invert(nm_dev*3,invV,NS*nb*2,'sancho')

!    call w90_inverse_bare_coulomb_full_device(invV2,0.0d0,length,eps_screen,20,1)    
!    open(unit=11,file='inv_V2.dat',status='unknown')
!    do i=1, size(invV2,1)
!        do j=1, size(invV2,2)
!            write(11,'(2I6,2E18.4)') i,j, dble(invV2(i,j)), aimag(invV2(i,j))
!        end do
!        write(11,*)
!    end do
!    close(11)    
    
    V2=V(nm_dev+1:nm_dev*2,nm_dev+1:nm_dev*2)
 !   invV2=invV(nm_dev+1:nm_dev*2,nm_dev+1:nm_dev*2)
    
    deallocate(V)
    deallocate(invV)

    !open(unit=11,file='inv_V.dat',status='unknown')
    !do i=1, size(invV2,1)
    !    do j=1, size(invV2,2)
    !        write(11,'(2I6,2E18.4)') i,j, dble(invV2(i,j)), aimag(invV2(i,j))
    !    end do
    !    write(11,*)
    !end do
    !close(11)    
    !
    call green_solve_gw_1D(niter,nm_dev,Lx,length,temps,tempd,mus,mud,&
        alpha_mix,nen,En,nb,ns,Ham,H00ld,H10ld,T,V2,&
        G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
        W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
        Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
    
    deallocate(pot)
    deallocate(H00ld)
    deallocate(H10ld) 
    deallocate(Ham)
    deallocate(V2)
    deallocate(invV2)
    deallocate(T)
    deallocate(G_retarded)
    deallocate(G_lesser)
    deallocate(G_greater)
    deallocate(P_retarded)
    deallocate(P_lesser)
    deallocate(P_greater)
    deallocate(W_retarded)
    deallocate(W_lesser)
    deallocate(W_greater)
    deallocate(Sig_retarded)
    deallocate(Sig_lesser)
    deallocate(Sig_greater)
    deallocate(en)
    
end if

deallocate(BHii)
deallocate(BH1i)
deallocate(phix)
deallocate(B)

call w90_free_memory

END PROGRAM main
