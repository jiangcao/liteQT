PROGRAM main
USE wannierHam, only : NB, w90_load_from_file, w90_free_memory,Ly, w90_MAT_DEF, CBM,VBM,eig,w90_MAT_DEF_ribbon_simple, w90_ribbon_add_peierls, w90_MAT_DEF_full_device, invert, Lx,w90_MAT_DEF_dot,w90_dot_add_peierls
use green, only : green_calc_g

implicit none
real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,i,j,k,axis,num_B,ib,ncpu,xyz(3),length
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),Hii(:,:),H1i(:,:,:),BHii(:,:),BH1i(:,:,:),Ham(:,:),BHam(:,:),G_retarded(:,:,:),sig(:,:),g00(:,:),gbb(:,:),s00(:,:),A(:,:),T(:,:,:),H00ld(:,:,:),H10ld(:,:,:),G_lesser(:,:,:),G_greater(:,:,:)
complex(8)::ldos,pdos,ndos
logical :: reorder_axis, ltrans, lreadpot, lqdot
integer :: nen
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)

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
read(10,*) ltrans
if (ltrans) then
    read(10,*) length
    read(10,*) emin, emax, nen
    read(10,*) lreadpot
end if
read(10,*) lqdot
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
    allocate(pot(length))
    allocate(H00ld(nb*NS,nb*NS,2))
    allocate(H10ld(nb*NS,nb*NS,2))
    allocate(T(nb*ns,nb*length,2))
    call w90_MAT_DEF(H00ld(:,:,1),H10ld(:,:,1),0.0d0, 0.0d0,NS)
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

    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin

    call w90_MAT_DEF_full_device(Ham,0.0d0,length)
        
    do j = 1,length
        do ib = 1,nb
            Ham((j-1)*nb+ib,(j-1)*nb+ib)=Ham((j-1)*nb+ib,(j-1)*nb+ib)+pot(j)
        end do
    end do
    
    do ib = 1,nb*NS
        H00ld(ib,ib,1)=H00ld(ib,ib,1)+pot(1)
        H00ld(ib,ib,2)=H00ld(ib,ib,2)+pot(length)
    end do
        
    call green_calc_g(nen,En,2,nb*length,(/nb*ns,nb*ns/),nb*ns,Ham,H00ld,H10ld,T,G_retarded,G_lesser,G_greater,(/-4.1d0,-3.5d0/),300.0d0)
        
    
    open(unit=11,file='ldos.dat',status='unknown')
    open(unit=12,file='pdos.dat',status='unknown')
    open(unit=13,file='ndos.dat',status='unknown')
    do i = 1,nen
        do j = 1,length
            ldos=0.0d0
            pdos=0.0d0
            ndos=0.0d0
            do ib=1,nb
                ldos = ldos+ G_retarded((j-1)*nb+ib,(j-1)*nb+ib,i)
                ndos = ndos+ G_lesser((j-1)*nb+ib,(j-1)*nb+ib,i)
                pdos = pdos+ G_greater((j-1)*nb+ib,(j-1)*nb+ib,i)
            end do
            write(11,'(3F15.4)') j*Lx, en(i) , -aimag(ldos)
            write(12,'(3F15.4)') j*Lx, en(i) , aimag(pdos)
            write(13,'(3F15.4)') j*Lx, en(i) , aimag(ndos)
        end do
        write(11,*)
        write(12,*)
        write(13,*)
    end do
    close(11)
    close(12)
    close(13)
    
    deallocate(pot)
    deallocate(H00ld)
    deallocate(H10ld) 
    deallocate(Ham)
    deallocate(T)
    deallocate(G_retarded)
    deallocate(G_lesser)
    deallocate(G_greater)
    deallocate(en)
    
end if

deallocate(BHii)
deallocate(BH1i)
deallocate(phix)
deallocate(B)

call w90_free_memory

END PROGRAM main
