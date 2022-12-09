PROGRAM main
USE wannierHam, only : NB, w90_load_from_file, w90_free_memory,Ly, w90_MAT_DEF, CBM,VBM,eig,w90_MAT_DEF_ribbon_simple, w90_ribbon_add_peierls, w90_MAT_DEF_full_device, invert, Lx,w90_MAT_DEF_dot,w90_dot_add_peierls, w90_bare_coulomb_full_device,kt_CBM
use green, only : green_calc_g, green_calc_polarization, green_calc_w, green_calc_gw_selfenergy
use mod_string, only : string
implicit none
real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,i,j,k,axis,num_B,ib,ncpu,xyz(3),length
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),Hii(:,:),H1i(:,:,:),BHii(:,:),BH1i(:,:,:),Ham(:,:),BHam(:,:),H00ld(:,:,:),H10ld(:,:,:),T(:,:,:),V(:,:)
complex(8), allocatable :: G_retarded(:,:,:),G_lesser(:,:,:),G_greater(:,:,:)
complex(8), allocatable :: P_retarded(:,:,:),P_lesser(:,:,:),P_greater(:,:,:)
complex(8), allocatable :: W_retarded(:,:,:),W_lesser(:,:,:),W_greater(:,:,:)
complex(8), allocatable :: Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:)
complex(8)::ldos,pdos,ndos,prT,plT,pgT,wrT,wlT,wgT,srT,slT,sgT
logical :: reorder_axis, ltrans, lreadpot, lqdot
integer :: nen
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)
integer :: nm_dev, iter

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
    allocate(V(nb*length,nb*length))
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
    
    Sig_retarded = dcmplx(0.0d0,0.0d0)
    Sig_lesser = dcmplx(0.0d0,0.0d0)
    Sig_greater = dcmplx(0.0d0,0.0d0)
    
    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
    ! device Ham matrix
    call w90_MAT_DEF_full_device(Ham,kt_CBM,length)
    ! Coulomb operator
    call w90_bare_coulomb_full_device(V,0.0d0,length,1.0d0)
    !
    open(unit=11,file='V.dat',status='unknown')
    do i=1, size(V,1)
        do j=1, size(V,2)
            write(11,'(2I6,2F15.4)') i,j, dble(V(i,j)), aimag(V(i,j))
        end do
        write(11,*)
    end do
    close(11)
    !
    open(unit=11,file='Ham.dat',status='unknown')
    do i=1, size(Ham,1)
        do j=1, size(Ham,2)
            write(11,'(2I6,2F15.4)') i,j, dble(Ham(i,j)), aimag(Ham(i,j))
        end do
        write(11,*)
    end do
    close(11)
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
    
    do iter = 0,1
        !
        print *, 'calc G'
        call green_calc_g(nen,En,2,nb*length,(/nb*ns,nb*ns/),nb*ns,Ham,H00ld,H10ld,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,(/VBM,CBM/),(/300.0d0,300.0d0/))
        !
        print *, 'calc P'
        call green_calc_polarization(nen,100,150,En,nb*length,G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater)    
        !
        print *, 'calc W'
        call green_calc_w(nen,100,150,En,nm_dev,NS*NB,V,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater)
        !
        print *, 'calc SigGW'
        call green_calc_gw_selfenergy(nen,100,150,En,nm_dev,G_retarded,G_lesser,G_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater)
        !
        open(unit=11,file='ldos_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=12,file='pdos_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=13,file='ndos_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=14,file='P_r_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=15,file='P_lesser_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=16,file='P_greater_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=17,file='W_r_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=18,file='W_lesser_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=19,file='W_greater_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=20,file='Sig_r_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=21,file='Sig_lesser_'//TRIM(STRING(iter))//'.dat',status='unknown')
        open(unit=22,file='Sig_greater_'//TRIM(STRING(iter))//'.dat',status='unknown')
        do i = 1,nen
            do j = 1,length
                ldos=0.0d0
                pdos=0.0d0
                ndos=0.0d0
                prT=0.0d0
                plT=0.0d0
                pgT=0.0d0
                wrT=0.0d0
                wlT=0.0d0
                wgT=0.0d0
                srT=0.0d0
                slT=0.0d0
                sgT=0.0d0
                do ib=1,nb
                    ldos = ldos+ G_retarded((j-1)*nb+ib,(j-1)*nb+ib,i)
                    ndos = ndos+ G_lesser((j-1)*nb+ib,(j-1)*nb+ib,i)
                    pdos = pdos+ G_greater((j-1)*nb+ib,(j-1)*nb+ib,i)
                    ! P
                    prT = prT+ P_retarded((j-1)*nb+ib,(j-1)*nb+ib,i)
                    plT = plT+ P_lesser((j-1)*nb+ib,(j-1)*nb+ib,i)
                    pgT = pgT+ P_greater((j-1)*nb+ib,(j-1)*nb+ib,i)
                    ! W
                    wrT = wrT+ W_retarded((j-1)*nb+ib,(j-1)*nb+ib,i)
                    wlT = wlT+ W_lesser((j-1)*nb+ib,(j-1)*nb+ib,i)
                    wgT = wgT+ W_greater((j-1)*nb+ib,(j-1)*nb+ib,i)
                    ! Sig
                    srT = srT+ Sig_retarded((j-1)*nb+ib,(j-1)*nb+ib,i)
                    slT = slT+ Sig_lesser((j-1)*nb+ib,(j-1)*nb+ib,i)
                    sgT = sgT+ Sig_greater((j-1)*nb+ib,(j-1)*nb+ib,i)
                end do
                write(11,'(3F15.4)') j*Lx, en(i) , -aimag(ldos)
                write(12,'(3F15.4)') j*Lx, en(i) , -aimag(pdos)
                write(13,'(3F15.4)') j*Lx, en(i) , aimag(ndos)
                !
                write(14,'(4F15.4)') j*Lx, en(i)-en(1) , dble(prT), aimag(prT)
                write(15,'(4F15.4)') j*Lx, en(i)-en(1) , dble(plT), aimag(plT)
                write(16,'(4F15.4)') j*Lx, en(i)-en(1) , dble(pgT), aimag(pgT)
                !
                write(17,'(4F15.4)') j*Lx, en(i)-en(1) , dble(wrT), aimag(wrT)
                write(18,'(4F15.4)') j*Lx, en(i)-en(1) , dble(wlT), aimag(wlT)
                write(19,'(4F15.4)') j*Lx, en(i)-en(1) , dble(wgT), aimag(wgT)
                !
                write(20,'(4F15.4)') j*Lx, en(i) , dble(srT), aimag(srT)
                write(21,'(4F15.4)') j*Lx, en(i) , dble(slT), aimag(slT)
                write(22,'(4F15.4)') j*Lx, en(i) , dble(sgT), aimag(sgT)
            end do
            write(11,*)
            write(12,*)
            write(13,*)
            write(14,*)
            write(15,*)
            write(16,*)
            write(17,*)
            write(18,*)
            write(19,*)
            write(20,*)
            write(21,*)
            write(22,*)
        end do
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        close(17)
        close(18)
        close(19)
        close(20)
        close(21)
        close(22)
    end do
    !
    iter = 2
    print *, 'calc G'
    call green_calc_g(nen,En,2,nb*length,(/nb*ns,nb*ns/),nb*ns,Ham,H00ld,H10ld,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,(/VBM,CBM/),(/300.0d0,300.0d0/))
    open(unit=11,file='ldos_'//TRIM(STRING(iter))//'.dat',status='unknown')
    open(unit=12,file='pdos_'//TRIM(STRING(iter))//'.dat',status='unknown')
    open(unit=13,file='ndos_'//TRIM(STRING(iter))//'.dat',status='unknown')
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
            write(12,'(3F15.4)') j*Lx, en(i) , -aimag(pdos)
            write(13,'(3F15.4)') j*Lx, en(i) , aimag(ndos)
        end do
        write(11,*)
        write(12,*)
        write(13,*)
    end do
    close(11)
    close(12)
    close(13)    
    !
    deallocate(pot)
    deallocate(H00ld)
    deallocate(H10ld) 
    deallocate(Ham)
    deallocate(V)
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
