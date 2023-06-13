PROGRAM main
USE wannierHam3d, only : NB, w90_load_from_file, w90_free_memory,Ly,Lz,CBM,VBM,eig,w90_MAT_DEF_full_device, invert, Lx, w90_bare_coulomb_full_device,kt_CBM,spin_deg,Eg,w90_MAT_DEF,w90_bare_coulomb_blocks,w90_momentum_blocks
use green, only : green_calc_g, green_solve_gw_3D
use green_rgf, only : green_rgf_solve_gw_1d,green_RGF_CMS,green_rgf_solve_gw_3d,green_rgf_solve_gw_ephoton_3d
implicit none
real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,nky,i,j,k,axis,num_B,ib,ncpu,xyz(3),length,num_vac,iky,ik
integer,allocatable:: orb_vac(:)
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),H01(:,:),tmpV(:,:)
complex(8), allocatable,dimension(:,:,:,:) :: Hii,H1i,Vii,V1i
complex(8), allocatable,dimension(:,:,:,:,:) :: Pii,P1i
complex(8), allocatable :: H00ld(:,:,:,:),H10ld(:,:,:,:),T(:,:,:,:),V(:,:,:),Ham(:,:,:),invV(:,:,:)
complex(8), allocatable,dimension(:,:,:,:) :: G_retarded,G_lesser,G_greater
complex(8), allocatable,dimension(:,:,:,:) :: P_retarded,P_lesser,P_greater
complex(8), allocatable,dimension(:,:,:,:) :: W_retarded,W_lesser,W_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded,Sig_lesser,Sig_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8), allocatable :: Pmn(:,:,:,:)

logical :: reorder_axis, ltrans, lreadpot, lqdot, lkz, lephot, lnogw, lrgf,ldiag,lrcoulomb,labs
integer :: nen
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)
integer :: nm_dev, iter, niter, nkz,ikz,ndiag,nk
real(8) :: eps_screen, mud,mus,temps,tempd, alpha_mix, dkz,kz, r0,potscale,encut(2),dky
real(8) :: intensity,hw,midgap(2),polaris(3)
num_vac=0 ! number of vacancies
open(unit=10,file='input',status='unknown')
read(10,*) ns
read(10,*) width
read(10,*) nkx
read(10,*) num_B
allocate(B(num_B))
read(10,*) B
if (B(1)>0) then
  num_vac=floor(B(1))
  allocate(orb_vac(num_vac))
  read(10,*) orb_vac ! orbit index 
endif
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
    if (lreadpot) then
      read(10,*) potscale
    endif
    read(10,*) niter
    read(10,*) eps_screen
    read(10,*) r0
    read(10,*) mus,mud
    read(10,*) temps, tempd
    read(10,*) alpha_mix
    read(10,*) encut(1:2) ! intraband interband cutoff energies for P and W
    read(10,*) lrcoulomb
    read(10,*) ldiag
    read(10,*) lrgf
    read(10,*) lkz
    if (lkz) then
      read(10,*) nky,nkz
    endif
    read(10,*) lephot
    if (lephot) then
      read(10,*) hw,intensity
      read(10,*) polaris(:)
      read(10,*) lnogw
      read(10,*) labs
    endif
end if
close(10)

open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10, reorder_axis,xyz)
close(10)

call omp_set_num_threads(ncpu)

if (ltrans) then    
  if (.not.lrgf) then    
    if (lkz) then ! 2d case
      print *, 'Build the full device H'
      print *, 'length=',length
      allocate(Ham(nb*length,nb*length,nky*nkz))
      allocate(  V(nb*length,nb*length,nky*nkz))      
      allocate(pot(length))
      allocate(H00ld(nb*NS,nb*NS,2,nky*nkz))
      allocate(H10ld(nb*NS,nb*NS,2,nky*nkz))
      allocate(T(nb*ns,nb*length,2,nky*nkz))
      T = dcmplx(0.0d0,0.0d0)
      if (nkz>1) then
        dkz=2.0d0*pi/Lz / dble(nkz)
      else
        dkz=pi/Lz
      endif
      if (nky>1) then
        dky=2.0d0*pi/Ly / dble(nky)
      else
        dky=pi/Ly
      endif
      ! contact Ham blocks
      open(unit=11,file='ek.dat',status='unknown')
      do iky=1,nky
        ky=-pi/Ly + dble(iky)*dky
        do ikz=1,nkz
          kz=-pi/Lz + dble(ikz)*dkz
          ik=ikz+(iky-1)*nkz
          call w90_MAT_DEF(H00ld(:,:,1,ik),H10ld(:,:,1,ik),0.0d0, ky,kz,NS)
          !
          H10ld(:,:,2,ik) = H10ld(:,:,1,ik)
          H10ld(:,:,1,ik) = transpose(conjg(H10ld(:,:,1,ik)))
          H00ld(:,:,2,ik) = H00ld(:,:,1,ik)              
          T(1:nb*ns,1:nb*ns,1,ik) = H10ld(:,:,1,ik)
          T(1:nb*ns,nb*(length-ns)+1:nb*length,2,ik) = H10ld(:,:,2,ik)    
          ! write bands into ek.dat        
          allocate(phix(nkx))
          allocate(ek(nb*ns,nkx))
          allocate(H01(nb*ns,nb*ns))
          allocate(H00(nb*ns,nb*ns))
          phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
          do i=1,nkx
              H00(:,:) = H00ld(:,:,1,ik)            
              H01(:,:) = transpose(H10ld(:,:,1,ik))
              H01(:,:) = conjg(H01(:,:))
              H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H10ld(:,:,1,ik)
              H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*H01               
              ek(:,i) = eig(NB*NS,H00)
          enddo    
          deallocate(H00,H01)    
          do j=1,nb*NS
              do i=1,nkx
                  write(11,'(4E18.8)') ky*Ly, kz*Lz, phix(i), ek(j,i)
              enddo
              write(11,*)
          enddo            
          deallocate(ek,phix)        
        enddo
      enddo
      close(11)
      pot(:) = 0.0d0
      if (lreadpot) then
          open(unit=10,file='pot_dat',status='unknown')
          do i = 1,length
              read(10,*) pot(i)
          enddo
          close(10)
      endif         
      allocate(en(nen))
      allocate(G_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(G_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(G_greater(nb*length,nb*length,nen,nky*nkz))      
      
      allocate(P_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(P_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(P_greater(nb*length,nb*length,nen,nky*nkz))    
      
      allocate(W_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(W_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(W_greater(nb*length,nb*length,nen,nky*nkz))    
      
      allocate(Sig_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_greater(nb*length,nb*length,nen,nky*nkz))
      
      allocate(Sig_retarded_new(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_lesser_new(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_greater_new(nb*length,nb*length,nen,nky*nkz))
      
      Sig_retarded = dcmplx(0.0d0,0.0d0)
      Sig_lesser = dcmplx(0.0d0,0.0d0)
      Sig_greater = dcmplx(0.0d0,0.0d0)
      
      en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin

      do iky=1,nky
        ky=-pi/Ly + dble(iky)*dky
        do ikz=1,nkz
          kz=-pi/Lz + dble(ikz)*dkz
          ik=ikz+(iky-1)*nkz
          ! device Ham matrix
          call w90_MAT_DEF_full_device(Ham(:,:,ik),ky,kz,length)      
          ! Coulomb operator
          call w90_bare_coulomb_full_device(V(:,:,ik),ky,kz,length,eps_screen,r0,ldiag)      
        enddo
      enddo
      open(unit=11,file='V.dat',status='unknown')
      do i=1, size(V,1)
          do j=1, size(V,2)
              write(11,'(2I6,2F15.4)') i,j, dble(V(i,j,nky*nkz/2+1)), aimag(V(i,j,nky*nkz/2+1))
          end do
          write(11,*)
      end do
      close(11)
      
      open(unit=11,file='Ham.dat',status='unknown')
      do i=1, size(Ham,1)
          do j=1, size(Ham,2)
              write(11,'(2I6,2F15.4)') i,j, dble(Ham(i,j,nky*nkz/2+1)), aimag(Ham(i,j,nky*nkz/2+1))
          end do
          write(11,*)
      end do
      close(11)
      
      ! add on potential    
      pot = pot*potscale
      do j = 1,length
          do ib = 1,nb
              Ham((j-1)*nb+ib,(j-1)*nb+ib,:)=Ham((j-1)*nb+ib,(j-1)*nb+ib,:) + pot(j)
          end do
      end do      
      do ib = 1,nb*NS
          H00ld(ib,ib,1,:)=H00ld(ib,ib,1,:)+pot(1)
          H00ld(ib,ib,2,:)=H00ld(ib,ib,2,:)+pot(length)
      end do
      nm_dev=nb*length    
      
      midgap=(/ (CBM+VBM)/2.0d0,(CBM+VBM)/2.0d0 /)
      midgap= midgap + (/ pot(1), pot(length)/)
      !
      call green_solve_gw_3D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
        alpha_mix,nen,En,nb,ns,nky,nkz,Ham,H00ld,H10ld,T,V,&
        G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
        W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
        Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag)

      deallocate(pot)
      deallocate(H00ld)
      deallocate(H10ld) 
      deallocate(Ham)    
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
      deallocate(Sig_retarded_new)
      deallocate(Sig_lesser_new)
      deallocate(Sig_greater_new)
      deallocate(en)    
      deallocate(V)
    endif
   else
    ! Long device, use RGF
    print *, '~~~~~~~~~~~~~~~~~ RGF ~~~~~~~~~~~~~~~~~'    
    print *, 'Build the full device H'
    print *, 'length=',length*NS,'uc =',length*NS*Lx/1.0d1,'(nm)'
    nm=NB*NS
    nk=nky*nkz
    allocate(Hii(nm,nm,length,nk))
    allocate(H1i(nm,nm,length,nk))
!    allocate(V(nm*length,nm*length,nky*nkz))   
    allocate(Vii(nm,nm,length,nk))
    allocate(V1i(nm,nm,length,nk))    
    allocate(Pii(nm,nm,3,length,nk))
    allocate(P1i(nm,nm,3,length,nk))    
    allocate(pot(length*NS))    
    if (nkz>1) then
      dkz=2.0d0*pi/Lz / dble(nkz)
    else
      dkz=pi/Lz
    endif
    if (nky>1) then
      dky=2.0d0*pi/Ly / dble(nky)
    else
      dky=pi/Ly
    endif
    open(unit=10,file='Hii.dat',status='unknown')
    open(unit=11,file='H1i.dat',status='unknown')
    open(unit=20,file='Vii.dat',status='unknown')
    open(unit=21,file='V1i.dat',status='unknown')
    do iky=1,nky
      ky=-pi/Ly + dble(iky)*dky
      do ikz=1,nkz
        kz=-pi/Lz + dble(ikz)*dkz
        ik=ikz+(iky-1)*nkz
        write(10,'(A,I6,2F15.4)') 'ik',ik,ky*Ly,kz*Lz
        write(11,'(A,I6,2F15.4)') 'ik',ik,ky*Ly,kz*Lz
        write(20,'(A,I6,2F15.4)') 'ik',ik,ky*Ly,kz*Lz
        write(21,'(A,I6,2F15.4)') 'ik',ik,ky*Ly,kz*Lz
        ! get Ham blocks
        call w90_MAT_DEF(Hii(:,:,1,ik),H1i(:,:,1,ik),0.0d0, ky,kz,NS)
        !
        call w90_bare_coulomb_blocks(Vii(:,:,1,ik),V1i(:,:,1,ik),0.0d0,ky,kz,eps_screen,r0,ns,ldiag)
        !
        call w90_momentum_blocks(Pii(:,:,:,1,ik),P1i(:,:,:,1,ik),0.0d0,ky,kz,NS,'approx')
        ! write Ham blocks
        do i=1,nm
          do j=1,nm
            write(10,'(2I6,2F15.6)') i,j,dble(Hii(i,j,1,ik)),aimag(Hii(i,j,1,ik))
            write(11,'(2I6,2F15.6)') i,j,dble(H1i(i,j,1,ik)),aimag(H1i(i,j,1,ik))
            write(20,'(2I6,2F15.6)') i,j,dble(Vii(i,j,1,ik)),aimag(Vii(i,j,1,ik))
            write(21,'(2I6,2F15.6)') i,j,dble(V1i(i,j,1,ik)),aimag(V1i(i,j,1,ik))
          enddo
        enddo
        write(10,*)
        write(11,*)
        write(20,*)
        write(21,*)
        ! build device Ham 
        do i=2,length
          Hii(:,:,i,ik)=Hii(:,:,1,ik)
          H1i(:,:,i,ik)=H1i(:,:,1,ik)
          !
          Vii(:,:,i,ik)=Vii(:,:,1,ik)
          V1i(:,:,i,ik)=V1i(:,:,1,ik)
          !
          Pii(:,:,:,i,ik)=Pii(:,:,:,1,ik)
          P1i(:,:,:,i,ik)=P1i(:,:,:,1,ik)
        enddo
      enddo
    enddo
    close(10)
    close(11)
    close(20)
    close(21)
    !
    open(unit=11,file='ek.dat',status='unknown')
    allocate(H00(nb*ns,nb*ns))
    allocate(phix(nkx))
    allocate(ek(nb*ns,nkx))
    do iky=1,nky
      ky=-pi/Ly + dble(iky)*dky
      do ikz=1,nkz
        kz=-pi/Lz + dble(ikz)*dkz
        ik=ikz+(iky-1)*nkz          
        ! write bands into ek.dat                                      
        phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
        do i=1,nkx
            H00(:,:) = Hii(:,:,1,ik)                          
            H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H1i(:,:,1,ik)
            H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*conjg(transpose(H1i(:,:,1,ik)))
            ek(:,i) = eig(NB*NS,H00)
        enddo              
        do j=1,nb*NS
            do i=1,nkx
                write(11,'(4E18.8)') ky*Ly, kz*Lz, phix(i), ek(j,i)
            enddo
            write(11,*)
        enddo                      
      enddo
    enddo
    deallocate(H00)    
    deallocate(ek,phix)            
    close(11)
    !
    pot(:) = 0.0d0
    if (lreadpot) then
      open(unit=10,file='pot_dat',status='unknown')
      do i = 1,length*NS
          read(10,*) pot(i)
      end do
      close(10)
      pot=pot*potscale
      ! add on potential    
      do ik=1,nk
        do j = 1,length
          do k = 1,NS
            do ib = 1,nb
              Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik)=Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik) + pot((j-1)*NS+k)
            enddo
          enddo
        enddo      
      enddo
    end if             
    allocate(en(nen))
    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
    
!    do iky=1,nky
!      ky=-pi/Ly + dble(iky)*dky
!      do ikz=1,nkz
!        kz=-pi/Lz + dble(ikz)*dkz
!        ik=ikz+(iky-1)*nkz
!        ! coulomb operator blocks
!        call w90_bare_coulomb_blocks(Vii(:,:,1,ik),V1i(:,:,1,ik),0.0d0,ky,kz,eps_screen,r0,ns,ldiag)
!        !
!        do i=2,length
!          Vii(:,:,i,ik)=Vii(:,:,1,ik)
!          V1i(:,:,i,ik)=V1i(:,:,1,ik)
!        enddo
!      enddo
!    enddo
!    do iky=1,nky
!      ky=-pi/Ly + dble(iky)*dky
!      do ikz=1,nkz
!          kz=-pi/Lz + dble(ikz)*dkz
!          ik=ikz+(iky-1)*nkz          
!          ! Coulomb operator
!          call w90_bare_coulomb_full_device(V(:,:,ik),ky,kz,length*NS,eps_screen,r0,ldiag)      
!      enddo
!    enddo

    if (ldiag) then
      ndiag=0
    else
      ndiag=NB*NS
    endif
    if (lephot) then
      call green_rgf_solve_gw_ephoton_3d(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg),Pii,P1i,polaris,intensity,hw,labs)
    else
      call green_rgf_solve_gw_3d(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg))
    endif
    deallocate(Hii,H1i)
    !deallocate(V)
    deallocate(Vii,V1i)
    deallocate(Pii,P1i)
    deallocate(pot)
    deallocate(en)
  endif        
end if
if (allocated(B)) deallocate(B)
if (allocated(orb_vac)) deallocate(orb_vac)
call w90_free_memory
print *, 'End of program'
END PROGRAM main
