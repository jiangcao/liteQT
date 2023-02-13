!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 11/2022

module green

implicit none 

private

public :: green_calc_g
public :: green_solve_gw_1D,green_solve_gw_2D
public :: green_solve_ephoton_freespace_1D

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
real(8), parameter :: m0=9.109d-31 ! kg
real(8), parameter :: eps0=8.854d-12 ! C/V/m 
real(8), parameter :: c0=2.998d8 ! m/s
real(8), parameter :: e0=1.6022d-19 ! C

CONTAINS

subroutine green_solve_ephoton_freespace_1D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,&
  Pmn,polarization,intensity,hw,&
  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
implicit none
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, Lx,spindeg,alpha_mix 
real(8), intent(in) :: hw ! hw is photon energy in eV
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: Pmn(nm_dev,nm_dev,3) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
real(8), intent(in) :: polarization(3) ! light polarization vector 
real(8), intent(in) :: intensity ! [W/m^2]
real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
!---------
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
integer::ie,nop,i,iter
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable::M(:,:) ! e-photon coupling matrix
real(8)::Nphot,mu(2)
  print *,'====== green_solve_ephoton_freespace_1D ======'
  allocate(tot_cur(nm_dev,nm_dev))
  allocate(tot_ecur(nm_dev,nm_dev))
  allocate(cur(nm_dev,nm_dev,nen))
  allocate(Ispec(nm_dev,nm_dev,nen))
  allocate(Itot(nm_dev,nm_dev))
  mu=(/ mus, mud /)
  allocate(siglead(NB*NS,NB*NS,nen,2))
  siglead=dcmplx(0.0d0,0.0d0)  
  allocate(M(nm_dev,nm_dev))
  M=dcmplx(0.0d0,0.0d0)
  print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
  do i=1,3
    M=M+ polarization(i) * Pmn(:,:,i) 
  enddo    
  !print *, 'pre_fact=', pre_fact
  print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  nop=floor(hw / (En(2)-En(1)))
  print *,'nop=',nop
  print '(a8,f15.4)','dE(meV)=',(En(2)-En(1))*1.0d3
  do iter=0,niter
    ! empty files for sancho 
    open(unit=101,file='sancho_gbb.dat',status='unknown')
    close(101)
    open(unit=101,file='sancho_g00.dat',status='unknown')
    close(101)
    open(unit=101,file='sancho_sig.dat',status='unknown')
    close(101)
    print *, 'calc G'  
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,mu,(/temps,tempd/))    
    call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    call write_current_spectrum('Jdens',iter,cur,nen,en,length,NB,Lx)
    call write_current('I',iter,tot_cur,length,NB,NS*2,Lx)
    call write_current('EI',iter,tot_ecur,length,NB,NS*2,Lx)
    call write_spectrum('ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0,-2.0/))
    call write_spectrum('ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
    call write_spectrum('pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0,-1.0/))
    !
    print *, 'calc Sig'
    call calc_sigma_ephoton_monochromatic(nm_dev,length,nen,En,nop,M,G_lesser,G_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
    Sig_retarded_new=Sig_retarded_new*pre_fact*intensity/hw**2
    Sig_greater_new=Sig_greater_new*pre_fact*intensity/hw**2
    Sig_lesser_new=Sig_lesser_new*pre_fact*intensity/hw**2
    ! mixing with the previous one
    Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
    Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
    Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
    ! get leads sigma
    siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
    siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
    call write_spectrum('SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0,1.0/))
    call write_spectrum('SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
    call write_spectrum('SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0,1.0/))
    ! calculate collision integral
    call calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
    call write_spectrum('Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0,1.0/))
  enddo
  deallocate(M,siglead)
  deallocate(cur,tot_cur,tot_ecur)
  deallocate(Ispec,Itot)
end subroutine green_solve_ephoton_freespace_1D
! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic(nm_dev,length,nen,En,nop,M,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater)
implicit none
integer,intent(in)::nm_dev,length,nen,nop
real(8),intent(in)::en(nen)
complex(8),intent(in),dimension(nm_dev,nm_dev)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
  Sig_lesser=0.0d0
  Sig_greater=0.0d0
  Sig_retarded=0.0d0  
  ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
  !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
  !$omp parallel default(none) private(ie,A,B) shared(nop,nen,nm_dev,G_lesser,G_greater,Sig_lesser,Sig_greater,M)
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev))  
  !$omp do
  do ie=1,nen
    ! Sig^<(E)
    A = 0.0d0
    if (ie-nop>=1) A =A+ G_lesser(:,:,ie-nop)
    if (ie+nop<=nen) A =A+ G_lesser(:,:,ie+nop)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,A,nm_dev,czero,B,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
    Sig_lesser(:,:,ie) = A    
    ! Sig^>(E)
    A = 0.0d0
    if (ie-nop>=1) A =A+ G_greater(:,:,ie-nop)
    if (ie+nop<=nen) A =A+ G_greater(:,:,ie+nop)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,A,nm_dev,czero,B,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
    Sig_greater(:,:,ie) = A    
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
end subroutine calc_sigma_ephoton_monochromatic

! 2D GW solver with one periodic direction (z)
! iterating G -> P -> W -> Sig 
subroutine green_solve_gw_2D(niter,nm_dev,Lx,length,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,nphiz,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
implicit none
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length, nphiz
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx
complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiz),H00lead(NB*NS,NB*NS,2,nphiz),H10lead(NB*NS,NB*NS,2,nphiz),T(NB*NS,nm_dev,2,nphiz)
complex(8), intent(in):: V(nm_dev,nm_dev,nphiz)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen,nphiz) ::  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
integer :: iter,ie,nopmax
integer :: i,j,nm,nop,l,h,iop,ndiag,ikz,iqz,ikzd
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
allocate(siglead(NB*NS,NB*NS,nen,2,nphiz))
siglead=dcmplx(0.0d0,0.0d0)
allocate(B(nm_dev,nm_dev))
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
do iter=0,niter
  ! empty files for sancho 
  open(unit=101,file='sancho_gbb.dat',status='unknown')
  close(101)
  open(unit=101,file='sancho_g00.dat',status='unknown')
  close(101)
  open(unit=101,file='sancho_sig.dat',status='unknown')
  close(101)
  print *, 'calc G'  
  do ikz=1,nphiz
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),G_retarded(:,:,:,ikz),G_lesser(:,:,:,ikz),G_greater(:,:,:,ikz),mu,(/temps,tempd/))
  enddo
  call write_spectrum_summed_over_kz('ldos',iter,G_retarded,nen,En,nphiz,length,NB,Lx,(/1.0,-2.0/))
  call write_spectrum_summed_over_kz('ndos',iter,G_lesser,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum_summed_over_kz('pdos',iter,G_greater,nen,En,nphiz,length,NB,Lx,(/1.0,-1.0/))
  !        
  print *, 'calc P'
  !
  nopmax=nen/2-10
  ndiag=NB*NS*2
  ! Pij^<>(hw,kz') = \int_dE Gij^<>(E,kz) * Gji^><(E-hw,kz-kz')
  ! Pij^r(hw,kz')  = \int_dE Gij^<(E,kz) * Gji^a(E-hw,kz-kz') + Gij^r(E,kz) * Gji^<(E-hw,kz-kz')
  !$omp parallel default(none) private(ndiag,l,h,iop,nop,ie,ikz,ikzd,iqz,dE,i,j) shared(nopmax,P_lesser,P_greater,P_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,nphiz)
  !$omp do
  do nop=-nopmax,nopmax
    do iqz=1,nphiz
      iop=nop+nen/2
      P_lesser(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)
      P_greater(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)    
      P_retarded(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)    
      do ie = max(nop+1,1),min(nen,nen+nop) 
        do ikz=1,nphiz
          dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )   
          do i = 1, nm_dev        
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)
              ikzd=ikz-iqz + nphiz/2
              if (ikzd<1) ikzd=ikzd+nphiz
              if (ikzd>nphiz) ikzd=ikzd-nphiz
              P_lesser(i,l:h,iop,iqz) = P_lesser(i,l:h,iop,iqz) + dE* G_lesser(i,l:h,ie,ikz) * G_greater(l:h,i,ie-nop,ikzd)
              P_greater(i,l:h,iop,iqz) = P_greater(i,l:h,iop,iqz) + dE* G_greater(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd)        
              P_retarded(i,l:h,iop,iqz) = P_retarded(i,l:h,iop,iqz) + &
                  & dE* (G_lesser(i,l:h,ie,ikz) * conjg(G_retarded(i,l:h,ie-nop,ikzd)) + G_retarded(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd))        
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
!  call write_spectrum_summed_over_kz('PR',iter,P_retarded,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('PL',iter,P_lesser  ,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('PG',iter,P_greater ,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc W'
  !
  do nop=-nopmax+nen/2,nopmax+nen/2   
    do iqz=1,nphiz
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V(:,:,iqz),nm_dev,P_retarded(:,:,nop,iqz),nm_dev,czero,B,nm_dev)   
      do i=1,nm_dev
        B(i,i) = 1.0d0 + B(i,i)
      enddo  
      ! calculate W^r = (I - V P^r)^-1 V
      call invert(B,nm_dev)
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V(:,:,iqz),nm_dev,czero,W_retarded(:,:,nop,iqz),nm_dev)     
      ! calculate W^< and W^> = W^r P^<> W^r dagger
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop,iqz),nm_dev,P_lesser(:,:,nop,iqz),nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop,iqz),nm_dev,czero,W_lesser(:,:,nop,iqz),nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop,iqz),nm_dev,P_greater(:,:,nop,iqz),nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop,iqz),nm_dev,czero,W_greater(:,:,nop,iqz),nm_dev)   
    enddo
  enddo
!  call write_spectrum_summed_over_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('WL',iter,W_lesser,  nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('WG',iter,W_greater, nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc SigGW'
  !
  ndiag=NB*NS*2
  nopmax=nen/2-10
  Sig_greater_new = dcmplx(0.0d0,0.0d0)
  Sig_lesser_new = dcmplx(0.0d0,0.0d0)
  Sig_retarded_new = dcmplx(0.0d0,0.0d0)
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
  ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
  !$omp parallel default(none) private(ndiag,l,h,nop,ie,i,j,iop,ikz,ikzd,iqz) shared(nopmax,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,dE,nphiz)
  !$omp do
  do ie=1,nen
    do ikz=1,nphiz
      do nop= -nopmax,nopmax    
        if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then      
          do iqz=1,nphiz
            iop=nop+nen/2
            ikzd=ikz-iqz + nphiz/2            
            if (ikzd<1) ikzd=ikzd+nphiz
            if (ikzd>nphiz) ikzd=ikzd-nphiz
            do i = 1,nm_dev   
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)       
              Sig_lesser_new(i,l:h,ie,ikz)=Sig_lesser_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h,iop,iqz)
              Sig_greater_new(i,l:h,ie,ikz)=Sig_greater_new(i,l:h,ie,ikz)+G_greater(i,l:h,ie-nop,ikzd)*W_greater(i,l:h,iop,iqz)
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h,iop,iqz) 
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_retarded(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h,iop,iqz) 
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_retarded(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h,iop,iqz)          
            enddo      
          enddo
        endif
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
  ! mixing with previous ones
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
  ! get leads sigma
  do iqz=1,nphiz
    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
  enddo
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    do iqz=1,nphiz
      call expand_size_bycopy(Sig_retarded(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_lesser(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_greater(:,:,ie,iqz),nm_dev,NB,2)
    enddo
  enddo
  call write_spectrum_summed_over_kz('SigR',iter,Sig_retarded,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
end do  
deallocate(siglead,B)
end subroutine green_solve_gw_2D


! driver for iterating G -> P -> W -> Sig 
subroutine green_solve_gw_1D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
implicit none
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
!----
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
integer :: iter,ie,nopmax
integer :: i,j,nm,nop,l,h,iop,ndiag
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
print *,'====== green_solve_gw_1D ======'
allocate(siglead(NB*NS,NB*NS,nen,2))
siglead=dcmplx(0.0d0,0.0d0)
allocate(B(nm_dev,nm_dev))
allocate(tot_cur(nm_dev,nm_dev))
allocate(tot_ecur(nm_dev,nm_dev))
allocate(cur(nm_dev,nm_dev,nen))
allocate(Ispec(nm_dev,nm_dev,nen))
allocate(Itot(nm_dev,nm_dev))
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
do iter=0,niter
  ! empty files for sancho 
  open(unit=101,file='sancho_gbb.dat',status='unknown')
  close(101)
  open(unit=101,file='sancho_g00.dat',status='unknown')
  close(101)
  open(unit=101,file='sancho_sig.dat',status='unknown')
  close(101)
  print *, 'calc G'  
  call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,mu,(/temps,tempd/))
 ! if (iter == 0) then     
 !   call calc_n_electron(G_lesser,G_greater,nen,En,NS,NB,nm_dev,nelec,pelec)    
 ! else    
 !   call calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,(/temps,tempd/),mu)    
 !   mu=(/ mus, mud /) - 0.2*sum(mu-(/ mus, mud /))/2.0d0 ! move Fermi level because Sig_GW shifts slightly the energies
 !   print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)    
 ! end if  
  call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
  call write_current_spectrum('Jdens',iter,cur,nen,en,length,NB,Lx)
  call write_current('I',iter,tot_cur,length,NB,NS,Lx)
  call write_current('EI',iter,tot_ecur,length,NB,NS,Lx)
  call write_spectrum('ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0,-2.0/))
  call write_spectrum('ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0,-1.0/))
  !call write_matrix_summed_overE('Gr',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
  !        
  print *, 'calc P'  
  nopmax=nen/2-10
  ndiag=NB*NS*2
  ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
  ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
  !$omp parallel default(none) private(ndiag,l,h,iop,nop,ie,dE,i,j) shared(nopmax,P_lesser,P_greater,P_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded)
  !$omp do
  do nop=-nopmax,nopmax
      iop=nop+nen/2
      P_lesser(:,:,iop) = dcmplx(0.0d0,0.0d0)
      P_greater(:,:,iop) = dcmplx(0.0d0,0.0d0)    
      P_retarded(:,:,iop) = dcmplx(0.0d0,0.0d0)    
      do ie = max(nop+1,1),min(nen,nen+nop) 
        if (ie.eq.1) then
          dE = dcmplx(0.0d0 , -1.0d0*( En(ie+1) - En(ie) ) / 2.0d0 / pi )	  
        else
          dE = dcmplx(0.0d0 , -1.0d0*( En(ie) - En(ie-1) ) / 2.0d0 / pi )	  
        endif
        do i = 1, nm_dev        
            l=max(i-ndiag,1)
            h=min(nm_dev,i+ndiag)
            P_lesser(i,l:h,iop) = P_lesser(i,l:h,iop) + dE* G_lesser(i,l:h,ie) * G_greater(l:h,i,ie-nop)
            P_greater(i,l:h,iop) = P_greater(i,l:h,iop) + dE* G_greater(i,l:h,ie) * G_lesser(l:h,i,ie-nop)        
            P_retarded(i,l:h,iop) = P_retarded(i,l:h,iop) + dE* (G_lesser(i,l:h,ie) * conjg(G_retarded(i,l:h,ie-nop)) & 
                                  & + G_retarded(i,l:h,ie) * G_lesser(l:h,i,ie-nop))        
        enddo
      enddo
  enddo
  !$omp end do
  !$omp end parallel
!  call write_spectrum('PR',iter,P_retarded,nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('PL',iter,P_lesser  ,nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('PG',iter,P_greater ,nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc W'  
  do nop=-nopmax+nen/2,nopmax+nen/2   
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded(:,:,nop),nm_dev,czero,B,nm_dev)   
    do i=1,nm_dev
      B(i,i) = 1.0d0 + B(i,i)
    enddo  
    ! calculate W^r = (I - V P^r)^-1 V
    call invert(B,nm_dev)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,W_retarded(:,:,nop),nm_dev)     
    ! calculate W^< and W^> = W^r P^<> W^r dagger
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop),nm_dev,P_lesser(:,:,nop),nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop),nm_dev,czero,W_lesser(:,:,nop),nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop),nm_dev,P_greater(:,:,nop),nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop),nm_dev,czero,W_greater(:,:,nop),nm_dev)   
  enddo
  call write_spectrum('WR',iter,W_retarded,nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('WL',iter,W_lesser,  nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('WG',iter,W_greater, nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
  call write_matrix_summed_overE('W_r',iter,W_retarded,nen,en,length,NB,(/1.0,1.0/))
  !
  print *, 'calc SigGW'
  ndiag=NB*NS*2
  nopmax=nen/2-10
  Sig_greater_new = dcmplx(0.0d0,0.0d0)
  Sig_lesser_new = dcmplx(0.0d0,0.0d0)
  Sig_retarded_new = dcmplx(0.0d0,0.0d0)
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
  ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
  !$omp parallel default(none) private(ndiag,l,h,nop,ie,i,j,iop) shared(nopmax,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,dE)
  !$omp do
  do ie=1,nen
    do nop= -nopmax,nopmax    
      if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then      
        iop=nop+nen/2
        do i = 1,nm_dev   
          l=max(i-ndiag,1)
          h=min(nm_dev,i+ndiag)       
          Sig_lesser_new(i,l:h,ie)=Sig_lesser_new(i,l:h,ie)+G_lesser(i,l:h,ie-nop)*W_lesser(i,l:h,iop)
          Sig_greater_new(i,l:h,ie)=Sig_greater_new(i,l:h,ie)+G_greater(i,l:h,ie-nop)*W_greater(i,l:h,iop)
          Sig_retarded_new(i,l:h,ie)=Sig_retarded_new(i,l:h,ie)+G_lesser(i,l:h,ie-nop)*W_retarded(i,l:h,iop) &
                                                            &  +G_retarded(i,l:h,ie-nop)*W_lesser(i,l:h,iop) &
                                                            &  +G_retarded(i,l:h,ie-nop)*W_retarded(i,l:h,iop)          
        enddo      
      endif
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
  ! mixing with the previous one
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
  ! get leads sigma
  siglead(:,:,:,1) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:)
  siglead(:,:,:,2) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:)    
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    call expand_size_bycopy(Sig_retarded(:,:,ie),nm_dev,NB,2)
    call expand_size_bycopy(Sig_lesser(:,:,ie),nm_dev,NB,2)
    call expand_size_bycopy(Sig_greater(:,:,ie),nm_dev,NB,2)
  enddo
  call write_spectrum('SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0,1.0/))
  !call write_matrix_summed_overE('Sigma_r',iter,Sig_retarded,nen,en,length,NB,(/1.0,1.0/))
  !!!! calculate collision integral
  call calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
  call write_spectrum('Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0,1.0/))
enddo                
deallocate(siglead)
deallocate(B,cur,tot_cur,tot_ecur)
deallocate(Ispec,Itot)
end subroutine green_solve_gw_1D


subroutine expand_size_bycopy(A,nm,nb,add)
complex(8),intent(inout)::A(nm,nm)
integer, intent(in)::nm,add,nb
integer::i,nm0,l,l2
nm0=nm-nb*add*2
A(1:add*nb,:)=0.0d0
A(:,1:add*nb)=0.0d0
A(add*nb+nm0+1:nm,:)=0.0d0
A(:,add*nb+nm0+1:nm)=0.0d0
do i=0,add-1
  A(i*nb+1:i*nb+nb,i*nb+1:i*nb+nm0)=A(add*nb+1:add*nb+nb,add*nb+1:add*nb+nm0)
  A(i*nb+1:i*nb+nm0,i*nb+1:i*nb+nb)=A(add*nb+1:add*nb+nm0,add*nb+1:add*nb+nb)
  l=add*nb+nm0+i*nb
  l2=add*nb+i*nb+nb
  A(l+1:l+nb,l2+1:l2+nm0)=A(add*nb+nm0-nb+1:add*nb+nm0,add*nb+1:add*nb+nm0)
  A(l2+1:l2+nm0,l+1:l+nb)=A(add*nb+1:add*nb+nm0,add*nb+nm0-nb+1:add*nb+nm0)  
enddo
end subroutine expand_size_bycopy

! calculate number of electrons and holes from G< and G> 
subroutine calc_n_electron(G_lesser,G_greater,nen,E,NS,NB,nm_dev,nelec,pelec)
complex(8), intent(in) :: G_lesser(nm_dev,nm_dev,nen)
complex(8), intent(in) :: G_greater(nm_dev,nm_dev,nen)
real(8), intent(in)    :: E(nen)
integer, intent(in)    :: NS,NB,nm_dev,nen
real(8), intent(out)   :: nelec(2),pelec(2)
integer::i,j
nelec=0.0d0
pelec=0.0d0
do i=1,nen
  do j=1,NS*NB
    nelec(1)=nelec(1)+aimag(G_lesser(j,j,i))*(E(2)-E(1))
    pelec(1)=pelec(1)-aimag(G_greater(j,j,i))*(E(2)-E(1))
  enddo
enddo
do i=1,nen
  do j=nm_dev-NS*NB+1,nm_dev
    nelec(2)=nelec(2)+aimag(G_lesser(j,j,i)*(E(2)-E(1)))
    pelec(1)=pelec(1)-aimag(G_greater(j,j,i))*(E(2)-E(1))
  enddo
enddo
end subroutine calc_n_electron

! determine the quasi-fermi level from the Gr and electron/hole number
subroutine calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,Temp,mu)
real(8),intent(in)::Temp(2)
real(8),intent(out)::mu(2)
complex(8), intent(in) :: G_retarded(nm_dev,nm_dev,nen)
real(8), intent(in)    :: En(nen)
integer, intent(in)    :: NS,NB,nm_dev,nen
real(8), intent(in)    :: nelec(2),pelec(2)
real(8)::dE,n(nen),p(nen),fd,mun,mup
real(8),allocatable::dos(:),K(:),Q(:),fermi_derivative(:,:)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
integer::nefermi,i,j,l
allocate(dos(nen))
allocate(K(nen))
allocate(Q(nen))
!nefermi=1+2*floor(10.0*boltz*maxval(temp)/(En(2)-En(1))) ! number of energy points for the fermi derivative function
nefermi=nen
dE=En(2)-En(1)
do i=1,2
  dos=0.0d0
  K=0.0d0
  Q=0.0d0
  n=0.0d0
  p=0.0d0
  do j=1,nen
    if (i==1) then
      do l=1,NS*NB
        dos(j) = dos(j)+aimag(G_retarded(l,l,j))
      enddo
    else
      do l=nm_dev-NS*NB+1,nm_dev
        dos(j) = dos(j)+aimag(G_retarded(l,l,j))
      enddo
    endif    
    dos(j)=-2.0d0*dos(j)
    if (j>1) K(j)=K(j-1)+dos(j)*dE    
    if (j>1) Q(nen-j+1)=Q(nen-j+2)+dos(j)*dE    
  enddo
  ! search for the Fermi level
!  if (dE<(BOLTZ*TEMP(i))) then
!    allocate(fermi_derivative(nefermi,2))
!    do j=1,nefermi
!      fd=ferm(dble(dE)*dble(j-nefermi/2-1)/(BOLTZ*TEMP(i)))
!      fermi_derivative(j,i) = -fd*(1.0d0-fd)/(BOLTZ*TEMP(i))    
!    enddo
!    do j=1,nen
!      n(j)=-sum(K(max(j-nefermi/2+1,1):min(j+nefermi/2,nen))*fermi_derivative(max(nefermi/2-j,1):min(nen-j,nefermi),i))*dE
!    enddo
!    deallocate(fermi_derivative)
!  else ! energy grid too coarse
!    ! approximate F-D to step-function
    n = K
    p = Q
!  endif
  n=n-nelec(i)  
  p=p-pelec(i)
  do j=2,nen
    if ((n(j)>=0.0).and.(n(j-1)<=0.0)) then
      mun=En(j)
      exit
    endif
  enddo
  do j=nen-1,1,-1
    if ((p(j)>=0.0).and.(p(j+1)<=0.0)) then
      mup=En(j)
      exit
    endif
  enddo
  mu(i)=(mun+mup)/2.0d0
enddo
deallocate(dos,K)
end subroutine calc_fermi_level


! calculate Gr and optionally G<>
subroutine green_calc_g(ne,E,num_lead,nm_dev,nm_lead,max_nm_lead,Ham,H00,H10,Siglead,T,Scat_Sig_retarded,Scat_Sig_lesser,Scat_Sig_greater,G_retarded,G_lesser,G_greater,mu,temp,mode)
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
complex(8), intent(in) :: Siglead(max_nm_lead,max_nm_lead,ne,num_lead) ! lead sigma_r scattering
complex(8), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
complex(8), intent(in) :: Scat_Sig_retarded(nm_dev,nm_dev,ne) ! scattering Selfenergy
complex(8), intent(in) :: Scat_Sig_lesser(nm_dev,nm_dev,ne)
complex(8), intent(in) :: Scat_Sig_greater(nm_dev,nm_dev,ne)
complex(8), intent(inout) :: G_retarded(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_lesser(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_greater(nm_dev,nm_dev,ne)
real(8), intent(in), optional :: mu(num_lead), temp(num_lead)
character(len=*),intent(in), optional :: mode
integer :: i,j,nm,ie
complex(8), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
real(8) :: fd
logical :: solve_Gr
allocate(sig(nm_dev,nm_dev))  
solve_Gr = .true.
if (present(mode).and.(mode=='use_gr')) then
  solve_Gr = .false.
endif
do ie = 1, ne
  if (solve_Gr) G_retarded(:,:,ie) = - Ham(:,:) - Scat_Sig_retarded(:,:,ie) 
  if ((present(G_lesser)).or.(present(G_greater))) then    
    if (ie .eq. 1) then
      allocate(sig_lesser(nm_dev,nm_dev))
      allocate(sig_greater(nm_dev,nm_dev))
      allocate(B(nm_dev,nm_dev))
      allocate(C(nm_dev,nm_dev))
    end if
    sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
    sig_greater(:,:) = dcmplx(0.0d0,0.0d0)      
  end if    
  ! compute and add contact self-energies    
  open(unit=101,file='sancho_gbb.dat',status='unknown',position='append')
  open(unit=102,file='sancho_g00.dat',status='unknown',position='append')
  open(unit=103,file='sancho_sig.dat',status='unknown',position='append')
  do i = 1,num_lead
    NM = nm_lead(i)    
    allocate(S00(nm,nm))
    allocate(G00(nm,nm))
    allocate(GBB(nm,nm))
    allocate(A(nm_dev,nm))    
    call identity(S00,nm)
    call sancho(NM,E(ie),S00,H00(1:nm,1:nm,i)+siglead(1:nm,1:nm,ie,i),H10(1:nm,1:nm,i),G00,GBB)
    call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
    write(101,'(i4,2E15.4)') i, E(ie), -aimag(trace(GBB,nm))*2.0d0
    write(102,'(i4,2E15.4)') i, E(ie), -aimag(trace(G00,nm))*2.0d0
    write(103,'(i4,2E15.4)') i, E(ie), -aimag(trace(sig,nm_dev))*2.0d0
    if (solve_Gr) G_retarded(:,:,ie) = G_retarded(:,:,ie) - sig(:,:)
    if ((present(G_lesser)).or.(present(G_greater))) then      
      fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))		
      B(:,:) = conjg(sig(:,:))
      C(:,:) = transpose(B(:,:))
      B(:,:) = sig(:,:) - C(:,:)
      sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd	        
      sig_greater(:,:) = sig_greater(:,:) + B(:,:)*(1.0d0-fd)	        
    end if
    deallocate(S00,G00,GBB,A)
  end do  
  close(101)
  close(102)
  close(103)
  if (solve_Gr) then
    do i = 1,nm_dev
      G_retarded(i,i,ie) = G_retarded(i,i,ie) + dcmplx(E(ie),0.0d0)
    end do
    !
    call invert(G_retarded(:,:,ie),nm_dev) 
    !
  endif
  if ((present(G_lesser)).or.(present(G_greater))) then    
    sig_lesser = sig_lesser + Scat_Sig_lesser(:,:,ie)
    sig_greater = sig_greater + Scat_Sig_greater(:,:,ie)     
    if (present(G_lesser)) then
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
      G_lesser(:,:,ie) = C
    end if
    if (present(G_greater)) then
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_greater,nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
      G_greater(:,:,ie) = C      
    end if      
  end if 
end do  
deallocate(sig)
if ((present(G_lesser)).or.(present(G_greater))) then      
  deallocate(B,C,sig_lesser,sig_greater)
end if
if (solve_Gr) G_retarded(:,:,:)=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
end subroutine green_calc_g

! calculate bond current using I_ij = H_ij G<_ji - H_ji G^<_ij
subroutine calc_bond_current(H,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
implicit none
complex(8),intent(in)::H(nm_dev,nm_dev),G_lesser(nm_dev,nm_dev,nen)
real(8),intent(in)::en(nen),spindeg
integer,intent(in)::nen,nm_dev
real(8),intent(out)::tot_cur(nm_dev,nm_dev) ! total bond current density
real(8),intent(out),optional::tot_ecur(nm_dev,nm_dev) ! total bond energy current density
real(8),intent(out),optional::cur(nm_dev,nm_dev,nen) ! energy resolved bond current density
!----
complex(8),allocatable::B(:,:)
integer::ie,io,jo
real(8),parameter::tpi=6.28318530718  
  allocate(B(nm_dev,nm_dev))
  tot_cur=0.0d0  
  tot_ecur=0.0d0
  do ie=1,nen
    do io=1,nm_dev
      do jo=1,nm_dev
        B(io,jo)=H(io,jo)*G_lesser(jo,io,ie) - H(jo,io)*G_lesser(io,jo,ie)
      enddo
    enddo    
    B=B*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
    if (present(cur)) cur(:,:,ie) = dble(B)
    if (present(tot_ecur)) tot_ecur=tot_ecur+ en(ie)*dble(B)
    tot_cur=tot_cur+ dble(B)          
  enddo
  deallocate(B)
end subroutine calc_bond_current

! calculate scattering collision integral from the self-energy
! I = Sig> G^< - Sig< G^>
subroutine calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,spindeg,nm_dev,I,Ispec)
implicit none
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_greater,G_lesser,Sig_lesser,Sig_greater
real(8),intent(in)::en(nen),spindeg
integer,intent(in)::nen,nm_dev
complex(8),intent(out)::I(nm_dev,nm_dev) ! collision integral
complex(8),intent(out),optional::Ispec(nm_dev,nm_dev,nen) ! collision integral spectrum
!----
complex(8),allocatable::B(:,:)
integer::ie
real(8),parameter::tpi=6.28318530718  
  allocate(B(nm_dev,nm_dev))
  I=dcmplx(0.0d0,0.0d0)
  do ie=1,nen
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,Sig_greater(:,:,ie),nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Sig_lesser(:,:,ie),nm_dev,G_greater(:,:,ie),nm_dev,cone,B,nm_dev) 
    I(:,:)=I(:,:)+B(:,:)
    if (present(Ispec)) Ispec(:,:,ie)=B(:,:)
  enddo
  deallocate(B)
end subroutine calc_collision


! write current into file 
subroutine write_current(dataset,i,cur,length,NB,NS,Lx)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: cur(:,:)
integer, intent(in)::i,length,NB,NS
real(8), intent(in)::Lx
integer:: j,ib,jb,ii
real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ii = 1,length-1
    tr=0.0d0          
    do ib=1,nb  
      do jb=1,nb       
        do j=ii,min(ii+NS-1,length-1)
          tr = tr+ cur((ii-1)*nb+ib,j*nb+jb)
        enddo
      enddo                        
    end do
    write(11,'(2E18.4)') dble(ii)*Lx, tr
  end do
end subroutine write_current

! write current spectrum into file (pm3d map)
subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: cur(:,:,:)
integer, intent(in)::i,nen,length,NB
real(8), intent(in)::Lx,en(nen)
integer:: ie,j,ib,jb
real(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length-1
        tr=0.0d0          
        do ib=1,nb  
          do jb=1,nb        
            tr = tr+ cur((j-1)*nb+ib,j*nb+jb,ie)
          enddo                        
        end do
        write(11,'(3E18.4)') dble(j)*Lx, en(ie), tr
    end do
    write(11,*)    
end do
close(11)
end subroutine write_current_spectrum

! write spectrum into file (pm3d map)
subroutine write_spectrum(dataset,i,G,nen,en,length,NB,Lx,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:)
integer, intent(in)::i,nen,length,NB
real(8), intent(in)::Lx,en(nen),coeff(2)
integer:: ie,j,ib
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length
        tr=0.0d0          
        do ib=1,nb
            tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie)            
        end do
        write(11,'(4E18.4)') j*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    end do
    write(11,*)    
end do
close(11)
end subroutine write_spectrum


! write a matrix summed over energy index into a file
subroutine write_matrix_summed_overE(dataset,i,G,nen,en,length,NB,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:)
integer, intent(in)::i,nen,length,NB
real(8), intent(in)::en(nen),coeff(2)
integer:: ie,j,ib,l
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do l=1,length*NB
  do j = 1,length*NB
    tr=0.0d0          
    do ie=1,nen  
        tr = tr+ G(l,j,ie)            
    end do
    tr=tr/dble(nen)
    write(11,'(2I8,2E18.4)') l,j, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
  end do
  write(11,*)    
end do
close(11)
end subroutine write_matrix_summed_overE



! write spectrum into file (pm3d map)
subroutine write_spectrum_summed_over_kz(dataset,i,G,nen,en,nkz,length,NB,Lx,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:,:)
integer, intent(in)::i,nen,length,NB,nkz
real(8), intent(in)::Lx,en(nen),coeff(2)
integer:: ie,j,ib,ikz
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length
        tr=0.0d0          
        do ib=1,nb
          do ikz=1,nkz
            tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz)            
          enddo
        enddo
        tr=tr/dble(nkz)
        write(11,'(4E18.4)') j*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    end do
    write(11,*)    
enddo
close(11)
end subroutine write_spectrum_summed_over_kz


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)
implicit none
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  integer i,j,k,nm,nmax
  COMPLEX(8) :: z
  real(8) :: E,error
  REAL(8) :: TOL=1.0D-100  ! [eV]
  COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
  COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)
  COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:)
  COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
  !COMPLEX(8), ALLOCATABLE :: WORK(:)
  !COMPLEX(8), EXTERNAL :: ZLANGE
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
    !
    call invert(A,nm)
    !
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
    tmp=H_SS
    IF ( abs(error) < TOL ) THEN	
      EXIT
    ELSE
    END IF
    IF (i .EQ. nmax) THEN
      write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error
      call abort
    END IF
  enddo
  G00 = z*S00 - H_SS
  !
  call invert(G00,nm)
  !
  GBB = z*S00 - H_BB
  !
  call invert(GBB,nm)
  !
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


subroutine subspace_invert(nm,A,add_nm,in_method)
integer, intent(in) :: nm ! dimension of subspace 
integer, intent(in) :: add_nm ! additional dimension 
complex(8), intent(inout) :: A(nm,nm)
character(len=*), intent(in),optional :: in_method
character(len=10) :: method
complex(8), dimension(:,:), allocatable :: invA,V00,V10,sigmal,sigmar,S00,G00,Gbb,B
integer::nm_lead,nm_dev
if (present(in_method)) then
  method=in_method
else
  method='direct'
endif

select case (trim(method))
  case('direct')    
    allocate(invA(nm+2*add_nm,nm+2*add_nm))
    invA=cmplx(0.0d0,0.0d0)
    invA(add_nm+1:add_nm+nm,add_nm+1:add_nm+nm)=A    
    invA(1:add_nm,add_nm+1:2*add_nm)=A(1:add_nm,add_nm+1:2*add_nm)
    invA(add_nm+1:2*add_nm,1:add_nm)=A(add_nm+1:2*add_nm,1:add_nm)    
    invA(add_nm+nm+1:2*add_nm+nm,nm+1:add_nm+nm)=A(nm-add_nm+1:nm,nm-2*add_nm+1:nm-add_nm)
    invA(nm+1:add_nm+nm,add_nm+nm+1:2*add_nm+nm)=A(nm-2*add_nm+1:nm-add_nm,nm-add_nm+1:nm)
    call invert(invA,nm+2*add_nm)
    A=invA(add_nm+1:add_nm+nm,add_nm+1:add_nm+nm)
    deallocate(invA)
  case('sancho')       
    allocate(V00(add_nm,add_nm))
    allocate(V10(add_nm,add_nm))
    allocate(S00(add_nm,add_nm))
    allocate(G00(add_nm,add_nm))
    allocate(GBB(add_nm,add_nm))
    allocate(B(add_nm,add_nm))
    allocate(sigmal(add_nm,add_nm))
    allocate(sigmar(add_nm,add_nm))
    nm_lead=add_nm
    nm_dev=nm
    ! get OBC on left  
    V00 = - A(1:nm_lead,1:nm_lead) 
    V10 = - A(nm_lead+1:2*nm_lead,1:nm_lead)
    call identity(S00,nm_lead)
    call sancho(nm_lead,0.0d0,S00,V00,V10,G00,GBB)
    call zgemm('n','n',nm_lead,nm_lead,nm_lead,cone,V10,nm_lead,G00,nm_lead,czero,B,nm_lead) 
    call zgemm('n','c',nm_lead,nm_lead,nm_lead,cone,B,nm_lead,V10,nm_lead,czero,sigmal,nm_lead)  
    ! get OBC on right
    call sancho(nm_lead,0.0d0,S00,V00,transpose(conjg(V10)),G00,GBB)
    call zgemm('c','n',nm_lead,nm_lead,nm_lead,cone,V10,nm_lead,G00,nm_lead,czero,B,nm_lead) 
    call zgemm('n','n',nm_lead,nm_lead,nm_lead,cone,B,nm_lead,V10,nm_lead,czero,sigmar,nm_lead)  
    !    
    A(1:nm_lead,1:nm_lead) = A(1:nm_lead,1:nm_lead) + sigmal
    A(nm_dev-nm_lead+1:nm_dev,nm_dev-nm_lead+1:nm_dev) = A(nm_dev-nm_lead+1:nm_dev,nm_dev-nm_lead+1:nm_dev) + sigmar
    !
    call invert(A,nm_dev)
    deallocate(V00,V10,sigmar,sigmal,G00,Gbb,B,S00)
  
end select  

end subroutine subspace_invert


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

Function trace(A,nn)
  implicit none      
  integer :: nn,i        
  complex(8), dimension(nn,nn),intent(in) :: A
  complex(8) :: trace, tr
  tr=dcmplx(0.0d0,0.0d0)
  do i=1,nn
    tr=tr+A(i,i)
  enddo
  trace=tr
end function trace

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


FUNCTION STRING(inn)
  IMPLICIT NONE
  INTEGER, PARAMETER :: POS= 4
  INTEGER, INTENT(IN) :: inn
  CHARACTER(LEN=POS) :: STRING
  !............................................................
  INTEGER :: cifra, np, mm, num  
  IF (inn > (10**POS)-1) stop "ERRORE: (inn > (10**3)-1)  in STRING"
  num= inn
  DO np= 1, POS
     mm= pos-np
     cifra= num/(10**mm)            
     STRING(np:np)= ACHAR(48+cifra)
     num= num - cifra*(10**mm)
  END DO
END FUNCTION STRING

end module green
