!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf
use green,only:get_OBC_blocks_for_W,get_dL_OBC_for_W    
use mpi

implicit none 

private

public :: green_rgf_cms, green_rgf_solve_gw_1d
public :: green_rgf_solve_gw_3d
public :: green_rgf_solve_gw_ephoton_3d
public :: green_rgf_solve_gw_ephoton_3d_mpi

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
real(8), parameter :: m0=9.109d-31 ! kg
real(8), parameter :: eps0=8.854d-12 ! C/V/m 
real(8), parameter :: c0=2.998d8 ! m/s
real(8), parameter :: e0=1.6022d-19 ! C
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: tpi = 3.14159265359d0*2.0d0

CONTAINS


! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic(nm,length,nen,En,nop,Mii,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,pre_fact,intensity,hw)
integer,intent(in)::nm,length,nen,nop
real(8),intent(in)::en(nen),pre_fact,intensity,hw
complex(8),intent(in),dimension(nm,nm,length)::Mii ! e-photon interaction matrix blocks
complex(8),intent(in),dimension(nm,nm,length,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm,nm,length,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie,ix
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
  Sig_lesser=czero
  Sig_greater=czero
  Sig_retarded=czero
  ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
  !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
  !$omp parallel default(none) private(ie,A,B,ix) shared(length,nop,nen,nm,G_lesser,G_greater,Sig_lesser,Sig_greater,Mii)
  allocate(B(nm,nm))
  allocate(A(nm,nm))  
  !$omp do
  do ie=1,nen
    do ix=1,length
      ! Sig^<(E)
      A = czero
      if (ie-nop>=1) A =A+ G_lesser(:,:,ix,ie-nop)
      if (ie+nop<=nen) A =A+ G_lesser(:,:,ix,ie+nop)
      call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
      call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
      Sig_lesser(:,:,ix,ie) = A 
      ! Sig^>(E)
      A = czero
      if (ie-nop>=1) A =A+ G_greater(:,:,ix,ie-nop)
      if (ie+nop<=nen) A =A+ G_greater(:,:,ix,ie+nop)
      call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
      call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
      Sig_greater(:,:,ix,ie) = A
    enddo
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Sig_greater=Sig_greater*pre_fact*intensity/hw**2
  Sig_lesser=Sig_lesser*pre_fact*intensity/hw**2
  Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
end subroutine calc_sigma_ephoton_monochromatic


subroutine green_rgf_solve_gw_ephoton_3d(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs)
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
  complex(8),allocatable,dimension(:,:,:)::Sii
  complex(8),allocatable,dimension(:,:,:,:)::Mii,M1i
  real(8)::tr(nen,nky*nkz),tre(nen,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk
  complex(8)::dE, B(nm,nm)
  !
  print *,'======== green_rgf_solve_gw_ephoton_3D ========'
  nk=nky*nkz
  allocate(g_r(nm,nm,nx,nen,nk))
  allocate(g_r_i1(nm,nm,nx,nen,nk))
  allocate(g_greater(nm,nm,nx,nen,nk))
  allocate(g_lesser(nm,nm,nx,nen,nk))
  allocate(cur(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_new(nm,nm,nx,nen,nk))
  allocate(sigma_greater_new(nm,nm,nx,nen,nk))
  allocate(sigma_r_new(nm,nm,nx,nen,nk))
  allocate(P_lesser(nm,nm,nx,nen,nk))
  allocate(P_greater(nm,nm,nx,nen,nk))
  allocate(P_retarded(nm,nm,nx,nen,nk))
!  allocate(P_lesser_1i(nm,nm,nx,nen,nk))
!  allocate(P_greater_1i(nm,nm,nx,nen,nk))
!  allocate(P_retarded_1i(nm,nm,nx,nen,nk))
  allocate(W_lesser(nm,nm,nx,nen,nk))
  allocate(W_greater(nm,nm,nx,nen,nk))
  allocate(W_retarded(nm,nm,nx,nen,nk))
!  allocate(W_lesser_i1(nm,nm,nx,nen,nk))
!  allocate(W_greater_i1(nm,nm,nx,nen,nk))
!  allocate(W_retarded_i1(nm,nm,nx,nen,nk))  
  allocate(Mii(nm,nm,nx,nk))
  Mii(:,:,:,:)=czero
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo   
  print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  nopphot=floor(hw / (En(2)-En(1)))
  print *,'nop photon=',nopphot
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'      
    do ik=1,nk
      print *, ' ik=', ik,'/',nk
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nx,en,mu,Hii,H1i,sigma_lesser_gw,sigma_greater_gw,sigma_r_gw,g_lesser,g_greater,g_r,tre,tr,cur,g_r_i1)    
      !$omp do 
      do ie=1,nen     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik),H1i(:,:,:,ik),sigma_lesser_gw(:,:,:,ie,ik),sigma_greater_gw(:,:,:,ie,ik),&
                        & sigma_r_gw(:,:,:,ie,ik),g_lesser(:,:,:,ie,ik),g_greater(:,:,:,ie,ik),g_r(:,:,:,ie,ik),tr(ie,ik),&
                        & tre(ie,ik),cur(:,:,:,ie,ik),g_r_i1(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    !
    call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum_summed_over_kz('gw_Jdens',iter,cur,nen,En,nx,NB*NS,Lx*NS,nk)  
    call write_transmission_spectrum_k('gw_trR',iter,tr*spindeg,nen,En,nk)
    call write_transmission_spectrum_k('gw_trL',iter,tre*spindeg,nen,En,nk)
    open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10      
    P_lesser = dcmplx(0.0d0,0.0d0)
    P_greater = dcmplx(0.0d0,0.0d0)    
    P_retarded = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i = dcmplx(0.0d0,0.0d0)
!    P_greater_1i = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nx
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
              do iky=1,nky
                do ikz=1,nkz              
                  ik=ikz + (iky-1)*nkz
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky                
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1
                  ikd=ikzd + (ikyd-1)*nkz
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser(i,l:h,ix,iop,iq) = P_lesser(i,l:h,ix,iop,iq) + G_lesser(i,l:h,ix,ie,ik) * G_greater(l:h,i,ix,ie-nop,ikd)
                    P_greater(i,l:h,ix,iop,iq) = P_greater(i,l:h,ix,iop,iq) + G_greater(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd)        
                    P_retarded(i,l:h,ix,iop,iq) = P_retarded(i,l:h,ix,iop,iq) + (G_lesser(i,l:h,ix,ie,ik) * conjg(G_r(i,l:h,ix,ie-nop,ikd)) & 
                                              & + G_r(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo      
    enddo          
    !$omp end do 
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum_summed_over_k('gw_PR',iter,P_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'                
    do iq=1,nky*nkz
      print *, ' iq=', iq,'/',nk
      !$omp parallel default(none) private(nop) shared(iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,w_lesser,w_greater,w_retarded)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    call write_spectrum_summed_over_k('gw_WR',iter,W_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nx
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new(i,l:h,ix,ie,ik)=Sigma_lesser_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq)
                      Sigma_greater_new(i,l:h,ix,ie,ik)=Sigma_greater_new(i,l:h,ix,ie,ik)+G_greater(i,l:h,ix,ie-nop,ikd)*W_greater(i,l:h,ix,iop,iq)
                      Sigma_r_new(i,l:h,ix,ie,ik)=Sigma_r_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! symmetrize the selfenergies
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new(:,:,ix,ie,ik))
          Sigma_r_new(:,:,ix,ie,ik) = (Sigma_r_new(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new(:,:,ix,ie,ik))
          Sigma_lesser_new(:,:,ix,ie,ik) = (Sigma_lesser_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new(:,:,ix,ie,ik))
          Sigma_greater_new(:,:,ix,ie,ik) = (Sigma_greater_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! 
    call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    if (iter>=(niter-5)) then
      print *, 'calc SigEPhoton'
      do ik=1,nk
        print *, ' ik=', ik,'/',nk
        call calc_sigma_ephoton_monochromatic(nm,nx,nen,En,nopphot,Mii(:,:,:,ik),G_lesser(:,:,:,:,ik),G_greater(:,:,:,:,ik),Sigma_r_new(:,:,:,:,ik),Sigma_lesser_new(:,:,:,:,ik),Sigma_greater_new(:,:,:,:,ik),pre_fact,intensity,hw)
      enddo
      Sigma_r_gw       = Sigma_r_gw + Sigma_r_new 
      Sigma_lesser_gw  = Sigma_lesser_gw + Sigma_lesser_new 
      Sigma_greater_gw = Sigma_greater_gw + Sigma_greater_new  
      call write_spectrum_summed_over_k('eph_SigR',iter,Sigma_r_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('eph_SigL',iter,Sigma_lesser_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('eph_SigG',iter,Sigma_greater_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    endif
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:,:)=Sigma_r_gw(:,:,3,:,:)
      Sigma_lesser_gw(:,:,ix,:,:)=Sigma_lesser_gw(:,:,3,:,:)
      Sigma_greater_gw(:,:,ix,:,:)=Sigma_greater_gw(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:,:)=Sigma_r_gw(:,:,nx-2,:,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:,:)=Sigma_lesser_gw(:,:,nx-2,:,:)
      Sigma_greater_gw(:,:,nx-ix+1,:,:)=Sigma_greater_gw(:,:,nx-2,:,:)
    enddo    
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
  deallocate(Mii)
end subroutine green_rgf_solve_gw_ephoton_3d



subroutine green_rgf_solve_gw_ephoton_3d_mpi(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs)
use, intrinsic :: iso_fortran_env
use mpi_f08
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser,cur
  complex(8),allocatable,dimension(:,:,:,:,:)::P_r,P_greater,P_lesser
  complex(8),allocatable,dimension(:,:,:,:,:)::W_r,W_greater,W_lesser
  complex(8),allocatable,dimension(:,:,:,:,:)::Sigma_r_gw,Sigma_lesser_gw,Sigma_greater_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r_local_k,g_greater_local_k,g_lesser_local_k,cur_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r_local_x,g_greater_local_x,g_lesser_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw_local_x,sigma_greater_gw_local_x,sigma_r_gw_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new_local_k,sigma_greater_new_local_k,sigma_r_new_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new_local_x,sigma_greater_new_local_x,sigma_r_new_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_local_k,P_greater_local_k,P_retarded_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_local_x,P_greater_local_x,P_retarded_local_x  
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_local_k,W_greater_local_k,W_retarded_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_local_x,W_greater_local_x,W_retarded_local_x  
  complex(8),allocatable,dimension(:,:,:,:,:)::sbuf,rbuf
  complex(8),dimension(nm,nm,nx,nky*nkz)::Mii!,M1i
  real(8),allocatable,dimension(:,:)::tr,tre,tr_local_k,tre_local_k
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk,ndata
  complex(8)::dE, B(nm,nm)
  ! --- mpi
  integer(kind=int32) :: rank,num_proc
  integer(kind=int32) :: ierror
  type(MPI_Status)   :: mpistatus
  integer::scount,rcount
  integer::nnode,inode,nnx
  ! Get the individual process (rank)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)  
  !
  call MPI_Barrier(MPI_COMM_WORLD)
  if (rank == 0) then
    print *,'======== green_rgf_solve_gw_ephoton_3D_MPI ========'
  endif
  !
  nnode=num_proc
  nk=nky*nkz/nnode  ! each node does nk k point
  nnx=nx/nnode ! each node does nnx x points
  !  
  !
  allocate(g_r_local_k(nm,nm,nx,nen,nk))  
  allocate(g_greater_local_k(nm,nm,nx,nen,nk))
  allocate(g_lesser_local_k(nm,nm,nx,nen,nk))      
  allocate(cur_local_k(nm,nm,nx,nen,nk))
  allocate(tr_local_k(nen,nk))  
  allocate(tre_local_k(nen,nk))  
  !
  allocate(sigma_lesser_gw_local_k(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw_local_k(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw_local_k(nm,nm,nx,nen,nk)) 
  sigma_lesser_gw_local_k = czero
  sigma_greater_gw_local_k = czero
  sigma_r_gw_local_k = czero
  !
  Mii(:,:,:,:)=czero  ! e-photon coupling matrix
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo   
  call MPI_Barrier(MPI_COMM_WORLD)
  if (rank == 0) then
    print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  endif
  nopphot=floor(hw / (En(2)-En(1)))
  if (rank == 0) then
    print *,'nop photon=',nopphot
  endif
  do iter=0,niter
    if (rank == 0) then
      print *,'+ iter=',iter
      !
      print *, 'calc G'      
    endif
    !
    do ik=1,nk
      call MPI_Barrier(MPI_COMM_WORLD)
      !print *, ' ik=', ik+(rank)*nk,'/',nky*nkz      
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nk,nx,rank,en,mu,Hii,H1i,sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k,g_lesser_local_k,g_greater_local_k,g_r_local_k,tre_local_k,tr_local_k,cur_local_k)    
      !$omp do 
      do ie=1,nen             
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik+(rank)*nk),H1i(:,:,:,ik+(rank)*nk),&
                        & sigma_lesser_gw_local_k(:,:,:,ie,ik),sigma_greater_gw_local_k(:,:,:,ie,ik),&
                        & sigma_r_gw_local_k(:,:,:,ie,ik),g_lesser_local_k(:,:,:,ie,ik),g_greater_local_k(:,:,:,ie,ik),&
                        & g_r_local_k(:,:,:,ie,ik),tr_local_k(ie,ik),&
                        & tre_local_k(ie,ik),cur_local_k(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel      
    enddo      
    call write_spectrum_summed_over_k('gw_ldos_k'//string(rank),iter,g_r_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum_summed_over_k('gw_ndos_k'//string(rank),iter,g_lesser_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum_summed_over_k('gw_pdos_k'//string(rank),iter,g_greater_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum_summed_over_kz('gw_Jdens_k'//string(rank),iter,cur_local_k,nen,En,nx,NB*NS,Lx*NS,nk)  
    call write_transmission_spectrum_k('gw_trR_k'//string(rank),iter,tr_local_k*spindeg,nen,En,nk)
    call write_transmission_spectrum_k('gw_trL_k'//string(rank),iter,tre_local_k*spindeg,nen,En,nk)
    !!! reduce to sum up k points
    allocate(g_r(nm,nm,nx,nen,nk))  
    allocate(g_greater(nm,nm,nx,nen,nk))
    allocate(g_lesser(nm,nm,nx,nen,nk))  
    allocate(cur(nm,nm,nx,nen,nk))
    allocate(tr(nen,nk))  
    allocate(tre(nen,nk))  
    g_r=czero
    g_lesser=czero
    g_greater=czero
    cur=czero
    tr=0.0d0
    tre=0.0d0
    if (rank == 0) then
      print *,' reduce'
    endif
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(cur_local_k, cur, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    cur=cur/dble(nnode)        
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(g_r_local_k, g_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_r=g_r/dble(nnode)        
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(g_lesser_local_k, g_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_lesser=g_lesser/dble(nnode)               
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(g_greater_local_k, g_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_greater=g_greater/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(tr_local_k, tr, nen*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tr=tr/dble(nnode)       
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(tre_local_k, tre, nen*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tre=tre/dble(nnode)       
    if (rank == 0) then
      call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
      call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
      call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
      call write_current_spectrum_summed_over_kz('gw_Jdens',iter,cur,nen,En,nx,NB*NS,Lx*NS,nk)  
      call write_transmission_spectrum_k('gw_trR',iter,tr*spindeg,nen,En,nk)
      call write_transmission_spectrum_k('gw_trL',iter,tre*spindeg,nen,En,nk)
      open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
      write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
      close(101)
    endif
    deallocate(g_r,g_greater,g_lesser,cur,tr,tre)
    !!!
    g_r_local_k = dcmplx( 0.0d0*dble(g_r_local_k), aimag(g_r_local_k))
    g_lesser_local_k = dcmplx( 0.0d0*dble(g_lesser_local_k), aimag(g_lesser_local_k))
    g_greater_local_k = dcmplx( 0.0d0*dble(g_greater_local_k), aimag(g_greater_local_k))
    !        
    allocate(g_r_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)  
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(g_greater_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(g_lesser_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    call MPI_Barrier(MPI_COMM_WORLD)
    !!! all-to-all communication to make G local in x
    !!! G(:,:,ix,:,ik=rank) --> buff(:,:,:) --> G(:,:,ix=rank,:,ik)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    allocate(sbuf(nm,nm,nen,nx,nk),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx),stat=ierror)    
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    !!!
    scount=nm*nm*nen*nnx*nk
    rcount=nm*nm*nen*nnx*nk
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_r_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD)    
    g_r_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !!
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_lesser_local_k(:,:,ix,:,ik)
      enddo    
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD)    
    g_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)        
    !!
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_greater_local_k(:,:,ix,:,ik)
      enddo    
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD)
    g_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)       
    call MPI_Barrier(MPI_COMM_WORLD) 
    deallocate(sbuf)
    deallocate(rbuf)
    !!!
    if (rank == 0) then     
      print *, 'calc P'    
    endif
    nopmax=nen/2-10              
    allocate(P_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(P_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(P_retarded_local_x(nm,nm,nen,nky*nkz,nnx))
    P_lesser_local_x = dcmplx(0.0d0,0.0d0)
    P_greater_local_x = dcmplx(0.0d0,0.0d0)    
    P_retarded_local_x = dcmplx(0.0d0,0.0d0)        
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)           
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iky,ikyd,ik,ikd,iqz,iqy,iq) shared(ndiag,nopmax,P_lesser_local_x,P_greater_local_x,P_retarded_local_x,nen,En,nm,G_lesser_local_x,G_greater_local_x,G_r_local_x,nx,nkz,nky,nk,nnx)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nnx  ! x parallelized by MPI
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
              do iky=1,nky
                do ikz=1,nkz              
                  ik=ikz + (iky-1)*nkz
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky                
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1
                  ikd=ikzd + (ikyd-1)*nkz                
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser_local_x(i,l:h,iop,iq,ix) = P_lesser_local_x(i,l:h,iop,iq,ix)  & 
                          & + G_lesser_local_x(i,l:h,ie,ik,ix) * G_greater_local_x(l:h,i,ie-nop,ikd,ix)
                    P_greater_local_x(i,l:h,iop,iq,ix) = P_greater_local_x(i,l:h,iop,iq,ix) &
                          & + G_greater_local_x(i,l:h,ie,ik,ix) * G_lesser_local_x(l:h,i,ie-nop,ikd,ix)        
                    P_retarded_local_x(i,l:h,iop,iq,ix) = P_retarded_local_x(i,l:h,iop,iq,ix) &
                          & + (G_lesser_local_x(i,l:h,ie,ik,ix) * conjg(G_r_local_x(i,l:h,ie-nop,ikd,ix)) & 
                          & + G_r_local_x(i,l:h,ie,ik,ix) * G_lesser_local_x(l:h,i,ie-nop,ikd,ix))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo                    
    !$omp end do 
    !$omp end parallel             
    deallocate(G_r_local_x,G_lesser_local_x,G_greater_local_x)
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nky*nkz)
    P_lesser_local_x  =P_lesser_local_x*dE
    P_greater_local_x =P_greater_local_x*dE
    P_retarded_local_x=P_retarded_local_x*dE
    !
    allocate(P_lesser_local_k(nm,nm,nx,nen,nk))
    allocate(P_greater_local_k(nm,nm,nx,nen,nk))
    allocate(P_retarded_local_k(nm,nm,nx,nen,nk))
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make P local in k
    !!! P(:,:,ix=rank,:,ik) --> buff(:,:,:) --> P(:,:,ix,:,ik=rank)
    allocate(sbuf(nm,nm,nen,nky*nkz,nnx))        
    allocate(rbuf(nm,nm,nen,nx,nk))  
    !!!    
    sbuf(:,:,:,:,:)=P_retarded_local_x(:,:,:,:,:)      
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    do ik=1,nk
      do ix=1,nx    
        P_retarded_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !    
    sbuf(:,:,:,:,:)=P_lesser_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    do ik=1,nk
      do ix=1,nx
        P_lesser_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=P_greater_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    do ik=1,nk
      do ix=1,nx
        P_greater_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)
    deallocate(sbuf,rbuf)
    !!!    
    deallocate(P_retarded_local_x,P_lesser_local_x,P_greater_local_x)  
    allocate(P_r(nm,nm,nx,nen,nk))  
    allocate(P_greater(nm,nm,nx,nen,nk))
    allocate(P_lesser(nm,nm,nx,nen,nk))
    P_r=czero
    P_lesser=czero
    P_greater=czero    
    if (rank == 0) then
      print *,' reduce'
    endif
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(P_retarded_local_k, P_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_r=P_r/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(P_lesser_local_k, P_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_lesser=P_lesser/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(P_greater_local_k, P_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_greater=P_greater/dble(nnode)   
    !
    call MPI_Barrier(MPI_COMM_WORLD)
    if (rank == 0) then  
      call write_spectrum_summed_over_k('gw_PR',iter,P_r,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    endif
    !
    deallocate(P_r,P_greater,P_lesser)
    if (rank == 0) then  
      print *, 'calc W'                
    endif
    allocate(W_lesser_local_k(nm,nm,nx,nen,nk))
    allocate(W_greater_local_k(nm,nm,nx,nen,nk))
    allocate(W_retarded_local_k(nm,nm,nx,nen,nk))
    do iq=1,nk            
      !$omp parallel default(none) private(nop) shared(rank,nk,iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser_local_k,p_greater_local_k,p_retarded_local_k,w_lesser_local_k,w_greater_local_k,w_retarded_local_k)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq+(rank)*nk),V1i(:,:,:,iq+(rank)*nk),p_lesser_local_k(:,:,:,nop,iq),p_greater_local_k(:,:,:,nop,iq),p_retarded_local_k(:,:,:,nop,iq),w_lesser_local_k(:,:,:,nop,iq),w_greater_local_k(:,:,:,nop,iq),w_retarded_local_k(:,:,:,nop,iq))
      enddo
      !$omp end do 
      !$omp end parallel      
    enddo
    deallocate(P_retarded_local_k,P_lesser_local_k,P_greater_local_k)
    allocate(W_r(nm,nm,nx,nen,nk))  
    allocate(W_greater(nm,nm,nx,nen,nk))
    allocate(W_lesser(nm,nm,nx,nen,nk))
    W_r=czero
    W_lesser=czero
    W_greater=czero   
    if (rank == 0) then
      print *,' reduce'
    endif 
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(W_retarded_local_k, W_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_r=W_r/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(W_lesser_local_k, W_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_lesser=W_lesser/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(W_greater_local_k, W_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_greater=W_greater/dble(nnode)       
    if (rank == 0) then  
        call write_spectrum_summed_over_k('gw_WR',iter,W_r       ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
        call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
        call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))    
    endif
    deallocate(W_r,W_greater,W_lesser)    
    !
    allocate(W_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(W_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(W_retarded_local_x(nm,nm,nen,nky*nkz,nnx))
    call MPI_Barrier(MPI_COMM_WORLD)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make W local in x
    !!! W(:,:,ix,:,ik=rank) --> buff(:,:,:) --> W(:,:,ix=rank,:,ik)
    allocate(sbuf(nm,nm,nen,nx,nk))
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx))        
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_retarded_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    W_retarded_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_lesser_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    W_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !   
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_greater_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    W_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)
    deallocate(sbuf,rbuf)
    deallocate(W_retarded_local_k,W_lesser_local_k,W_greater_local_k)
    !!!
    call MPI_Barrier(MPI_COMM_WORLD)
    if (rank == 0) then 
      print *, 'calc SigGW'
    endif
    !        
    allocate(g_r_local_x(nm,nm,nen,nky*nkz,nnx))  
    allocate(g_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(g_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make G local in x
    !!! G(:,:,ix,:,ik=rank) --> buff(:,:,:) --> G(:,:,ix=rank,:,ik)
    allocate(sbuf(nm,nm,nen,nx,nk))
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx))        
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_r_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    g_r_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_lesser_local_k(:,:,ix,:,ik)
      enddo
    enddo    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    g_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_greater_local_k(:,:,ix,:,ik)
      enddo
    enddo    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    g_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    call MPI_Barrier(MPI_COMM_WORLD)
    deallocate(sbuf,rbuf)
    !!!
    allocate(Sigma_lesser_new_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(Sigma_greater_new_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(Sigma_r_new_local_x(nm,nm,nen,nky*nkz,nnx))
    nopmax=nen/2-10
    Sigma_greater_new_local_x = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new_local_x = dcmplx(0.0d0,0.0d0)
    Sigma_r_new_local_x = dcmplx(0.0d0,0.0d0)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser_local_x,w_greater_local_x,w_retarded_local_x,sigma_lesser_new_local_x,sigma_greater_new_local_x,sigma_r_new_local_x,nen,En,nm,G_lesser_local_x,G_greater_local_x,G_r_local_x,nx,nkz,nky,nk,nnx)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nnx ! x via MPI parallel
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new_local_x(i,l:h,ie,ik,ix)=Sigma_lesser_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_lesser_local_x(i,l:h,ie-nop,ikd,ix) * W_lesser_local_x(i,l:h,iop,iq,ix)
                      Sigma_greater_new_local_x(i,l:h,ie,ik,ix)=Sigma_greater_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_greater_local_x(i,l:h,ie-nop,ikd,ix) * W_greater_local_x(i,l:h,iop,iq,ix)
                      Sigma_r_new_local_x(i,l:h,ie,ik,ix)=Sigma_r_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_lesser_local_x(i,l:h,ie-nop,ikd,ix) * W_retarded_local_x(i,l:h,iop,iq,ix) &
                          & + G_r_local_x(i,l:h,ie-nop,ikd,ix) * W_lesser_local_x(i,l:h,iop,iq,ix) &
                          & + G_r_local_x(i,l:h,ie-nop,ikd,ix) * W_retarded_local_x(i,l:h,iop,iq,ix)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nky*nkz) 
    Sigma_lesser_new_local_x = Sigma_lesser_new_local_x  * dE
    Sigma_greater_new_local_x= Sigma_greater_new_local_x * dE
    Sigma_r_new_local_x=Sigma_r_new_local_x* dE
    Sigma_r_new_local_x = dcmplx( dble(Sigma_r_new_local_x), aimag(Sigma_greater_new_local_x-Sigma_lesser_new_local_x)/2.0d0 )    
    deallocate(W_retarded_local_x,W_lesser_local_x,W_greater_local_x)
    deallocate(G_r_local_x,G_lesser_local_x,G_greater_local_x)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make Sigma local in k
    !!! Sig(:,:,ix=rank,:,ik) --> buff(:,:,:) --> Sig(:,:,ix,:,ik=rank)
    allocate(sigma_lesser_new_local_k(nm,nm,nx,nen,nk))
    allocate(sigma_greater_new_local_k(nm,nm,nx,nen,nk))
    allocate(sigma_r_new_local_k(nm,nm,nx,nen,nk)) 
    allocate(rbuf(nm,nm,nen,nx,nk))
    allocate(sbuf(nm,nm,nen,nky*nkz,nnx))        
    !    
    sbuf(:,:,:,:,:)=Sigma_r_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)    
    do ik=1,nk
      do ix=1,nx      
        Sigma_r_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=Sigma_lesser_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    do ik=1,nk
      do ix=1,nx
        Sigma_lesser_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=Sigma_greater_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD)
    do ik=1,nk    
      do ix=1,nx
        Sigma_greater_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD)
    deallocate(rbuf,sbuf)
    deallocate(Sigma_r_new_local_x,Sigma_lesser_new_local_x,Sigma_greater_new_local_x)
    !!!
    ! symmetrize the selfenergies
    call MPI_Barrier(MPI_COMM_WORLD)
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new_local_k(:,:,ix,ie,ik))
          Sigma_r_new_local_k(:,:,ix,ie,ik) = (Sigma_r_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new_local_k(:,:,ix,ie,ik))
          Sigma_lesser_new_local_k(:,:,ix,ie,ik) = (Sigma_lesser_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new_local_k(:,:,ix,ie,ik))
          Sigma_greater_new_local_k(:,:,ix,ie,ik) = (Sigma_greater_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw_local_k= Sigma_r_gw_local_k+ alpha_mix * (Sigma_r_new_local_k -Sigma_r_gw_local_k)
    Sigma_lesser_gw_local_k  = Sigma_lesser_gw_local_k+ alpha_mix * (Sigma_lesser_new_local_k -Sigma_lesser_gw_local_k)
    Sigma_greater_gw_local_k = Sigma_greater_gw_local_k+ alpha_mix * (Sigma_greater_new_local_k -Sigma_greater_gw_local_k)  
    ! 
    allocate(sigma_r_gw(nm,nm,nx,nen,nk))  
    allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
    allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))      
    sigma_r_gw=czero
    sigma_lesser_gw=czero
    sigma_greater_gw=czero    
    if (rank == 0) then
      print *,' reduce'
    endif
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(sigma_r_gw_local_k, sigma_r_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_r_gw=sigma_r_gw/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(sigma_lesser_gw_local_k, sigma_lesser_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_lesser_gw=sigma_lesser_gw/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD)
    call MPI_Reduce(sigma_greater_gw_local_k, sigma_greater_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_greater_gw=sigma_greater_gw/dble(nnode)
    !
    if ( rank==0 ) then
      call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    endif
    !
    if (iter>=(niter-5)) then
      if ( rank==0 ) then
        print *, 'calc SigEPhoton'
      endif      
      do ik=1,nk
        call MPI_Barrier(MPI_COMM_WORLD)
        !print *, ' ik=', ik+(rank)*nk,'/',nky*nkz
        call calc_sigma_ephoton_monochromatic(nm,nx,nen,En,nopphot,Mii(:,:,:,ik+(rank)*nk),G_lesser_local_k(:,:,:,:,ik),G_greater_local_k(:,:,:,:,ik),Sigma_r_new_local_k(:,:,:,:,ik),Sigma_lesser_new_local_k(:,:,:,:,ik),Sigma_greater_new_local_k(:,:,:,:,ik),pre_fact,intensity,hw)
      enddo      
      Sigma_r_gw_local_k       = Sigma_r_gw_local_k + Sigma_r_new_local_k 
      Sigma_lesser_gw_local_k  = Sigma_lesser_gw_local_k + Sigma_lesser_new_local_k 
      Sigma_greater_gw_local_k = Sigma_greater_gw_local_k + Sigma_greater_new_local_k 
      ! 
      sigma_r_gw=czero
      sigma_lesser_gw=czero
      sigma_greater_gw=czero    
      if (rank == 0) then
        print *,' reduce'
      endif
      call MPI_Barrier(MPI_COMM_WORLD)
      call MPI_Reduce(sigma_r_new_local_k, sigma_r_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_r_gw=sigma_r_gw/dble(nnode)
      call MPI_Barrier(MPI_COMM_WORLD)
      call MPI_Reduce(sigma_lesser_new_local_k, sigma_lesser_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_lesser_gw=sigma_lesser_gw/dble(nnode)
      call MPI_Barrier(MPI_COMM_WORLD)
      call MPI_Reduce(sigma_greater_new_local_k, sigma_greater_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_greater_gw=sigma_greater_gw/dble(nnode) 
      if ( rank==0 ) then
          call write_spectrum_summed_over_k('eph_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
          call write_spectrum_summed_over_k('eph_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
          call write_spectrum_summed_over_k('eph_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      endif
    endif
    !!! make sure self-energy is continuous near leads by copying edge block
    do ix=1,2
      Sigma_r_gw_local_k(:,:,ix,:,:)=Sigma_r_gw_local_k(:,:,3,:,:)
      Sigma_lesser_gw_local_k(:,:,ix,:,:)=Sigma_lesser_gw_local_k(:,:,3,:,:)
      Sigma_greater_gw_local_k(:,:,ix,:,:)=Sigma_greater_gw_local_k(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_r_gw_local_k(:,:,nx-2,:,:)
      Sigma_lesser_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_lesser_gw_local_k(:,:,nx-2,:,:)
      Sigma_greater_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_greater_gw_local_k(:,:,nx-2,:,:)
    enddo    
    deallocate(sigma_r_new_local_k,sigma_lesser_new_local_k,sigma_greater_new_local_k)     
    deallocate(sigma_r_gw,sigma_lesser_gw,sigma_greater_gw)
  enddo  
  deallocate(g_r_local_k,g_lesser_local_k,g_greater_local_k,cur_local_k)  
  deallocate(sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k)       
end subroutine green_rgf_solve_gw_ephoton_3d_mpi




subroutine green_rgf_solve_gw_3d(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg)
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  ! -------- local variables
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
  complex(8),allocatable,dimension(:,:,:)::Sii
  real(8)::tr(nen,nky*nkz),tre(nen,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk
  complex(8)::dE, B(nm,nm)
  !
  print *,'======== green_rgf_solve_gw_3D ========'
  nk=nky*nkz
  allocate(g_r(nm,nm,nx,nen,nk))
  allocate(g_r_i1(nm,nm,nx,nen,nk))
  allocate(g_greater(nm,nm,nx,nen,nk))
  allocate(g_lesser(nm,nm,nx,nen,nk))
  allocate(cur(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_new(nm,nm,nx,nen,nk))
  allocate(sigma_greater_new(nm,nm,nx,nen,nk))
  allocate(sigma_r_new(nm,nm,nx,nen,nk))
  allocate(P_lesser(nm,nm,nx,nen,nk))
  allocate(P_greater(nm,nm,nx,nen,nk))
  allocate(P_retarded(nm,nm,nx,nen,nk))
!  allocate(P_lesser_1i(nm,nm,nx,nen,nk))
!  allocate(P_greater_1i(nm,nm,nx,nen,nk))
!  allocate(P_retarded_1i(nm,nm,nx,nen,nk))
  allocate(W_lesser(nm,nm,nx,nen,nk))
  allocate(W_greater(nm,nm,nx,nen,nk))
  allocate(W_retarded(nm,nm,nx,nen,nk))
!  allocate(W_lesser_i1(nm,nm,nx,nen,nk))
!  allocate(W_greater_i1(nm,nm,nx,nen,nk))
!  allocate(W_retarded_i1(nm,nm,nx,nen,nk))  
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'      
    do ik=1,nk
      print *, ' ik=', ik,'/',nk
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nx,en,mu,Hii,H1i,sigma_lesser_gw,sigma_greater_gw,sigma_r_gw,g_lesser,g_greater,g_r,tre,tr,cur,g_r_i1)    
      !$omp do 
      do ie=1,nen     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik),H1i(:,:,:,ik),sigma_lesser_gw(:,:,:,ie,ik),sigma_greater_gw(:,:,:,ie,ik),&
                        & sigma_r_gw(:,:,:,ie,ik),g_lesser(:,:,:,ie,ik),g_greater(:,:,:,ie,ik),g_r(:,:,:,ie,ik),tr(ie,ik),&
                        & tre(ie,ik),cur(:,:,:,ie,ik),g_r_i1(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    !
    call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum_summed_over_kz('gw_Jdens',iter,cur,nen,En,nx,NB*NS,Lx*NS,nk)          
    call write_transmission_spectrum_k('gw_trR',iter,tr*spindeg,nen,En,nk)
    call write_transmission_spectrum_k('gw_trL',iter,tre*spindeg,nen,En,nk)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10      
    P_lesser = dcmplx(0.0d0,0.0d0)
    P_greater = dcmplx(0.0d0,0.0d0)    
    P_retarded = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i = dcmplx(0.0d0,0.0d0)
!    P_greater_1i = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nx
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
              do iky=1,nky
                do ikz=1,nkz              
                  ik=ikz + (iky-1)*nkz
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky                
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1
                  ikd=ikzd + (ikyd-1)*nkz
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser(i,l:h,ix,iop,iq) = P_lesser(i,l:h,ix,iop,iq) + G_lesser(i,l:h,ix,ie,ik) * G_greater(l:h,i,ix,ie-nop,ikd)
                    P_greater(i,l:h,ix,iop,iq) = P_greater(i,l:h,ix,iop,iq) + G_greater(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd)        
                    P_retarded(i,l:h,ix,iop,iq) = P_retarded(i,l:h,ix,iop,iq) + (G_lesser(i,l:h,ix,ie,ik) * conjg(G_r(i,l:h,ix,ie-nop,ikd)) & 
                                              & + G_r(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo      
    enddo          
    !$omp end do 
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum_summed_over_k('gw_PR',iter,P_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'                
    do iq=1,nky*nkz
      print *, ' iq=', iq,'/',nk
      !$omp parallel default(none) private(nop) shared(iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,w_lesser,w_greater,w_retarded)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    call write_spectrum_summed_over_k('gw_WR',iter,W_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nx
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new(i,l:h,ix,ie,ik)=Sigma_lesser_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq)
                      Sigma_greater_new(i,l:h,ix,ie,ik)=Sigma_greater_new(i,l:h,ix,ie,ik)+G_greater(i,l:h,ix,ie-nop,ikd)*W_greater(i,l:h,ix,iop,iq)
                      Sigma_r_new(i,l:h,ix,ie,ik)=Sigma_r_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! symmetrize the selfenergies
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new(:,:,ix,ie,ik))
          Sigma_r_new(:,:,ix,ie,ik) = (Sigma_r_new(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new(:,:,ix,ie,ik))
          Sigma_lesser_new(:,:,ix,ie,ik) = (Sigma_lesser_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new(:,:,ix,ie,ik))
          Sigma_greater_new(:,:,ix,ie,ik) = (Sigma_greater_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:,:)=Sigma_r_gw(:,:,3,:,:)
      Sigma_lesser_gw(:,:,ix,:,:)=Sigma_lesser_gw(:,:,3,:,:)
      Sigma_greater_gw(:,:,ix,:,:)=Sigma_greater_gw(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:,:)=Sigma_r_gw(:,:,nx-2,:,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:,:)=Sigma_lesser_gw(:,:,nx-2,:,:)
      Sigma_greater_gw(:,:,nx-ix+1,:,:)=Sigma_greater_gw(:,:,nx-2,:,:)
    enddo
    call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
end subroutine green_rgf_solve_gw_3d


subroutine green_rgf_solve_gw_1d(alpha_mix,niter,NB,NS,nm,nx,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg)
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx)::Hii,H1i,Vii,V1i
  ! -------- local variables
  complex(8),allocatable,dimension(:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
  complex(8),allocatable,dimension(:,:,:)::Sii
  real(8)::tr(nen),tre(nen)
  integer::ie,iter,i,ix,nopmax,nop,iop,l,h,io
  complex(8)::dE
  !
  print *,'====== green_rgf_solve_gw_1D ======'
  allocate(g_r(nm,nm,nx,nen))
  allocate(g_r_i1(nm,nm,nx,nen))
  allocate(g_greater(nm,nm,nx,nen))
  allocate(g_lesser(nm,nm,nx,nen))
  allocate(cur(nm,nm,nx,nen))
  allocate(sigma_lesser_gw(nm,nm,nx,nen))
  allocate(sigma_greater_gw(nm,nm,nx,nen))
  allocate(sigma_r_gw(nm,nm,nx,nen))
  allocate(sigma_lesser_new(nm,nm,nx,nen))
  allocate(sigma_greater_new(nm,nm,nx,nen))
  allocate(sigma_r_new(nm,nm,nx,nen))
  allocate(P_lesser(nm,nm,nx,nen))
  allocate(P_greater(nm,nm,nx,nen))
  allocate(P_retarded(nm,nm,nx,nen))
  !allocate(P_lesser_1i(nm,nm,nx,nen))
  !allocate(P_greater_1i(nm,nm,nx,nen))
  !allocate(P_retarded_1i(nm,nm,nx,nen))
  allocate(W_lesser(nm,nm,nx,nen))
  allocate(W_greater(nm,nm,nx,nen))
  allocate(W_retarded(nm,nm,nx,nen))
  !allocate(W_lesser_i1(nm,nm,nx,nen))
  !allocate(W_greater_i1(nm,nm,nx,nen))
  !allocate(W_retarded_i1(nm,nm,nx,nen))  
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'      
    do ie=1,nen     
      if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii,H1i,sigma_lesser_gw(:,:,:,ie),sigma_greater_gw(:,:,:,ie),&
                      & sigma_r_gw(:,:,:,ie),g_lesser(:,:,:,ie),g_greater(:,:,:,ie),g_r(:,:,:,ie),tr(ie),&
                      & tre(ie),cur(:,:,:,ie),g_r_i1(:,:,:,ie)) 
    enddo
    !
    call write_spectrum('gw_ldos',iter,g_r,nen,En,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum('gw_ndos',iter,g_lesser,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum('gw_pdos',iter,g_greater,nen,En,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum('gw_Jdens',iter,cur,nen,En,nx,NB*NS,Lx*NS)          
    call write_transmission_spectrum('gw_trR',iter,tr(:)*spindeg,nen,En)
    call write_transmission_spectrum('gw_trL',iter,tre(:)*spindeg,nen,En)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10  
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg   
    P_lesser(:,:,:,:) = dcmplx(0.0d0,0.0d0)
    P_greater(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
    P_retarded(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)
!    P_greater_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx)    
    !$omp do 
    do ix=1,nx
      do i=1,nm      
        l=max(i-ndiag,1)
        h=min(nm,i+ndiag)   
        do nop=-nopmax,nopmax
          iop=nop+nen/2  
          do ie = max(nop+1,1),min(nen,nen+nop)           
              P_lesser(i,l:h,ix,iop) = P_lesser(i,l:h,ix,iop) + G_lesser(i,l:h,ix,ie) * G_greater(l:h,i,ix,ie-nop)
              P_greater(i,l:h,ix,iop) = P_greater(i,l:h,ix,iop) + G_greater(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop)        
              P_retarded(i,l:h,ix,iop) = P_retarded(i,l:h,ix,iop) + (G_lesser(i,l:h,ix,ie) * conjg(G_r(i,l:h,ix,ie-nop)) & 
                                       & + G_r(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop))        
          enddo
        enddo
      enddo
    enddo
    !$omp end do 
    !$omp end parallel
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum('gw_PR',iter,P_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_PG',iter,P_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'
    do ie=1,nen
      if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      !call green_calc_w_full(0,nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie))
      !call green_rgf_calc_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie),w_lesser_i1(:,:,:,ie),w_greater_i1(:,:,:,ie),w_retarded_i1(:,:,:,ie))
      !call green_rgf_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie))
    enddo
    call write_spectrum('gw_WR',iter,W_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_WG',iter,W_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx)    
    !$omp do
    do ix=1,nx
      do i=1,nm
        l=max(i-ndiag,1)
        h=min(nm,i+ndiag)  
        do ie=1,nen      
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then  
              iop=nop+nen/2                
              Sigma_lesser_new(i,l:h,ix,ie)=Sigma_lesser_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_lesser(i,l:h,ix,iop)
              Sigma_greater_new(i,l:h,ix,ie)=Sigma_greater_new(i,l:h,ix,ie)+G_greater(i,l:h,ix,ie-nop)*W_greater(i,l:h,ix,iop)
              Sigma_r_new(i,l:h,ix,ie)=Sigma_r_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_retarded(i,l:h,ix,iop) &
                                                            &  +G_r(i,l:h,ix,ie-nop)*W_lesser(i,l:h,ix,iop) &
                                                            &  +G_r(i,l:h,ix,ie-nop)*W_retarded(i,l:h,ix,iop)    
            endif
          enddo
        enddo
      enddo    
    enddo
    !$omp end do
    !$omp end parallel
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:)=Sigma_r_gw(:,:,3,:)
      Sigma_lesser_gw(:,:,ix,:)=Sigma_lesser_gw(:,:,3,:)
      Sigma_greater_gw(:,:,ix,:)=Sigma_greater_gw(:,:,3,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:)=Sigma_r_gw(:,:,nx-2,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:)=Sigma_lesser_gw(:,:,nx-2,:)
      Sigma_greater_gw(:,:,nx-ix+1,:)=Sigma_greater_gw(:,:,nx-2,:)
    enddo
    call write_spectrum('gw_SigR',iter,Sigma_r_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_SigL',iter,Sigma_lesser_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_SigG',iter,Sigma_greater_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
end subroutine green_rgf_solve_gw_1d



subroutine full2block(fullA,A,NB,length,offdiag)
complex(8),intent(in)::fullA(NB*length,NB*length)
complex(8),intent(out)::A(NB,NB,length)
integer,intent(in)::NB,length
integer,intent(in),optional::offdiag
integer::i,l,h
A=0.0d0
do i=1,length
  l=(i-1)*NB+1
  h=i*NB
  if ((present(offdiag)).and.((i+offdiag) <= length)) then ! n-th off-diagonal blocks
    A(:,:,i)=fullA(l:h,l+nb*offdiag:h+nb*offdiag)
  else ! diagonal blocks
    A(:,:,i)=fullA(l:h,l:h)
  endif
enddo
end subroutine full2block


subroutine block2full(A,fullA,NB,length,offdiag)
complex(8),intent(in)::A(NB,NB,length)
complex(8),intent(out)::fullA(NB*length,NB*length)
integer,intent(in)::NB,length
integer,intent(in),optional::offdiag
integer::i,l,h
do i=1,length
  l=(i-1)*NB+1
  h=i*NB  
  if (present(offdiag)) then
    if ((i+offdiag) <= length) then ! n-th off-diagonal blocks
      fullA(l:h,l+nb*offdiag:h+nb*offdiag) = transpose(conjg(A(:,:,i)))
      fullA(l+nb*offdiag:h+nb*offdiag,l:h) = A(:,:,i)
    endif
  else ! diagonal blocks
    fullA(l:h,l:h)=A(:,:,i)
  endif
enddo
end subroutine block2full


subroutine green_calc_w_full(NBC,nm,nx,Vii,V1i,PLii,PGii,PRii,WLii,WGii,WRii)
integer,intent(in)::nm,nx,NBC
!complex(8),intent(in),dimension(nm*nx,nm*nx) :: V
complex(8),intent(in),dimension(nm,nm,nx) :: PLii,PGii,PRii
complex(8),intent(in),dimension(nm,nm,nx) :: Vii, V1i
!complex(8),intent(in),dimension(nm,nm,nx) :: PL1i,PG1i,PR1i
complex(8),intent(inout),dimension(nm,nm,nx) :: WLii,WGii,WRii!,WLi1,WGi1,WRi1
! --------- local
complex(8),allocatable,dimension(:,:)::B,M,V
complex(8),allocatable,dimension(:,:)::PL,PR,PG, WL,WR,WG
integer::i,NT
NT=nm*nx! total size
allocate(V(NT,NT))
allocate(PR(NT,NT))
allocate(PL(NT,NT))
allocate(PG(NT,NT))
allocate(WR(NT,NT))
allocate(WL(NT,NT))
allocate(WG(NT,NT))
allocate(B(NT,NT))
allocate(M(NT,NT))  
PR=czero
PL=czero
PG=czero
V=czero
call block2full(Vii,V,nm,nx)  
call block2full(V1i,V,nm,nx,1)  
call block2full(PRii,PR,nm,nx)  
call block2full(PLii,PL,nm,nx)  
call block2full(PGii,PG,nm,nx)  
!! M = I - V P^r
call zgemm('n','n',NT,NT,NT,-cone,V,NT,PR,NT,czero,M,NT)
!
do i=1,NT
   M(i,i) = 1.0d0 + M(i,i)
enddo    
!!!! calculate W^r = (I - V P^r)^-1 V    
call invert(M,NT) ! M -> xR
call zgemm('n','n',NT,NT,NT,cone,M,NT,V,NT,czero,WR,NT)           
! calculate W^< and W^> = W^r P^<> W^r dagger
call zgemm('n','n',NT,NT,NT,cone,WR,NT,PL,NT,czero,B,NT) 
call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WL,NT) 
call zgemm('n','n',NT,NT,NT,cone,WR,NT,PG,NT,czero,B,NT) 
call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WG,NT)  
call full2block(WR,WRii,nm,nx)
call full2block(WL,WLii,nm,nx)
call full2block(WG,WGii,nm,nx)
deallocate(B,M,V)  
deallocate(PR,PL,PG,WR,WL,WG)
end subroutine green_calc_w_full




subroutine identity(A,n)
  integer, intent(in) :: n        
  complex(8), dimension(n,n), intent(inout) :: A
  integer :: i
  A = dcmplx(0.0d0,0.0d0)
  do i = 1,n
    A(i,i) = dcmplx(1.0d0,0.0d0)
  end do
end subroutine identity



!!!! RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_RS(TEMP,nm,nx,E,mu,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,GRi1,GLi1,GGi1)     
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: cur,ldos,ndens,pdens ! diag blocks of GFs
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    ! ------- 
    COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    COMPLEX(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr
    COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
    REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
    mul=mu(1)
    mur=mu(2)
    templ=temp(1)
    tempr=temp(2)
    !
    z=cmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    S00=0.0d0
    do i=1,nm
        S00(i,i)=1.0d0
    enddo
    ! self energy on the left contact        
    H00(:,:)=Hii(:,:,1)+sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)      
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+sigma_lesser_ph(:,:,1)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    Do l=2,nx-1                
        H00(:,:)=Hii(:,:,l)+sigma_r_ph(:,:,l)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)             
        C(:,:)=C(:,:)+sigma_lesser_ph(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact        
    H00(:,:)=Hii(:,:,nx)+sigma_r_ph(:,:,nx)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)!cmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01    
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+sigma_lesser_ph(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        do j=i,i!1,nm
            tim=tim-dble(B(i,j)-C(i,j))
        enddo
    enddo
    tr=tim
    ! transmission
    !-------------------------
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
        !-------------------------
!        A(:,:)=Gp(:,:)
!        call zgemm('c','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
!        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
!        !
!        Gp(:,:)= Glp(:,:,l) + C(:,:)
!        A(:,:)=Glp(:,:,l)
!        call zgemm('c','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
!        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
!        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
        pdens(:,:,l)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
    deallocate(Gl)
    deallocate(Gln)
    deallocate(Glp)
end subroutine green_RGF_RS

!!!! Coupled-Mode-Space RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_CMS(TEMP,nm,nx,E,mu,Sii,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,GRi1,GLi1,GGi1)         
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Sii,Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: cur,ldos,ndens,pdens ! diag blocks of GFs
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    ! ------- 
    COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    COMPLEX(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr
    COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
    REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
    mul=mu(1)
    mur=mu(2)
    templ=temp(1)
    tempr=temp(2)
    !
    z=cmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    ! self energy on the left contact
    S00(:,:)=Sii(:,:,1)
    call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,1),nm,S00,nm,beta,B,nm)
    H00(:,:)=Hii(:,:,1)+B(:,:)!sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)  
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,1),nm,S00,nm,beta,B,nm)
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+B(:,:)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    Do l=2,nx-1
        S00(:,:)=Sii(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,l),nm,S00,nm,beta,B,nm)
        H00(:,:)=Hii(:,:,l)+B(:,:)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,l),nm,S00,nm,beta,B,nm)
        C(:,:)=C(:,:)+B(:,:)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact
    S00(:,:)=Sii(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,nx),nm,S00,nm,beta,B,nm)
    H00(:,:)=Hii(:,:,nx)+B(:,:)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)!cmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,nx),nm,S00,nm,beta,B,nm)
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+B(:,:)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        do j=i,i!1,nm
            tim=tim-dble(B(i,j)-C(i,j))
        enddo
    enddo
    tr=tim
    ! transmission
    !-------------------------
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
        !-------------------------
        A(:,:)=Gp(:,:)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        !
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        !
        Gp(:,:)= Glp(:,:,l) + C(:,:)
        A(:,:)=Glp(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        !
        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
        pdens(:,:,l)=Gp(:,:)
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
    deallocate(Gl)
    deallocate(Gln)
    deallocate(Glp)
end subroutine green_RGF_CMS


! find the inverse of a band matrix A by solving a system of linear equations
! on exit, A contains the band matrix of inv(A)
subroutine invert_banded(A,nn,nb)
integer,intent(in)::nn,nb
complex(8),intent(inout)::A(3*nb+1,nn)
complex(8),allocatable::work(:),B(:,:),X(:,:)
integer,allocatable::ipiv(:)
integer::info,lda,lwork,ldb,i,nrhs
allocate(ipiv(nn))
allocate(work(nn*nn))
lda=3*nb+1
call zgbtrf(nn,nn,nb,nb,A,lda,ipiv,info)
if (info.ne.0) then
  print*,'SEVERE warning: zgbtrf failed, info=',info
  call abort
endif
ldb=1
allocate(B(ldb,nn))
allocate(X(lda,nn))
nrhs=ldb
do i=1,nn
  B=0.0d0
  B(1,i)=1.0d0
  call zgbtrs('N',nn,nb,nb,nrhs,A,lda,ipiv,B,ldb,info)
  if (info.ne.0) then
    print*,'SEVERE warning: zgbtrs failed, info=',info
    call abort
  endif
  X(1:nb*2+1,i)=B(1,i-nb:i+nb)
enddo
A=X
deallocate(B,work,ipiv,X)
end subroutine invert_banded


subroutine invert(A,nn)  
  integer :: info,lda,lwork,nn      
  integer, dimension(:), allocatable :: ipiv
  complex(8), dimension(nn,nn),intent(inout) :: A
  complex(8), dimension(:), allocatable :: work
  allocate(work(nn*nn))
  allocate(ipiv(nn))
  call zgetrf(nn,nn,A,nn,ipiv,info)
  if (info.ne.0) then
    print*,'SEVERE warning: zgetrf failed, info=',info
    A=czero
  else
    call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
    if (info.ne.0) then
      print*,'SEVERE warning: zgetri failed, info=',info
      A=czero
    endif
  endif
  deallocate(work)
  deallocate(ipiv)
end subroutine invert


Function ferm(a)
    Real (8) a,ferm
    ferm=1.0d0/(1.0d0+Exp(a))
End Function ferm


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)    
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


! write spectrum into file (pm3d map)
subroutine write_spectrum(dataset,i,G,nen,en,length,NB,NS,Lx,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(NB*NS,NB*NS,length,nen)
integer, intent(in)::i,nen,length,NB,NS
real(8), intent(in)::Lx,en(nen),coeff(2)
integer:: ie,j,ib,k
complex(8)::tr
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length
      do k=1,NS         
        tr=0.0d0         
        do ib=1,nb
            tr = tr+ G(ib+(k-1)*NB,ib+(k-1)*NB,j,ie)            
        end do
        write(11,'(4E18.4)') (j-1)*Lx*NS+(k-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
      enddo
    end do
    write(11,*)    
end do
close(11)
end subroutine write_spectrum


! write current spectrum into file (pm3d map)
subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: cur(:,:,:,:)
  integer, intent(in)::i,nen,length,NB
  real(8), intent(in)::Lx,en(nen)
  integer:: ie,j,ib,jb
  real(8)::tr
  open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
    do j = 1,length-1
      tr=0.0d0                
      do ib=1,nb  
        do jb=1,nb        
          tr = tr+ cur(ib,jb,j,ie)
        enddo                        
      enddo      
      write(11,'(3E18.4)') dble(j-1)*Lx, en(ie), dble(tr)
    enddo
    write(11,*)    
  enddo
  close(11)
end subroutine write_current_spectrum

! write current spectrum into file (pm3d map)
subroutine write_current_spectrum_summed_over_kz(dataset,i,cur,nen,en,length,NB,Lx,nk)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: cur(:,:,:,:,:)
  integer, intent(in)::i,nen,length,NB,nk
  real(8), intent(in)::Lx,en(nen)
  integer:: ie,j,ib,jb,ik
  real(8)::tr
  open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
    do j = 1,length-1
      tr=0.0d0          
      do ik=1,nk
        do ib=1,nb  
          do jb=1,nb        
            tr = tr+ cur(ib,jb,j,ie,ik)
          enddo                        
        enddo
      enddo
      write(11,'(3E18.4)') dble(j-1)*Lx, en(ie), dble(tr)
    enddo
    write(11,*)    
  enddo
  close(11)
end subroutine write_current_spectrum_summed_over_kz


! write spectrum into file (pm3d map)
subroutine write_spectrum_summed_over_k(dataset,i,G,nen,en,nk,length,NB,NS,Lx,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(NB*NS,NB*NS,length,nen,nk)
integer, intent(in)::i,nen,length,NB,NS,nk
real(8), intent(in)::Lx,en(nen),coeff(2)
integer:: ie,j,ib,k,ik
complex(8)::tr
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
  do j = 1,length
    do k=1,NS         
      tr=0.0d0         
      do ik=1,nk
        do ib=1,nb
          tr = tr+ G(ib+(k-1)*NB,ib+(k-1)*NB,j,ie,ik)            
        end do
      enddo
      tr=tr/dble(nk)
      write(11,'(4E18.4)') (j-1)*Lx*NS+(k-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    enddo
  end do
  write(11,*)    
end do
close(11)
end subroutine write_spectrum_summed_over_k

! write transmission spectrum into file
subroutine write_transmission_spectrum(dataset,i,tr,nen,en)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: tr(:)
integer, intent(in)::i,nen
real(8), intent(in)::en(nen)
integer:: ie,j,ib
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen    
  write(11,'(2E18.4)') en(ie), dble(tr(ie))      
end do
close(11)
end subroutine write_transmission_spectrum


! write transmission spectrum into file
subroutine write_transmission_spectrum_k(dataset,i,tr,nen,en,nk)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: tr(:,:)
integer, intent(in)::i,nen,nk
real(8), intent(in)::en(nen)
integer:: ie,j,ib,ik
real(8):: sumtr(nen)
sumtr=0.0d0
do ik=1,nk
  !open(unit=11,file=trim(dataset)//'kz'//TRIM(STRING(ik))//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  !do ie = 1,nen    
  !  write(11,'(2E18.4)') en(ie), dble(tr(ie,ik))      
  !end do
  !close(11)
  sumtr=sumtr+tr(:,ik)
enddo
sumtr=sumtr/dble(nk)
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen    
  write(11,'(2E18.4)') en(ie), dble(sumtr(ie))      
end do
close(11)
end subroutine write_transmission_spectrum_k


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

subroutine open_boundary_conditions(nm,M00,M01,M10,V10,xR,dM,dV,cond)
integer,intent(in)::nm
complex(8),intent(in),dimension(nm,nm)::M00,M01,M10,V10
complex(8),intent(out),dimension(nm,nm)::xR,dM,dV
real(8),intent(out)::cond
complex(8),dimension(nm,nm)::tmp1
call surface_function(nm,M00,M01,M10,xR,cond);
!dM=M01*xR*M10
call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dM,nm)
!dV=M01*xR*V10
call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,V10,nm,czero,dV,nm)
end subroutine open_boundary_conditions


! a slightly modified version of sancho
subroutine surface_function(nm,M00,M01,M10,SF,cond)
integer,intent(in)::nm
complex(8),intent(in),dimension(nm,nm)::M00,M01,M10
complex(8),intent(out),dimension(nm,nm)::SF
real(8),intent(out)::cond
real(8)::cond_limit
integer::max_iteration,IC
complex(8),dimension(:,:),allocatable::alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta
allocate(alpha(nm,nm))
allocate(beta(nm,nm))
allocate(Eps(nm,nm))
allocate(Eps_surf(nm,nm))
allocate(inv_element(nm,nm))
allocate(i_alpha(nm,nm))
allocate(i_beta(nm,nm))
allocate(a_i_b(nm,nm))
allocate(b_i_a(nm,nm))
cond=1.0d10;
cond_limit=1.0d-10;
max_iteration=5000;
IC=1;
alpha=M01
beta=M10
Eps=M00
Eps_surf=M00
do while ((cond>cond_limit).and.(IC<max_iteration))      
    inv_element=Eps
    call invert(inv_element,nm)
    i_alpha=matmul(inv_element,alpha)
    i_beta=matmul(inv_element,beta)
    a_i_b=matmul(alpha,i_beta)
    b_i_a=matmul(beta,i_alpha)
    Eps=Eps-a_i_b-b_i_a
    Eps_surf=Eps_surf-a_i_b
    alpha=matmul(alpha,i_alpha)
    beta=matmul(beta,i_beta)
    !
    cond=sum(abs(alpha)+abs(beta))/2.0d0;
    !
    IC=IC+1;
end do
if (cond>cond_limit) then 
  write(*,*) 'SEVERE warning: nmax reached in surface function!!!',cond
endif
call invert(Eps_surf,nm)
SF=Eps_surf
deallocate(alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta)
end subroutine surface_function


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
!
CALL zheev( 'V','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )
!
deallocate(work,rwork)
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
   call abort
endif
eigv(:)=W(:)
END FUNCTION eigv

end module green_rgf
