!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf
use green,only:get_OBC_blocks_for_W,get_dL_OBC_for_W    

implicit none 

private

public :: green_rgf_cms, green_rgf_solve_gw_1d, green_rgf_calc_w
public :: green_rgf_solve_gw_3d
public :: green_rgf_solve_gw_ephoton_3d

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
    call write_transmission_spectrum_k('gw_trL',iter,tr*spindeg,nen,En,nk)
    call write_transmission_spectrum_k('gw_trR',iter,tre*spindeg,nen,En,nk)
    open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
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
    call write_spectrum_summed_over_k('gw_Jdens',iter,cur,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))         
    call write_transmission_spectrum_k('gw_trL',iter,tr*spindeg,nen,En,nk)
    call write_transmission_spectrum_k('gw_trR',iter,tre*spindeg,nen,En,nk)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
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
    call write_spectrum('gw_Jdens',iter,cur,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))         
    call write_transmission_spectrum('gw_trL',iter,tr(:)*spindeg,nen,En)
    call write_transmission_spectrum('gw_trR',iter,tre(:)*spindeg,nen,En)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg), sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
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


subroutine prepare_rgf_W(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,Hii,H1i,Hi1,sigma_lesser,sigma_greater)
integer,intent(in)::nm,nx
complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,P_lesser,P_greater,P_retarded
complex(8),intent(in),dimension(nm,nm,nx) :: P_lesser_1i,P_greater_1i,P_retarded_1i
complex(8),intent(inout),dimension(nm,nm,nx) :: Hii,H1i,Hi1,sigma_lesser,sigma_greater
!----------
complex(8),dimension(nm,nm) :: A
integer::i
Hii=czero
H1i=czero
Hi1=czero
do i=1,nx
  ! Hi,i = Vi,i Pi,i + Vi,i+1 Pi+1,i + Vi,i-1 Pi-1,i  
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_retarded(:,:,i),nm,cone,Hii(:,:,i),nm) 
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_retarded_1i(:,:,i),nm,cone,Hii(:,:,i),nm) 
  if (i==1) then
    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,i),nm,p_retarded_1i(:,:,i),nm,cone,Hii(:,:,i),nm) 
  else
    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,i-1),nm,p_retarded_1i(:,:,i-1),nm,cone,Hii(:,:,i),nm) 
  endif
  ! Hi+1,i = Vi+1,i+1 Pi+1,i + Vi+1,i Pi,i
  if (i==nx) then
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_retarded_1i(:,:,i),nm,cone,H1i(:,:,i),nm) 
  else    
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,p_retarded_1i(:,:,i),nm,cone,H1i(:,:,i),nm) 
  endif
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_retarded(:,:,i),nm,cone,H1i(:,:,i),nm) 
  ! Hi,i+1 = Vi,i Pi,i+1 + Vi,i+1 Pi+1,i+1
  call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,p_retarded_1i(:,:,i),nm,cone,Hi1(:,:,i),nm)
  if (i==nx) then 
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_retarded(:,:,i),nm,cone,Hi1(:,:,i),nm) 
  else
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_retarded(:,:,i+1),nm,cone,Hi1(:,:,i),nm) 
  endif
  ! Sig<> = V P<> V'
  ! Sig<_i,i
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser(:,:,i),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
!  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_lesser(:,:,min(i+1,nx)),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
!  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_lesser_1i(:,:,i),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
!  !
!  call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser_1i(:,:,i),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
!  !
!  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser_1i(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_lesser(:,:,i),nm)
!  !
!  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_lesser(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_lesser(:,:,i),nm)
!  !
!  call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_lesser_1i(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  ! Sig>_i,i
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater(:,:,i),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
!  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_greater(:,:,min(i+1,nx)),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
!  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_greater_1i(:,:,i),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
!  !
!  call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater_1i(:,:,i),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
!  !
!  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater_1i(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_greater(:,:,i),nm)
!  !
!  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_greater(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_greater(:,:,i),nm)
!  !
!  call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_greater_1i(:,:,max(i-1,1)),nm,czero,A,nm)
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
enddo
end subroutine prepare_rgf_W


! RGF for diagonal blocks of screened Coulomb GF W^r,<,>
! W^r = [I-VPr]^-1 V
! W^<> = W^r P<> W^a
! compare to usual RGF: 
!     V.Pr as Ham, and V.P<>.V' as Sig<>
!     Wr = Gr.V  
!     W<> = G<>
subroutine green_rgf_w(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,w_lesser,w_greater,w_retarded)
  integer,intent(in)::nm,nx
  complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,P_lesser,P_greater,P_retarded
  complex(8),intent(in),dimension(nm,nm,nx) :: P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),intent(inout),dimension(nm,nm,nx) :: W_lesser,W_greater,W_retarded
  ! --------
  COMPLEX(8) :: S00(nm,nm),H00(nm,nm),H10(nm,nm),H01(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),&
      GN0(nm,nm),G0N(nm,nm),Gn(nm,nm),Gp(nm,nm),G00(nm,nm),sigmal(nm,nm),sigmar(nm,nm),sig(nm,nm)       
  COMPLEX(8) :: z
  integer::i,j,k,l
  COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected G
  COMPLEX(8), allocatable,dimension(:,:,:) :: Hii,H1i,Hi1,sigma_lesser,sigma_greater
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  !
  allocate(Gl(nm,nm,nx))
  allocate(Gln(nm,nm,nx))
  allocate(Glp(nm,nm,nx))
  allocate(sigma_lesser(nm,nm,nx))
  allocate(sigma_greater(nm,nm,nx))    
  allocate(Hii(nm,nm,nx))    
  allocate(Hi1(nm,nm,nx))    
  allocate(H1i(nm,nm,nx))    
  !
  Gln=0.0D0
  Glp=0.0D0
  Gl=0.0D0
  w_lesser=0.0d0
  w_greater=0.0d0
  w_retarded=0.0d0    
  !
  S00=0.0d0
  do i=1,nm
      S00(i,i)=1.0d0
  enddo
  ! prepare
  ! H = V.Pr
  ! Sig<> = V.P<>.V'
  !
  call prepare_rgf_W(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,Hii,H1i,Hi1,sigma_lesser,sigma_greater)
  !
  ! start RGF
  z=dcmplx(1.0d0, 0.0)
  ! on the left contact        
  H00(:,:)=Hii(:,:,1)
  H10(:,:)=H1i(:,:,1)  
  !
  A=z*S00-H00
  !                
  call invert(A,nm)
  Gl(:,:,1)=A(:,:)
  !
  sig=sigma_lesser(:,:,1)
  call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
  call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
  Gln(:,:,1)=C(:,:)
  !
!  call zgemm('n','n',nm,nm,nm,alpha,sigma_greater(:,:,1),nm,S00,nm,beta,B,nm)
!  sig=B
!  call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
!  call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
!  Glp(:,:,1)=C(:,:)
  Do l=1,nx-1                
    H00(:,:)=Hii(:,:,l)
    H10(:,:)=H1i(:,:,l)
    H01(:,:)=Hi1(:,:,l)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,max(l-1,1)),nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)
    A=z*S00-H00-C   
    !
    call invert(A,nm)
    Gl(:,:,l)=A(:,:)
    !
    sig=Gln(:,:,max(l-1,1))
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)     
    !    
    C(:,:)=C(:,:)+sigma_lesser(:,:,l)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
    Gln(:,:,l)=Gn(:,:)
!    !
!    sig=Glp(:,:,max(l-1,1))
!    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
!    call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)     
!    call zgemm('n','n',nm,nm,nm,alpha,sigma_greater(:,:,l),nm,S00,nm,beta,B,nm)
!    C(:,:)=C(:,:)+B(:,:)
!    call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
!    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gp,nm)
!    Glp(:,:,l)=Gp(:,:)
  enddo
  !
  H00(:,:)=Hii(:,:,nx)
  H10(:,:)=H1i(:,:,nx)
  H01(:,:)=Hi1(:,:,nx)    
  call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
  call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)
  G00=z*S00-H00-C   
  !
  call invert(G00,nm)
  !
  call zgemm('n','n',nm,nm,nm,alpha,G00,nm,Vii(:,:,nx),nm,beta,w_retarded(:,:,nx),nm)
  !
  sig=Gln(:,:,nx-1)
  call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
  call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)  ! C=H10 Gl< H01    
  sig(:,:)=C(:,:)+sigma_lesser(:,:,nx)
  call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
  call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
  ! G<00 = G00 sig< G00'
  w_lesser(:,:,nx)=Gn(:,:)    
  !
!  sig=Glp(:,:,nx-1)
!  call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
!  call zgemm('n','n',nm,nm,nm,alpha,B,nm,H01,nm,beta,C,nm)  ! C=H10 Gl> H01
!  call zgemm('n','n',nm,nm,nm,alpha,sigma_greater(:,:,nx),nm,S00,nm,beta,B,nm)
!  ! B=Sig> S00
!  sig(:,:)=C(:,:)+B(:,:)
!  call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
!  call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gp,nm) 
!  ! G>00 = G00 sig> G00'
  Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
  w_greater(:,:,nx)=Gp(:,:)      
  !-------------------------
  do l=nx-1,1,-1
    H10(:,:)=H1i(:,:,l)
    H01(:,:)=Hi1(:,:,l)
    A=Gn
    call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
    A=Gln(:,:,l)          
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
    B=C+A
    call zgemm('n','n',nm,nm,nm,alpha,H01,nm,B,nm,beta,A,nm)      !!! G<_i+1,i                     
    D(:,:)= Gl(:,:,l)
    call zgemm('n','n',nm,nm,nm,alpha,D,nm,H01,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
    !
    call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,V1i(:,:,l),nm,alpha,w_retarded(:,:,l),nm)        
    !
    call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
    G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,Vii(:,:,l),nm,alpha,w_retarded(:,:,l),nm)                
    !
    D(:,:)=Gl(:,:,l+1)
    call zgemm('n','n',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,G0N,nm)      !!! G_i+1,i
    !
    call zgemm('n','c',nm,nm,nm,alpha,G0N,nm,V1i(:,:,l),nm,alpha,w_retarded(:,:,l+1),nm)
    !-------------------------
    A(:,:)=Gn(:,:)     
    call zgemm('n','n',nm,nm,nm,alpha,D,nm,H01,nm,beta,B,nm)  
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
    call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
    Gn(:,:)= Gln(:,:,l) + C(:,:)
    A(:,:)=Gln(:,:,l)
    call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
    Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H01,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,G0N,nm,beta,C,nm)          
    Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
    !-------------------------
    A(:,:)=Gp(:,:)
    call zgemm('n','n',nm,nm,nm,alpha,D,nm,H01,nm,beta,B,nm)  
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
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H01,nm,beta,B,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,B,nm,G0N,nm,beta,C,nm)     
    Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
    !-------------------------    
    w_lesser(:,:,l)=Gn(:,:)
    w_greater(:,:,l)=Gp(:,:)
  enddo
  !    
  deallocate(Gl,Gln,Glp)
  deallocate(sigma_lesser,sigma_greater)
  deallocate(Hii,H1i,Hi1)
end subroutine green_rgf_w


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


!!!! RGF for diagonal blocks of screened Coulomb GF W^r,<,>
! W^r = [I-VPr]^-1 V
! W^<> = W^r P<> W^a
! compare to usual RGF: 
!     V.Pr as Ham, and V.P<>.V' as Sig<>
!     Wr  ~ Gr.V  
!     W<> ~ G<>
subroutine green_rgf_calc_w(NBC,nm,nx,Vii,V1i,PL,PG,PR,PL1i,PG1i,PR1i,WL,WG,WR,WLi1,WGi1,WRi1)
integer,intent(in)::nm,nx
complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,PL,PG,PR
complex(8),intent(in),dimension(nm,nm,nx) :: PL1i,PG1i,PR1i
complex(8),intent(inout),dimension(nm,nm,nx) :: WL,WG,WR,WLi1,WGi1,WRi1
integer,intent(in) :: NBC
! -------- local
COMPLEX(8)::A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),Wn(nm,nm),Wp(nm,nm),sig(nm,nm),AL(nm,nm),BL(nm,nm)
COMPLEX(8)::S00(nm,nm),H00(nm,nm),H10(nm,nm),H01(nm,nm),GN0(nm,nm),G0N(nm,nm),G00(nm,nm)
complex(8),allocatable,dimension(:,:,:)::S,M,LL,LG,VV,M1i,Mi1,S1i,Si1
complex(8),allocatable,dimension(:,:,:)::xlr,Wln,Wlp,xR,Wlr
complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
    PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
    MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
integer::i,NL,NR,LBsize,RBsize
integer::ix,l
real(8)::condL,condR

S00=czero ! overlap matrix
do i=1,nm
    S00(i,i)=1.0d0
enddo
if (NBC>0) then ! BC correction
  NL=nm ! left contact block size
  NR=nm ! right contact block size
  LBsize=NL*NBC
  RBsize=NR*NBC
  allocate(M(nm,nm,nx))
  allocate(M1i(nm,nm,nx))
  allocate(Mi1(nm,nm,nx))
  allocate(S(nm,nm,nx))
  allocate(S1i(nm,nm,nx))
  allocate(Si1(nm,nm,nx))
  allocate(LL(nm,nm,nx))
  allocate(LG(nm,nm,nx))
  LL=czero
  LG=czero
  allocate(V00 (LBsize,LBsize))
  allocate(V01 (LBsize,LBsize))
  allocate(V10 (LBsize,LBsize))
  allocate(M00 (LBsize,LBsize))
  allocate(M01 (LBsize,LBsize))
  allocate(M10 (LBsize,LBsize))
  allocate(PR00(LBsize,LBsize))
  allocate(PR01(LBsize,LBsize))
  allocate(PR10(LBsize,LBsize))
  allocate(PG00(LBsize,LBsize))
  allocate(PG01(LBsize,LBsize))
  allocate(PG10(LBsize,LBsize))
  allocate(PL00(LBsize,LBsize))
  allocate(PL01(LBsize,LBsize))
  allocate(PL10(LBsize,LBsize))
  allocate(LG00(LBsize,LBsize))
  allocate(LG01(LBsize,LBsize))
  allocate(LG10(LBsize,LBsize))
  allocate(LL00(LBsize,LBsize))
  allocate(LL01(LBsize,LBsize))
  allocate(LL10(LBsize,LBsize))
  allocate(dM11(LBsize,LBsize))
  allocate(xR11(LBsize,LBsize))
  allocate(dV11(LBsize,LBsize))
  allocate(dLL11(LBsize,LBsize))
  allocate(dLG11(LBsize,LBsize))
  !
  allocate(VNN  (RBsize,RBsize))
  allocate(VNN1 (RBsize,RBsize))
  allocate(Vn1n (RBsize,RBsize))
  allocate(Mnn  (RBsize,RBsize))
  allocate(Mnn1 (RBsize,RBsize))
  allocate(Mn1n (RBsize,RBsize))
  allocate(PRnn (RBsize,RBsize))
  allocate(PRnn1(RBsize,RBsize))
  allocate(PRn1n(RBsize,RBsize))
  allocate(PGnn (RBsize,RBsize))
  allocate(PGnn1(RBsize,RBsize))
  allocate(PGn1n(RBsize,RBsize))
  allocate(PLnn (RBsize,RBsize))
  allocate(PLnn1(RBsize,RBsize))
  allocate(PLn1n(RBsize,RBsize))
  allocate(LGnn (RBsize,RBsize))
  allocate(LGnn1(RBsize,RBsize))
  allocate(LGn1n(RBsize,RBsize))
  allocate(LLnn (RBsize,RBsize))
  allocate(LLnn1(RBsize,RBsize))
  allocate(LLn1n(RBsize,RBsize))
  allocate(dMnn (RBsize,RBsize))
  allocate(xRnn (RBsize,RBsize))
  allocate(dLLnn(RBsize,RBsize))
  allocate(dLGnn(RBsize,RBsize))
  allocate(dVnn(RBsize,RBsize))
  ! left
  call get_OBC_blocks_for_W(NL,Vii(:,:,1),V1i(:,:,1),PR(:,:,1),PR1i(:,:,1),&
      PL(:,:,1),PL1i(:,:,1),PG(:,:,1),PG(:,:,1),NBC,&
      V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
      LL00,LL01,LL10,LG00,LG01,LG10)
  ! right   
  call get_OBC_blocks_for_W(NR,Vii(:,:,nx),transpose(conjg(V1i(:,:,nx))),PR(:,:,nx),&
      transpose(PR1i(:,:,nx)),PL(:,:,nx),-transpose(conjg(PL1i(:,:,nx))),&
      PG(:,:,nx),-transpose(conjg(PG1i(:,:,nx))),NBC,&
      VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
      LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
  !
  !! S = V P^r
  ! Si,i = Vi,i Pi,i + Vi,i+1 Pi+1,i + Vi,i-1 Pi-1,i
  do ix=1,nx
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,ix),nm,PR(:,:,ix),nm,czero,S(:,:,ix),nm)
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
    if (ix==1) then
      call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
    else
      call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,ix-1),nm,PR1i(:,:,ix-1),nm,cone,S(:,:,ix),nm) 
    endif
  
  enddo
  !! Correct first and last block to account for elements in the contacts
  S(:,:,1)=S(:,:,1) + matmul(V10,PR01)
  S(:,:,nx)=S(:,:,nx) + matmul(VNN1,PRN1N)
  !
  do i=1,nx
    ! Si+1,i = Vi+1,i+1 Pi+1,i + Vi+1,i Pi,i
    if (i==nx) then
      call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PR1i(:,:,i),nm,cone,S1i(:,:,i),nm) 
    else    
      call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PR1i(:,:,i),nm,cone,S1i(:,:,i),nm) 
    endif
    call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i),nm,cone,S1i(:,:,i),nm) 
    ! Si,i+1 = Vi,i Pi,i+1 + Vi,i+1 Pi+1,i+1
    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PR1i(:,:,i),nm,cone,Si1(:,:,i),nm)
    if (i==nx) then 
      call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i),nm,cone,Si1(:,:,i),nm) 
    else
      call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i+1),nm,cone,Si1(:,:,i),nm) 
    endif
  enddo
  M = -S
  M1i=-S1i
  Mi1=-Si1
  do ix=1,nx
    do i=1,nm
       M(i,i,ix) = dcmplx(1.0d0, 1.0d-3) + M(i,i,ix) 
    enddo
  enddo
  deallocate(S,S1i,Si1)
  !! LL=V P^l V'
  !! LLi,i
  do i=1,nx
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL(:,:,min(i+1,nx)),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PL(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LL(:,:,i),nm)
    !
    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
  enddo
  !!Correct first and last block to account for elements in the contacts
  LL(:,:,1)=LL(:,:,1) + matmul(matmul(V10,PL00),V01) + &
    matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
  LL(:,:,nx)=LL(:,:,nx) + &
    matmul(matmul(VNN,PLNN1),VN1N) + matmul(matmul(VNN1,PLN1N),VNN) + &
    matmul(matmul(VNN1,PLNN),VN1N)
  !
  !! LG=V P^g V'    
  !! LGi,i
  do i=1,nx
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG(:,:,min(i+1,nx)),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PG(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LG(:,:,i),nm)
    !
    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
  enddo
  !!Correct first and last block to account for elements in the contacts
  LG(:,:,1)=LG(:,:,1) + matmul(matmul(V10,PG00),V01) + &
    matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
  LG(:,:,nx)=LG(:,:,nx) + &
    matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + matmul(matmul(VNN1,PGNN),VN1N)
    
  ! WR/WL/WG OBC Left
  call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
  !! WR/WL/WG OBC right
  call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
  allocate(VV(nm,nm,nx))
  VV = Vii
  if (condL<1.0d-6) then   
      !
      call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
     !
      M(:,:,1)=M(:,:,1) - dM11
      VV(:,:,1)=VV(:,:,1) - dV11    
      LL(:,:,1)=LL(:,:,1) + dLL11
      LG(:,:,1)=LG(:,:,1) + dLG11    
  endif
  if (condR<1.0d-6) then    
      !
      call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
      !
      M(:,:,nx)=M(:,:,nx) - dMNN
      VV(:,:,nx)=VV(:,:,nx)- dVNN
      LL(:,:,nx)=LL(:,:,nx) + dLLNN
      LG(:,:,nx)=LG(:,:,nx) + dLGNN    
  endif
  
  deallocate(V00,V01,V10)
  deallocate(M00,M01,M10)
  deallocate(PR00,PR01,PR10)
  deallocate(PG00,PG01,PG10)
  deallocate(PL00,PL01,PL10)
  deallocate(LG00,LG01,LG10)
  deallocate(LL00,LL01,LL10)
  deallocate(VNN,VNN1,Vn1n)
  deallocate(Mnn,Mnn1,Mn1n)
  deallocate(PRnn,PRnn1,PRn1n)
  deallocate(PGnn,PGnn1,PGn1n)
  deallocate(PLnn,PLnn1,PLn1n)
  deallocate(LGnn,LGnn1,LGn1n)
  deallocate(LLnn,LLnn1,LLn1n)
  deallocate(dM11,xR11,dLL11,dLG11,dV11)
  deallocate(dMnn,xRnn,dLLnn,dLGnn,dVnn)

endif
!
allocate(xlr(nm,nm,nx)) ! left-connected xR
allocate(wlr(nm,nm,nx)) ! left-connected Wr
allocate(wln(nm,nm,nx)) ! left-connected W<
allocate(wlp(nm,nm,nx)) ! left-connected W>
allocate(xR(nm,nm,nx))  ! fully connected xR
WR=czero
WL=czero
WG=czero
!
! Start RGF



!! on the left contact        
!H00(:,:)=M(:,:,1)
!H10(:,:)=M1i(:,:,1)
!!
!A=H00
!!                
!call invert(A,nm)
!xlr(:,:,1)=A(:,:)
!!
!call zgemm('n','n',nm,nm,nm,cone,LL(:,:,1),nm,S00,nm,czero,B,nm)
!sig=B
!call zgemm('n','n',nm,nm,nm,cone,A,nm,sig,nm,czero,B,nm) 
!call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm) 
!wln(:,:,1)=C(:,:)
!!
!call zgemm('n','n',nm,nm,nm,cone,LG(:,:,1),nm,S00,nm,czero,B,nm)
!sig=B
!call zgemm('n','n',nm,nm,nm,cone,A,nm,sig,nm,czero,B,nm) 
!call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm) 
!wlp(:,:,1)=C(:,:)
!Do l=2,nx-1                
!    H00(:,:)=M(:,:,l)
!    H10(:,:)=M1i(:,:,l-1)
!    H01(:,:)=Mi1(:,:,l-1)
!    call zgemm('n','n',nm,nm,nm,cone,H10,nm,xlr(:,:,max(l-1,1)),nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,H01,nm,czero,C,nm)
!    A=H00-C   
!    !
!    call invert(A,nm)
!    xlr(:,:,l)=A(:,:)
!    !
!    sig=wln(:,:,max(l-1,1))
!    call zgemm('n','n',nm,nm,nm,cone,H10,nm,sig,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,H01,nm,czero,C,nm)     
!    call zgemm('n','n',nm,nm,nm,cone,LL(:,:,l),nm,S00,nm,czero,B,nm)
!    C(:,:)=C(:,:)+B(:,:)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
!    call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wn,nm)
!    wln(:,:,l)=Wn(:,:)
!    !
!    sig=wlp(:,:,max(l-1,1))
!    call zgemm('n','n',nm,nm,nm,cone,H10,nm,sig,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,H01,nm,czero,C,nm)     
!    call zgemm('n','n',nm,nm,nm,cone,LG(:,:,l),nm,S00,nm,czero,B,nm)
!    C(:,:)=C(:,:)+B(:,:)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
!    call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wp,nm)
!    wlp(:,:,l)=Wp(:,:)
!enddo
!! on the right contact
!H00(:,:)=M(:,:,nx)
!H10(:,:)=M1i(:,:,nx-1)
!H01(:,:)=Mi1(:,:,nx-1)    
!call zgemm('n','n',nm,nm,nm,cone,H10,nm,xlr(:,:,nx-1),nm,czero,B,nm) 
!call zgemm('n','n',nm,nm,nm,cone,B,nm,H01,nm,czero,C,nm)
!G00=H00-C   
!!
!call invert(G00,nm)
!!
!call zgemm('n','n',nm,nm,nm,cone,G00,nm,Vii(:,:,nx),nm,czero,WR(:,:,nx),nm)
!!
!sig=wln(:,:,nx-1)
!call zgemm('n','n',nm,nm,nm,cone,H10,nm,sig,nm,czero,B,nm) 
!call zgemm('n','n',nm,nm,nm,cone,B,nm,H01,nm,czero,C,nm)  ! C=H10 Gl< H01
!call zgemm('n','n',nm,nm,nm,cone,LL(:,:,nx),nm,S00,nm,czero,B,nm)
!! B=Sig< S00
!sig(:,:)=C(:,:)+B(:,:)
!call zgemm('n','n',nm,nm,nm,cone,G00,nm,sig,nm,czero,B,nm) 
!call zgemm('n','c',nm,nm,nm,cone,B,nm,G00,nm,czero,Wn,nm) 
!! G<00 = G00 sig< G00'
!WL(:,:,nx)=Wn(:,:)
!WG(:,:,nx)=Wn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
!Wp(:,:)=WG(:,:,nx)
!!-------------------------
!do l=nx-1,1,-1
!    H10(:,:)=M1i(:,:,l)
!    H01(:,:)=Mi1(:,:,l)
!    A=Wn
!    call zgemm('n','c',nm,nm,nm,cone,H10,nm,xlr(:,:,l),nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,B,nm,czero,C,nm) 
!    A=wln(:,:,l)          
!    call zgemm('n','n',nm,nm,nm,cone,H10,nm,A,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,G00,nm,B,nm,czero,A,nm)
!    B=C+A
!    call zgemm('n','n',nm,nm,nm,cone,H01,nm,B,nm,czero,A,nm)      !!! G<_i+1,i                     
!    D(:,:)= xlr(:,:,l)
!    call zgemm('n','n',nm,nm,nm,cone,D,nm,H01,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,G00,nm,czero,GN0,nm)      !!! G_i,i+1
!    !
!    call zgemm('n','n',nm,nm,nm,cone,GN0,nm,V1i(:,:,l),nm,cone,WR(:,:,l),nm)    ! G_i,i+1 V_i+1,i -> WR_i,i    
!    !
!    call zgemm('n','n',nm,nm,nm,cone,GN0,nm,H10,nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,D,nm,czero,C,nm)     
!    G00(:,:)=xlr(:,:,l)+C(:,:)                                       !!! G_i,i
!    xR(:,:,l)=G00(:,:)
!    !
!    call zgemm('n','n',nm,nm,nm,cone,G00,nm,Vii(:,:,l),nm,cone,WR(:,:,l),nm)     ! G_i,i V_i,i  -> WR_i,i   
!    !
!    call zgemm('n','n',nm,nm,nm,cone,xlr(:,:,l+1),nm,H10,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,G00,nm,czero,G0N,nm)      !!! G_i+1,i
!    !
!    call zgemm('n','c',nm,nm,nm,cone,G0N,nm,V1i(:,:,l),nm,cone,WR(:,:,l+1),nm)  ! G_i+1,i V_i,i+1  -> WR_i+1,i+1
!    !-------------------------
!    A(:,:)=Wn(:,:)     
!    call zgemm('n','n',nm,nm,nm,cone,D,nm,H01,nm,czero,B,nm)  
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm)     
!    call zgemm('n','n',nm,nm,nm,cone,C,nm,H10,nm,czero,A,nm)
!    call zgemm('n','c',nm,nm,nm,cone,A,nm,D,nm,czero,C,nm)
!    Wn(:,:)= wln(:,:,l) + C(:,:)
!    A(:,:)=wln(:,:,l)
!    call zgemm('n','n',nm,nm,nm,cone,GN0,nm,H10,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm)         
!    Wn(:,:)= Wn(:,:)+C(:,:)!                     			 
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,H01,nm,czero,B,nm) 
!    call zgemm('n','c',nm,nm,nm,cone,B,nm,GN0,nm,czero,C,nm)          
!    Wn(:,:)= Wn(:,:)+C(:,:)!     					 !!! G<_i,i
!    !-------------------------
!    A(:,:)=Wp(:,:)
!    call zgemm('n','c',nm,nm,nm,cone,D,nm,H10,nm,czero,B,nm)  
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm)     
!    !
!    call zgemm('n','n',nm,nm,nm,cone,C,nm,H10,nm,czero,A,nm)
!    call zgemm('n','c',nm,nm,nm,cone,A,nm,D,nm,czero,C,nm)
!    !
!    Wp(:,:)= wlp(:,:,l) + C(:,:)
!    A(:,:)=wlp(:,:,l)
!    call zgemm('n','n',nm,nm,nm,cone,GN0,nm,H10,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,A,nm,czero,C,nm)     
!    !
!    Wp(:,:)= Wp(:,:)+C(:,:)!                     			 
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,H01,nm,czero,B,nm) 
!    call zgemm('n','n',nm,nm,nm,cone,B,nm,G0N,nm,czero,C,nm)     
!    Wp(:,:)= Wp(:,:)+C(:,:)!     					 !!! G>_i,i
!    !-------------------------    
!    WL(:,:,l)=Wn(:,:)
!    WG(:,:,l)=Wp(:,:)
!enddo
        




! first pass, from right to left
A=M(:,:,nx)
call invert(A,nm)
xlr(:,:,nx)=A
call zgemm('n','n',nm,nm,nm,cone,A,nm,VV(:,:,nx),nm,czero,Wlr(:,:,nx),nm) 
!
call zgemm('n','n',nm,nm,nm,cone,A,nm,LL(:,:,nx),nm,czero,B,nm) 
call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wln(:,:,nx),nm) 
!
call zgemm('n','n',nm,nm,nm,cone,A,nm,LG(:,:,nx),nm,czero,B,nm) 
call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wlp(:,:,nx),nm) 
!
do ix=nx-1,1,-1
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,Xlr(:,:,ix+1),nm,czero,B,nm) ! B -> MxR
  call zgemm('n','c',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,C,nm)
  A = M(:,:,ix) - C
  call invert(A,nm)
  Xlr(:,:,ix)=A     
  call zgemm('n','n',nm,nm,nm,cone,B,nm,V1i(:,:,ix),nm,czero,C,nm) 
  D=VV(:,:,nx) - C
  call zgemm('n','n',nm,nm,nm,cone,A,nm,D,nm,czero,Wlr(:,:,ix),nm) 
  !
  !AL=MxR*LL10
  !AG=MxR*LG10
  !
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,Wln(:,:,ix+1),nm,czero,B,nm) 
  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)       
  C=LL(:,:,ix)+sig !-(AL-transpose(conjg(AL)))
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wn,nm)     
  Wln(:,:,ix)=Wn
  !
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,Wlp(:,:,ix+1),nm,czero,B,nm) 
  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)       
  C=LG(:,:,ix)+sig !-(AG-transpose(conjg(AG)))
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wp,nm)     
  Wlp(:,:,ix)=Wp    
enddo
! second pass, from left to right
WR(:,:,1)=Wlr(:,:,1)
call zgemm('n','n',nm,nm,nm,cone,WR(:,:,1),nm,Mi1(:,:,1),nm,czero,B,nm)
call zgemm('n','c',nm,nm,nm,cone,B,nm,xlr(:,:,2),nm,czero,B,nm) 
WRi1(:,:,1)=transpose(conjg(V1i(:,:,1)))-B  
xR(:,:,1)=xlr(:,:,1)
WL(:,:,1)=Wln(:,:,1)
WG(:,:,1)=Wlp(:,:,1)
do ix=2,nx
  call zgemm('n','n',nm,nm,nm,cone,xlr(:,:,ix),nm,M1i(:,:,ix-1),nm,czero,B,nm)
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WRi1(:,:,ix-1),nm,czero,B,nm) 
  WR(:,:,ix)=Wlr(:,:,ix)  - B
  call zgemm('n','n',nm,nm,nm,cone,xlr(:,:,ix),nm,M1i(:,:,ix-1),nm,czero,B,nm) ! B -> xlr * M10
  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix-1),nm,czero,A,nm) 
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix-1),nm,xlr(:,:,ix),nm,czero,C,nm)   ! C -> M01 * xlr
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)   
  xR(:,:,ix)=xlr(:,:,ix)  + D
  !   
  ! AL=xlr LL10 XR' M01 xlr
  ! BL=xlr M10 XR M01 wln
  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix),nm,czero,D,nm)   
  call zgemm('n','n',nm,nm,nm,cone,D,nm,Mi1(:,:,ix-1),nm,czero,A,nm)   
  call zgemm('n','n',nm,nm,nm,cone,A,nm,wln(:,:,ix),nm,czero,BL,nm)   
  WL(:,:,ix)=wln(:,:,ix) + (BL-transpose(conjg(BL)))
  !
  ! AG=xlr LG10 XR' M01 xlr
  ! BG=xlr M10 XR M01 wlp    
  call zgemm('n','n',nm,nm,nm,cone,A,nm,wlp(:,:,ix),nm,czero,BL,nm)   
  !
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WL(:,:,ix-1),nm,czero,A,nm) 
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)   
  WL(:,:,ix)=WL(:,:,ix) + D
  !
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WG(:,:,ix-1),nm,czero,A,nm) 
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)   
  !       
  WG(:,:,ix)=wlp(:,:,ix)  + (BL-transpose(conjg(BL))) + D
enddo

deallocate(M,M1i,Mi1,LL,LG)
deallocate(VV)
deallocate(wln,wlp,xlr,Xr,wlr)
end subroutine green_rgf_calc_w


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
    real(8),intent(inout) :: tr,tre
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
    real(8),intent(inout) :: tr,tre
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
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
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
subroutine write_current_spectrum_summed_over_kz(dataset,i,cur,nen,en,length,NB,Lx,nk)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: cur(:,:,:,:,:)
  integer, intent(in)::i,nen,length,NB,nk
  real(8), intent(in)::Lx,en(nen)
  integer:: ie,j,ib,jb,ik
  real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
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
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
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
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
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
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
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
