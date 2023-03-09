!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf

implicit none 

private

public :: green_rgf_cms, green_rgf_solve_gw_1d

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
real(8), parameter :: m0=9.109d-31 ! kg
real(8), parameter :: eps0=8.854d-12 ! C/V/m 
real(8), parameter :: c0=2.998d8 ! m/s
real(8), parameter :: e0=1.6022d-19 ! C
REAL(8), PARAMETER :: pi = 3.14159265359d0

CONTAINS

subroutine green_rgf_solve_gw_1d(alpha_mix,niter,NB,NS,nm,nx,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i)
  implicit none
  integer,intent(in)::nm,nx,nen,niter,NB,NS
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix
  complex(8),intent(in),dimension(nm,nm,nx)::Hii,H1i,Vii,V1i
  ! --------
  complex(8),allocatable,dimension(:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:)::P_lesser,P_greater,P_retarded
  complex(8),allocatable,dimension(:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:)::W_lesser,W_greater,W_retarded
  complex(8),allocatable,dimension(:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
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
  allocate(P_lesser_1i(nm,nm,nx,nen))
  allocate(P_greater_1i(nm,nm,nx,nen))
  allocate(P_retarded_1i(nm,nm,nx,nen))
  allocate(W_lesser(nm,nm,nx,nen))
  allocate(W_greater(nm,nm,nx,nen))
  allocate(W_retarded(nm,nm,nx,nen))
  allocate(W_lesser_i1(nm,nm,nx,nen))
  allocate(W_greater_i1(nm,nm,nx,nen))
  allocate(W_retarded_i1(nm,nm,nx,nen))
  allocate(Sii(nm,nm,nx))
  Sii=0.0d0
  do i=1,nm
    Sii(i,i,:)=1.0d0
  enddo
  do iter=0,niter
    !
    print *, 'calc G'  
    do ie=1,nen
      call green_RGF_CMS(TEMP,nm,nx,En(ie),mu,Sii,Hii,H1i,sigma_lesser_gw(:,:,:,ie),sigma_greater_gw(:,:,:,ie),&
                      & sigma_r_gw(:,:,:,ie),g_lesser(:,:,:,ie),g_greater(:,:,:,ie),g_r(:,:,:,ie),tr(ie),tre(ie),cur(:,:,:,ie),g_r_i1(:,:,:,ie)) 
    enddo
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    call write_spectrum('gw_ldos',iter,g_r,nen,En,nx,NB,NS,Lx,(/1.0,-2.0/))  
    call write_spectrum('gw_ndos',iter,g_lesser,nen,En,nx,NB,NS,Lx,(/1.0,1.0/))       
    call write_spectrum('gw_pdos',iter,g_greater,nen,En,nx,NB,NS,Lx,(/1.0,-1.0/))       
    !        
    print *, 'calc P'    
    nopmax=nen/2-10  
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	  
    P_lesser(:,:,:,:) = dcmplx(0.0d0,0.0d0)
    P_greater(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
    P_retarded(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
    P_lesser_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)
    P_greater_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
    P_retarded_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i) shared(nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx)
    l=1
    h=nm
    !$omp do 
    do io=1,nx*nm
      ix=(io-1)/nm+1
      i=io-(ix-1)*nm
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
    !$omp end do 
    !$omp end parallel
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum('gw_PR',iter,P_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    !  call write_spectrum('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    !  call write_spectrum('gw_PG',iter,P_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    !
    print *, 'calc W'
    do ie=1,nen
      call green_rgf_calc_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie),w_lesser_i1(:,:,:,ie),w_greater_i1(:,:,:,ie),w_retarded_i1(:,:,:,ie))
    enddo
    call write_spectrum('gw_WR',iter,W_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    call write_spectrum('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    call write_spectrum('gw_WG',iter,W_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0,1.0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i) shared(nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx)
    l=1
    h=nm
    !$omp do
    do ix=1,nx
    do i=1,nm
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
    call write_spectrum('gw_SigR',iter,Sigma_r_gw,nen,En,nx,NB,NS,Lx,(/1.0,1.0/))
    !call write_spectrum('gw_SigL',iter,Sigma_lesser_gw,nen,En,nx,NB,NS,Lx,(/1.0,1.0/))
    !call write_spectrum('gw_SigG',iter,Sigma_greater_gw,nen,En,nx,NB,NS,Lx,(/1.0,1.0/))
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)
  deallocate(Sii)         
end subroutine green_rgf_solve_gw_1d


subroutine green_rgf_calc_w(nm,nx,Vii,V1i,PL,PG,PR,PL1i,PG1i,PR1i,WL,WG,WR,WLi1,WGi1,WRi1)
implicit none
integer,intent(in)::nm,nx
complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,PL,PG,PR
complex(8),intent(in),dimension(nm,nm,nx) :: PL1i,PG1i,PR1i
complex(8),intent(inout),dimension(nm,nm,nx) :: WL,WG,WR,WLi1,WGi1,WRi1
! --------
COMPLEX(8)::A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),Wn(nm,nm),Wp(nm,nm),sig(nm,nm),AL(nm,nm),BL(nm,nm)
complex(8),allocatable,dimension(:,:,:)::S,M,LL,LG,VV,M1i,Mi1,S1i,Si1
complex(8),allocatable,dimension(:,:,:)::xlr,wlr,wln,wlp,xR
complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
    PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
    MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
integer::i,NL,NR,NT,LBsize,RBsize,NBC
integer::ix
real(8)::condL,condR
NBC=2
NL=nm ! left contact block size
NR=nm ! right contact block size
NT=nm*nx! total size
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
! Correct first and last block to account for elements in the contacts
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
     M(i,i,ix) = 1.0d0 + M(i,i,ix)
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
!Correct first and last block to account for elements in the contacts
LL(:,:,1)=LL(:,:,1) + matmul(matmul(V10,PL00),V01) + &
  matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
!  
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
!Correct first and last block to account for elements in the contacts
LG(:,:,1)=LG(:,:,1) + matmul(matmul(V10,PG00),V01) + &
  matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
LG(:,:,nx)=LG(:,:,nx) + &
  matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + matmul(matmul(VNN1,PGNN),VN1N)
  
! WR/WL/WG OBC Left
call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
! WR/WL/WG OBC right
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

allocate(xlr(nm,nm,nx)) ! left-connected xR
allocate(wlr(nm,nm,nx)) ! left-connected Wr
allocate(wln(nm,nm,nx)) ! left-connected W<
allocate(wlp(nm,nm,nx)) ! left-connected W>
allocate(xR(nm,nm,nx)) ! fully connected xR

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
  
  !
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WL(:,:,ix-1),nm,czero,A,nm) 
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)   
  WL(:,:,ix)=wln(:,:,ix)  + D + (BL-transpose(conjg(BL)))
  ! BG=xlr M10 XR M01 wlp  
  call zgemm('n','n',nm,nm,nm,cone,A,nm,wlp(:,:,ix),nm,czero,BL,nm)   
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WG(:,:,ix-1),nm,czero,A,nm) 
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)   
  !   
  ! AG=xlr LG10 XR' M01 xlr
  
  
  WG(:,:,ix)=wlp(:,:,ix)  + D
enddo
!
deallocate(M,M1i,Mi1,LL,LG,VV)
deallocate(wln,wlp,wlr,xlr,Xr)
end subroutine green_rgf_calc_w




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




!!!! RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_CMS(TEMP,nm,nx,E,mu,Sii,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,Gi1) 
    implicit none    
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of GFs
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Sii,Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: cur,ldos,ndens,pdens
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: Gi1
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
    pdens(:,:,nx)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
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
        cur(:,:,l)=2.0d0*dble(A(:,:))     
        !-------------------------
        !     A=Gn
        !     call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        !     call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)  
        !     A=Gln(:,:,l)
        !     call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        !     call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) 
        !     B=C+A
        !     call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)    !!! G<_i,i+1
        !     cur(:,:,l)=cur(:,:,l)-A(:,:)
        !     cur(:,:,l)=dble(cur(:,:,l))
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(Gi1)) then 
          Gi1(:,:,l)=GN0
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
implicit none
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
    implicit none      
    integer :: info,lda,lwork,nn      
    integer, dimension(:), allocatable :: ipiv
    complex(8), dimension(nn,nn),intent(inout) :: A
    complex(8), dimension(:), allocatable :: work
    allocate(work(nn*nn))
    allocate(ipiv(nn))
    call zgetrf(nn,nn,A,nn,ipiv,info)
    if (info.ne.0) then
      print*,'SEVERE warning: zgetrf failed, info=',info
      call abort
    endif
    call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
    if (info.ne.0) then
      print*,'SEVERE warning: zgetri failed, info=',info
      call abort
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
        write(11,'(4E18.4)') j*Lx*NS+(k-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
      enddo
    end do
    write(11,*)    
end do
close(11)
end subroutine write_spectrum

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
implicit none
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


! calculate matrix blocks for the Open Boundary Condition of W
subroutine get_OBC_blocks_for_W(n,v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01,NBC,&
    V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
    LL00,LL01,LL10,LG00,LG01,LG10)
implicit none
integer,intent(in)::n,NBC
complex(8),intent(in),dimension(n,n)::v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01
complex(8),intent(out),dimension(n*NBC,n*NBC)::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
    LL00,LL01,LL10,LG00,LG01,LG10
complex(8),dimension(n*NBC,n*NBC)::II
!
select case (NBC)
  !
  case(1)
    !
    V00=v_00
    V01=v_01
    V10=transpose(conjg(V01))
    !
    PR00=pR_00
    PR01=pR_01
    PR10=transpose(PR01)
    !
    PL00=pL_00;
    PL01=pL_01;
    PL10=-transpose(conjg(PL01))
    !
    PG00=pG_00
    PG01=pG_01
    PG10=-transpose(conjg(PG01))
        
  case(2)
    !
    V00(1:n,1:n)=v_00 
    V00(1:n,n+1:2*n)=v_01
    V00(n+1:2*n,1:n)=transpose(conjg(v_01))
    V00(n+1:2*n,n+1:2*n)= v_00
    V01=czero
    V01(n+1:2*n,1:n)=v_01
    V10=transpose(conjg(V01))
    !
    PR00(1:n,1:n)=pR_00 
    PR00(1:n,n+1:2*n)=pR_01
    PR00(n+1:2*n,1:n)=transpose(pR_01)
    PR00(n+1:2*n,n+1:2*n)= pR_00
    PR01=czero
    PR01(n+1:2*n,1:n)=pR_01
    PR10=transpose(PR01)
    !
    PG00(1:n,1:n)=pG_00 
    PG00(1:n,n+1:2*n)=pG_01
    PG00(n+1:2*n,1:n)=-transpose(conjg(pG_01))
    PG00(n+1:2*n,n+1:2*n)= pG_00
    PG01=czero
    PG01(n+1:2*n,1:n)=pG_01
    PG10=-transpose(conjg(PG01))
    !
    PL00(1:n,1:n)=pL_00 
    PL00(1:n,n+1:2*n)=pL_01
    PL00(n+1:2*n,1:n)=-transpose(conjg(pL_01))
    PL00(n+1:2*n,n+1:2*n)= pL_00
    PL01=czero
    PL01(n+1:2*n,1:n)=pL_01
    PL10=-transpose(conjg(PL01))    

!    case(3)

!        V00=[[v_00 v_01 sparse(N,N)];[v_01' v_00 v_01];[sparse(N,N),v_01' v_00]];
!        V01=[sparse(2*N,3*N);[v_01 sparse(N,2*N)]];
!        V10=V01';
    
!        PR00=[[pR_00 pR_01 sparse(N,N)];[conj(pR_01') pR_00 pR_01];[sparse(N,N),conj(pR_01') pR_00]];
!        PR01=[sparse(2*N,3*N);[pR_01 sparse(N,2*N)]];
!        PR10=conj(PR01');
    
!        PL00=[[pL_00 pL_01 sparse(N,N)];[-pL_01' pL_00 pL_01];[sparse(N,N),-pL_01' pL_00]];
!        PL01=[sparse(2*N,3*N);[pL_01 sparse(N,2*N)]];
!        PL10=-PL01';
    
!        PG00=[[pG_00 pG_01 sparse(N,N)];[-pG_01' pG_00 pG_01];[sparse(N,N),-pG_01' pG_00]];
!        PG01=[sparse(2*N,3*N);[pG_01 sparse(N,2*N)]];
!        PG10=-PG01';
end select

call identity(II,NBC*N)
M00=II*dcmplx(1.0d0,1d-10)-matmul(V10,PR01)
M00=M00-matmul(V00,PR00)-matmul(V01,PR10)
M01=-matmul(V00,PR01)-matmul(V01,PR00)
M10=-matmul(V10,PR00)-matmul(V00,PR10)
!
LL00=matmul(matmul(V10,PL00),V01)+matmul(matmul(V10,PL01),V00)
LL00=LL00+matmul(matmul(V00,PL10),V01)+matmul(matmul(V00,PL00),V00)
LL00=LL00+matmul(matmul(V00,PL01),V10)
LL00=LL00+matmul(matmul(V01,PL10),V00)+matmul(matmul(V01,PL00),V10)
!
LL01=matmul(matmul(V10,PL01),V01)+matmul(matmul(V00,PL00),V01)
LL01=LL01+matmul(matmul(V00,PL01),V00)+matmul(matmul(V01,PL10),V01)
LL01=LL01+matmul(matmul(V01,PL00),V00)
!
LL10=-transpose(conjg(LL01))
!
LG00=matmul(matmul(V10,PG00),V01)+matmul(matmul(V10,PG01),V00)
LG00=LG00+matmul(matmul(V00,PG10),V01)+matmul(matmul(V00,PG00),V00)
LG00=LG00+matmul(matmul(V00,PG01),V10)
LG00=LG00+matmul(matmul(V01,PG10),V00)+matmul(matmul(V01,PG00),V10)
!
LG01=matmul(matmul(V10,PG01),V01)
LG01=LG01+matmul(matmul(V00,PG00),V01)+matmul(matmul(V00,PG01),V00)
LG01=LG01+matmul(matmul(V01,PG10),V01)+matmul(matmul(V01,PG00),V00)
!
LG10=-transpose(conjg(LG01))  
end subroutine get_OBC_blocks_for_W


! calculate corrections to the L matrix blocks for the Open Boundary Condition
subroutine get_dL_OBC_for_W(nm,xR,LL00,LL01,LG00,LG01,M10,typ, dLL11,dLG11)
implicit none
integer,intent(in)::nm
character(len=*),intent(in)::typ
complex(8),intent(in),dimension(nm,nm)::xR,LL00,LL01,LG00,LG01,M10
complex(8),intent(out),dimension(nm,nm)::dLL11,dLG11
! -----
complex(8),dimension(nm,nm)::AL,AG,FL,FG,A,V,iV,yL_NN,wL_NN,yG_NN,wG_NN,tmp1,tmp2
complex(8),dimension(nm)::E
integer::i,j
!!!! AL=M10*xR*LL01;
!!!! AG=M10*xR*LG01;
call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LL01,nm,czero,AL,nm)
call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LG01,nm,czero,AG,nm)
!!!! FL=xR*(LL00-(AL-AL'))*xR';
!!!! FG=xR*(LG00-(AG-AG'))*xR';
call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LL00-(AL-transpose(conjg(AL)))),nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FL,nm)
call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LG00-(AG-transpose(conjg(AG)))),nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FG,nm)
!
call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,V,nm)  
do i=1,nm
  V(i,i)=V(i,i)+dcmplx(0.0d0,1.0d-4)  ! 1i*1e-4 added to stabilize matrix
enddo
E=eigv(nm,V)
iV=V
call invert(iV,nm)
!lesser component
call zgemm('n','n',nm,nm,nm,cone,iV,nm,FL,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yL_NN,nm)
yL_NN=yL_NN/(1.0d0 - sum(E*conjg(E)))
call zgemm('n','n',nm,nm,nm,cone,V,nm,yL_NN,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wL_NN,nm)
!refinement iteration
call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wL_NN,nm,czero,tmp2,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
wL_NN=FL+tmp2
!
call zgemm('n','n',nm,nm,nm,cone,M10,nm,wL_NN,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLL11,nm)
dLL11=dLL11-(AL-transpose(conjg(AL)))
!greater component
call zgemm('n','n',nm,nm,nm,cone,iV,nm,FG,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yG_NN,nm)
yG_NN=yG_NN/(1.0d0 - sum(E*conjg(E)))
call zgemm('n','n',nm,nm,nm,cone,V,nm,yG_NN,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wG_NN,nm)
!refinement iteration
call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wG_NN,nm,czero,tmp2,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
wG_NN=FG+tmp2
!
call zgemm('n','n',nm,nm,nm,cone,M10,nm,wG_NN,nm,czero,tmp1,nm)
call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLG11,nm)
dLG11=dLG11-(AG-transpose(conjg(AG)))
end subroutine get_dL_OBC_for_W


! a slightly modified version of sancho
subroutine surface_function(nm,M00,M01,M10,SF,cond)
implicit none
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
