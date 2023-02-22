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
      call green_rgf_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie))
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
  deallocate(Sii)         
end subroutine green_rgf_solve_gw_1d

! RGF for diagonal blocks of screened Coulomb GF W^r,<,>
! W^r = [I-VPr]^-1 V
! W^<> = W^r P<> W^a
! compare to usual RGF: 
!     V.Pr as Ham, and V.P<>.V' as Sig<>
!     Wr = Gr.V  
!     W<> = G<>
subroutine green_rgf_w(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,w_lesser,w_greater,w_retarded)
    implicit none
    integer,intent(in)::nm,nx
    complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,P_lesser,P_greater,P_retarded
    complex(8),intent(in),dimension(nm,nm,nx) :: P_lesser_1i,P_greater_1i,P_retarded_1i
    complex(8),intent(inout),dimension(nm,nm,nx) :: W_lesser,W_greater,W_retarded
    ! --------
    COMPLEX(8) :: S00(nm,nm),H00(nm,nm),H10(nm,nm),H01(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),GN0(nm,nm),G0N(nm,nm),Gn(nm,nm),Gp(nm,nm),G00(nm,nm),sigmal(nm,nm),sigmar(nm,nm),sig(nm,nm)       
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
    ! Sig = V.P.V
    call prepare_rgf_W(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,Hii,H1i,Hi1,sigma_lesser,sigma_greater)
    ! start RGF
    z=dcmplx(1.0d0, 1.0d-3)
    ! on the left contact        
    H00(:,:)=Hii(:,:,1)
    H10(:,:)=H1i(:,:,1)
    sigmal=0.0d0
    !
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser(:,:,1),nm,S00,nm,beta,B,nm)
    sig=B
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
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
        call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser(:,:,l),nm,S00,nm,beta,B,nm)
        C(:,:)=C(:,:)+B(:,:)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
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
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser(:,:,nx),nm,S00,nm,beta,B,nm)
    ! B=Sig< S00
    sig(:,:)=C(:,:)+B(:,:)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    w_lesser(:,:,nx)=Gn(:,:)
    w_greater(:,:,nx)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))    
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
        w_retarded(:,:,l)=G00(:,:)
        !
        call zgemm('n','n',nm,nm,nm,alpha,Gl(:,:,l+1),nm,H10,nm,beta,B,nm) 
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

subroutine prepare_rgf_W(nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,p_lesser_1i,p_greater_1i,p_retarded_1i,Hii,H1i,Hi1,sigma_lesser,sigma_greater)
implicit none
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
  ! Sig<> = V P<> V
  ! Sig<_i,i
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_lesser(:,:,min(i+1,nx)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_lesser_1i(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser_1i(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_lesser_1i(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_lesser(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_lesser(:,:,i),nm)
  !
  call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_lesser_1i(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_lesser(:,:,i),nm)
  ! Sig>_i,i
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_greater(:,:,min(i+1,nx)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,p_greater_1i(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater_1i(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,p_greater_1i(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_greater(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,sigma_greater(:,:,i),nm)
  !
  call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,p_greater_1i(:,:,max(i-1,1)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,sigma_greater(:,:,i),nm)
enddo
end subroutine prepare_rgf_W

! B = I-VP 
subroutine calc_OBC_blocks(P_retarded,B,V,NB,NS,nm_dev,nen,nop)
implicit none
integer,intent(in)::NB,NS,nm_dev,nen,nop
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::P_retarded
complex(8),intent(inout),dimension(nm_dev,nm_dev)::B,V
! -------
complex(8),allocatable::dB(:,:),S00(:,:),B00(:,:),B10(:,:),X00(:,:),XBB(:,:)
integer::i
allocate(S00(NB*NS,NB*NS))
allocate(dB(NB*NS,NB*NS))
allocate(B00(NB*NS,NB*NS))
allocate(B10(NB*NS,NB*NS))
allocate(X00(NB*NS,NB*NS))
allocate(XBB(NB*NS,NB*NS))
S00=0.0d0
do i=1,NB*NS
  S00(i,i)=1.0d0
enddo
! add OBC left
B00=-B(1:NB*NS,1:NB*NS)
B10=-B(NB*NS+1:2*NB*NS,1:NB*NS)
call sancho(NB*NS,1.0,S00,B00,B10,X00,XBB)
B(1:NB*NS,1:NB*NS)=B(1:NB*NS,1:NB*NS)-X00
! add OBC right
B00=-B(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev)
B10=-B(nm_dev-2*NB*NS+1:nm_dev-NB*NS,nm_dev-NB*NS+1:nm_dev)
call sancho(NB*NS,1.0,S00,B00,B10,X00,XBB)
B(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev)=B(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev)-X00
!   
call zgemm('n','n',nb*ns,nb*ns,nb*ns,-cone,V(NB*NS+1:2*NB*NS,1:NB*NS),nb*ns,P_retarded(1:NB*NS,NB*NS+1:2*NB*NS,nop),nb*ns,czero,dB,nb*ns)
B(1:NB*NS,1:NB*NS)=B(1:NB*NS,1:NB*NS)-dB
call zgemm('n','n',nb*ns,nb*ns,nb*ns,-cone,V(nm_dev-2*NB*NS+1:nm_dev-NB*NS,nm_dev-NB*NS+1:nm_dev),nb*ns,P_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-2*NB*NS+1:nm_dev-NB*NS,nop),nb*ns,czero,dB,nb*ns)
B(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev)=B(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev)-dB    
!
deallocate(B00,B10,X00,XBB)
deallocate(dB,S00)
end subroutine calc_OBC_blocks


! RGF for diagonal blocks of G^r,<,>
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
subroutine invert_band(A,nn,nb)
implicit none
integer,intent(in)::nn,nb
complex(8),intent(inout)::A(3*nb+1,nn)
complex(8),allocatable::work(:),B(:,:),X(:,:)
integer,allocatable::ipiv(:)
integer::info,lda,lwork,ldb,i
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
end subroutine invert_band


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

end module green_rgf
