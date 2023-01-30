!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 11/2022

module green

implicit none 

private

public :: green_calc_g, green_rgf_cms
public :: green_subspace_invert,green_solve_gw_1D

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)

CONTAINS


! RGF for G^r,<,>
subroutine green_RGF_CMS(TEMP,nm,nx,ny,EKZ,E,mul,mur,Sii,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,Jdens) 
  integer,intent(in)::nm,nx,ny
  real(8),intent(in) :: E,EKZ,mul,mur,temp
  COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph
  COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Sii,Hii,H1i
  COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: cur,ldos,ndens,pdens,Jdens
  real(8),intent(inout) :: tr,tre

  COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
  COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)!,Gl(nm,nm,nx),Gln(nm,nm,nx)!,Glp(nm,nm,nx)
  COMPLEX(8) :: z
  integer::i,j,k,l
  real(8)::tim
  COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:)
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
  !
  z=cmplx(E,0.0d-6)
  !
  allocate(Gl(nm,nm,nx))
  allocate(Gln(nm,nm,nx))
  !
  Gln=0.0D0
!!$  Glp=0.0D0
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
  sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMP))+B(:,:)
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
     H00(:,:)=Hii(:,:,l)+B(:,:)!sigma_r_ph(:,:,l)
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
     C(:,:)=C(:,:)+B(:,:)!sigma_lesser_ph(:,:,l)
     call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
     call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
     Gln(:,:,l)=Gn(:,:)
  enddo
  ! self energy on the right contact
  S00(:,:)=Sii(:,:,nx)
  call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,nx),nm,S00,nm,beta,B,nm)
  H00(:,:)=Hii(:,:,nx)+B(:,:)!sigma_r_ph(:,:,nx)
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
  call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
  call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,nx),nm,S00,nm,beta,B,nm)
  sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMP))+C(:,:)+B(:,:)
  call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
  call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
  ndens(:,:,nx)=Gn(:,:)
  pdens(:,:,nx)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))!Gp(:,:)!
  Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
  A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMP))
  call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
  A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMP))-1.0d0)
  call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
  tim=0.0d0
  do i=1,nm
      do j=i,i!1,nm
          tim=tim-dble(B(i,j)-C(i,j))
      enddo
  enddo
  tr=tim

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
     Jdens(:,:,l)=2.0d0*dble(A(:,:))
     cur(1,1,l)=2.0d0*dble(Hii(1,2,l+1)*Gn(2,1))        
     cur(2,2,l)=2.0d0*dble(Hii(2,1,l+1)*Gn(1,2))

!     A=Gn
!     call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
!     call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)  
!     A=Gln(:,:,l)
!     call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!     call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) 
!     B=C+A
!     call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)    
!     cur(:,:,l)=cur(:,:,l)-A(:,:)
!     cur(:,:,l)=dble(cur(:,:,l))

!-------------------------

     D(:,:)= Gl(:,:,l)
     call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
     call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
     call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
     call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
     G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
     ldos(:,:,l)=G00(:,:)!cmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))

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
!!$     A(:,:)=Gp(:,:)
!!$     
!!$     call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
!!$     call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!!$
!!$     call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
!!$     call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
!!$
!!$     Gp(:,:)= Glp(:,:,l) + C(:,:)
!!$     A(:,:)=Glp(:,:,l)
!!$     call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
!!$     call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!!$     
!!$     Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
!!$     call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!!$     call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
!!$     
!!$     Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G<_i,i
!-------------------------    
     ndens(:,:,l)=Gn(:,:)
     pdens(:,:,l)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))!Gp(:,:)
  enddo
  !
  Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
  A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMP))
  call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
  A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMP))-1.0d0)
  call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
  tim=0.0d0
  do i=1,nm
     tim=tim+dble(B(i,i)-C(i,i))
  enddo
  tre=tim
  deallocate(Gl)
  deallocate(Gln)
end subroutine green_RGF_CMS




! driver for iterating G -> P -> W -> Sig 
subroutine green_solve_gw_1D(niter,nm_dev,Lx,length,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
implicit none
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
integer :: iter,ie,nopmax
integer :: i,j,nm,nop,l,h,iop,ndiag
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
allocate(siglead(NB*NS,NB*NS,nen,2))
siglead=dcmplx(0.0d0,0.0d0)
allocate(B(nm_dev,nm_dev))
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
do iter=0,niter
  print *, 'calc G'  
  call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,mu,(/temps,tempd/))
  if (iter == 0) then     
    call calc_n_electron(G_lesser,G_greater,nen,En,NS,NB,nm_dev,nelec,pelec)    
  else    
    call calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,(/temps,tempd/),mu)    
    mu=(/ mus, mud /) - 0.2*sum(mu-(/ mus, mud /))/2.0d0 ! move Fermi level because Sig_GW shifts slightly the energies
    print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)    
  end if  
  call write_spectrum('ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0,-2.0/))
  call write_spectrum('ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum('pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0,-1.0/))
  !        
  print *, 'calc P'
  !call green_calc_polarization(nen,nen/2-10,En,nm_dev,G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,NB*NS*2)    
  nopmax=nen/2-10
  ndiag=NB*NS*2
  ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E+hw)
  ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E+hw) + Gij^r(E) * Gji^<(E+hw)
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
            P_retarded(i,l:h,iop) = P_retarded(i,l:h,iop) + dE* (G_lesser(i,l:h,ie) * conjg(G_retarded(i,l:h,ie-nop)) + G_retarded(i,l:h,ie) * G_lesser(l:h,i,ie-nop))        
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
  !call green_calc_w(nen,nen/2-10,En,nm_dev,NB*NS*2,V,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater)  
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
!  call write_spectrum('WR',iter,W_retarded,nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('WL',iter,W_lesser,  nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('WG',iter,W_greater, nen,En-en(nen/2),length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc SigGW'
  !call green_calc_gw_selfenergy(nen,nen/2-10,En,nm_dev,G_retarded,G_lesser,G_greater,W_retarded,W_lesser,W_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new,NB*NS*2)
  ndiag=NB*NS*2
  nopmax=nen/2-10
  Sig_greater = dcmplx(0.0d0,0.0d0)
  Sig_lesser = dcmplx(0.0d0,0.0d0)
  Sig_retarded = dcmplx(0.0d0,0.0d0)
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
  ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
  !$omp parallel default(none) private(ndiag,l,h,nop,ie,i,j,iop) shared(nopmax,Sig_lesser,Sig_greater,Sig_retarded,W_lesser,W_greater,W_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,dE)
  !$omp do
  do ie=1,nen
    do nop= -nopmax,nopmax    
      if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then      
        iop=nop+nen/2
        do i = 1,nm_dev   
          l=max(i-ndiag,1)
          h=min(nm_dev,i+ndiag)       
          Sig_lesser(i,l:h,ie)=Sig_lesser(i,l:h,ie)+dE*G_lesser(i,l:h,ie-nop)*W_lesser(i,l:h,iop)
          Sig_greater(i,l:h,ie)=Sig_greater(i,l:h,ie)+dE*G_greater(i,l:h,ie-nop)*W_greater(i,l:h,iop)
          Sig_retarded(i,l:h,ie)=Sig_retarded(i,l:h,ie)+dE*G_lesser(i,l:h,ie-nop)*W_retarded(i,l:h,iop) 
          Sig_retarded(i,l:h,ie)=Sig_retarded(i,l:h,ie)+dE*G_retarded(i,l:h,ie-nop)*W_lesser(i,l:h,iop) 
          Sig_retarded(i,l:h,ie)=Sig_retarded(i,l:h,ie)+dE*G_retarded(i,l:h,ie-nop)*W_retarded(i,l:h,iop)          
        enddo      
      endif
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  Sig_retarded(:,:,:) = dcmplx( dble(Sig_retarded(:,:,:)), aimag(Sig_greater(:,:,:)-Sig_lesser(:,:,:))/2.0d0 )
  ! mixing with previous ones
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
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
enddo                
deallocate(siglead)
deallocate(B)
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
    ! approximate F-D to step-function
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


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)
implicit none
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  integer i,j,k,nm,nmax
  COMPLEX(8) :: z
  real(8) :: E,error
  REAL(8) :: TOL=1.0D-10  ! [eV]
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


subroutine green_subspace_invert(nm,A,add_nm,in_method)
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

end subroutine green_subspace_invert


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
