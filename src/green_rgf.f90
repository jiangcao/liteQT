!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf

implicit none 

private

public :: green_rgf_cms

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
  ! ------- 
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
!!  Glp=0.0D0
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

end module green_
