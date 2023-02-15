!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf

implicit none 

private

public :: green_rgf_cms

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)

CONTAINS

    ! RGF for diagonal blocks of W^r,<,>
    subroutine green_rgf_w(nm,nx,Vii,p_lesser,p_greater,p_retarded,w_lesser,w_greater,w_retarded)
        implicit none
        integer,intent(in)::nm,nx
        complex(8),intent(in),dimension(nm,nm,nx) :: Vii,P_lesser,P_greater,P_retarded
        complex(8),intent(inout),dimension(nm,nm,nx) :: W_lesser,W_greater,W_retarded
        ! --------
    
    end subroutine green_rgf_w


    ! RGF for diagonal blocks of G^r,<,>
    subroutine green_RGF_CMS(TEMP,nm,nx,E,mu,Sii,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur) 
        implicit none    
        integer,intent(in)::nm,nx
        real(8),intent(in) :: E,mu(2),temp(2)
        COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph
        COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Sii,Hii,H1i
        COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: cur,ldos,ndens,pdens
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
            pdens(:,:,l)=Gp(:,:)!Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))!
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


end module green_rgf
