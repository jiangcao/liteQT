!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 01/2023
! Phonon dynamics
MODULE phononDyn

IMPLICIT NONE

private
real(8),allocatable::C(:,:,:,:,:,:,:)  ! 2nd order IFC, (dir,dir,atom1,atom2,cell_x,cell_y,cell_z) [eV/ang^2]
real(8),allocatable::W(:,:,:,:,:,:,:)  ! mass scaled dynamic matrix (dir,dir,atom1,atom2,cell_x,cell_y,cell_z) [1/sec^2]
real(8),allocatable::atom_pos(:,:) ! atom positions (unit cell) [ang]
real(8),allocatable::atom_mass(:) ! atom mass (a.u.)
real(8),allocatable::Ze(:,:,:)   ! Effective charge tensor
REAL(8), DIMENSION(3) :: alpha,beta,gamm,xhat,yhat,b1,b2
REAL(8) :: Lx, Ly,cell(3,3) ! in Ang
INTEGER :: zmin,zmax,ymin,ymax,xmin,xmax,nx,ny,nz,NA,NB
! neighbor index range 
! nb is the number of phonon bands = number of atoms NA*3
COMPLEX(8), PARAMETER :: zzero = dcmplx(0.0d0,0.0d0)
COMPLEX(8), PARAMETER :: z1j = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0


CONTAINS

SUBROUTINE phon_free_memory()
if (allocated(C)) deallocate(C)
if (allocated(atom_pos)) deallocate(atom_pos)
END SUBROUTINE phon_free_memory

! load from a force-constant file
SUBROUTINE phon_load_from_file()

END SUBROUTINE phon_load_from_file

END MODULE phononDyn

