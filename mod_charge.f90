      module mod_charge
      !
      ! Author: Vesselin Kolev <vesso.kolev@gmail.com>
      ! Version: 2016020102
      !
      ! The module mod_charge handles the atomic charges participating
      ! in the computation of the electrostatic potential.
      !
      use mod_t

      implicit none

      contains

      subroutine init_charge_matrix(atoms,chargeM)
      !
      ! Computes the matrix of the atomic charges, each element of
      ! which is the product qi*qj. The atomic charges are supplied
      ! as a type property of pdbatoms. See mod_t.f90 for the complete
      ! definitions of the type atom_t.
      !
      ! Because the matrix qi*qj is symmetric it can be represented as
      ! an 1D array with N*(N-1)/2 element. To switch from indexes
      ! i and j to the index of 1D array, k, the following conversion
      ! should be implemented:
      !
      ! k=N(i-1)-i(i+1)/2+j
      !
      ! The implementation of the conversion formula is the function
      ! get_index defined in mod_t.f90
      !
      ! Note: Don't look for a subroutine to update the charges matrix
      ! during the sumulations. The atomic charges are constants and
      ! they do not change here.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       chargeM
      integer(kind=int32) :: N
      integer(kind=int32) :: i,j

      N=size(atoms,1)

      allocate(chargeM(N*(N-1)/2))

      do i=1,N-1
         do j=i+1,N
            chargeM(get_index(i,j,N))=atoms(i)%charge*&
                                      atoms(j)%charge
         end do
      end do

      end subroutine init_charge_matrix


      subroutine init_charge_matrix_el(atoms,grid,chargeM_el)
      !
      ! Computes the matrix of the atomic charges, each element of
      ! which is the product qi*qj. Here index i denotes atoms inside
      ! the box while j is the index of the corresponding dummy atom on
      ! the electrode. Here it is not possible to reduce the 2D matrix
      ! to 1D as it is done in init_charge_matrix, because the charges
      ! are separated into two groups of atoms.
      !
      ! Note: Don't look for a subroutine to update the charges matrix
      ! during the sumulations. The atomic charges are constants and
      ! they do not change here.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      type(el_atom_t), dimension(:), allocatable :: grid
      real(kind=real32), dimension(:,:), allocatable, intent(out) :: &
       chargeM_el
      integer(kind=int32) :: N1,N2
      integer(kind=int32) :: i,j

      N1=size(atoms)
      N2=size(grid)

      allocate(chargeM_el(N1,N2))

      do i=1,N1
         do j=1,N2
            chargeM_el(i,j)=atoms(i)%charge*&
                            grid(j)%charge
         end do
      end do

      end subroutine init_charge_matrix_el


      end module mod_charge
