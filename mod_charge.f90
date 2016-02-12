      module mod_charge
      !
      ! Author: Vesselin Kolev <vesso.kolev@gmail.com>
      ! Version: 2016021001
      !
      ! The module mod_charge handles the atomic charges participating
      ! in the computation of the electrostatic potential.
      !
      use mod_t

      implicit none

      contains

      subroutine init_charge_matrix(atoms,params,chargeM)
      !
      ! Computes the matrix of the atomic charges product, each element
      ! of which is computed (for a certain atomic indexes i and j) as:
      !
      ! qi*qj
      !
      ! Because the matrix qi*qj is symmetric it can be represented as
      ! an 1D array with N1*(N1+1)/2 element (along with the elements
      ! of the diagonal). To switch from indexes i and j to the index
      ! of the 1D array, k, the function get_index_d introduced in
      ! mod_t.f90 is implemented.
      !
      ! NOTE: Don't look for a subroutine that updates the charges
      ! matrix during the simulation because the atomic charges remain
      ! constant.
      !
      !
      ! Arguments:
      !
      ! atoms - 1D array of atom_t type elements representing the atoms
      !         defined in the system. The properties of the atoms are
      !         described and the properties of the atom_t type in the
      !         module mod_t.f90.
      !
      ! params - the parameters of the simulations. More details can be
      !          found in the description of the type param_t in the
      !          module mod_t.f90.
      !
      ! chargeM - 1D array of size N1*(N1+1)/2 used to store both
      !           diagonal elements and those above the diagonal of the
      !           atomic charge product matrix.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       chargeM
      integer(kind=int32) :: i,j

      allocate(chargeM(params%N1*(params%N1+1)/2))

      do i=1,params%N1
         do j=i,params%N1
            chargeM(get_index_d(i,j,params%N1))=atoms(i)%charge*&
                                                atoms(j)%charge
         end do
      end do

      end subroutine init_charge_matrix


      subroutine init_charge_matrix_el(atoms,grid,params,chargeM_el)
      !
      ! Computes the analogue of chargeM matrix (see the subroutune
      ! init_charge_matrix introduced above) in the case of two groups
      ! of atoms charges - the first one contains the atomic charges
      ! of the atoms inside the simulation box, and the second group
      ! represents the atomic charges of the grid atoms (placed on the
      ! planes of the electrodes). In this particular case the matrix
      ! of the atomic charges qi*qj (here i is the index of the
      ! elements of the first group, and j - of the second one) is not
      ! symmetric. So it is not possible to apply the simplification
      ! used for chargeM matrix (above).
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
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:,:), allocatable, intent(out) :: &
       chargeM_el
      integer(kind=int32) :: i,j

      allocate(chargeM_el(params%N1,2))

      do i=1,params%N1
         do j=1,params%N2
            chargeM_el(i,grid(j)%plane_id)=atoms(i)%charge*&
            params%dummy_atom_charge(grid(j)%plane_id)
         end do
      end do

      end subroutine init_charge_matrix_el


      end module mod_charge
