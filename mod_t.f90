      module mod_t
      !
      ! Created by Vesselin Kolev <vesso.kolev@gmail.com>
      ! 2016020101
      !
      ! Defines the kind of the real and integer parameters, and the
      ! types used to pack the atom properties into the memory. Some
      ! important constant are also described here.
      !
      ! Note: if one uses gfortran of Intel Fortran compilers it is
      ! possible to get these numbers as they are efinied by the module
      ! iso_fortran_env. It might be called as:
      !
      ! use iso_fortran_env, only: real32, real64, int8, int16, int32,&
      !                            int64
      !
      ! But PGI Fortran compiler does not provided iso_fortran_env.
      ! That's why for the sake of compatibility of the code these
      ! parameters are defined in the code bellow.
      !
      implicit none

      integer, parameter :: real32=4
      integer, parameter :: real64=8

      integer, parameter :: int8=1
      integer, parameter :: int16=2
      integer, parameter :: int32=4
      integer, parameter :: int64=8

      type atom_t
         integer(kind=int64) :: anum
         character(3) :: resname
         character(4) :: aname
         real(kind=real32) :: x
         real(kind=real32) :: y
         real(kind=real32) :: z
         real(kind=real32) :: sigma
         real(kind=real32) :: eps
         integer(kind=int32) :: spatial_grp
         real(kind=real32) :: charge
      end type atom_t

      type el_atom_t
         real(kind=real32) :: x
         real(kind=real32) :: y
         real(kind=real32) :: z
         real(kind=real32) :: charge
      end type el_atom_t

      type param_t
         real(kind=real32), dimension(3) :: box
         real(kind=real32) :: left_el_charge
         real(kind=real32) :: right_el_charge
         real(kind=real32) :: el_grid_dens_x
         real(kind=real32) :: el_grid_dens_y
         real(kind=real32) :: core_dist_el
         real(kind=real32) :: core_dist_mem
         integer(kind=int64) :: accept_proposals
      end type param_t

      real(kind=real32), parameter :: &
      pi=3.141592653589793116_real32
      real(kind=real32), parameter :: &
      pix2=6.283185307179586232_real32
      real(kind=real32), parameter :: &
      pix3div2=4.712388980384689674_real32

      !
      ! The following parameters should be defined:
      !
      ! - in the module mod_input.f90 (see subroutine readInput there):
      ! atoms_s=size(atoms,1)
      ! bonds_s=size(bonds,1)
      ! pdihs_s=size(pdihs,1)
      !
      integer(kind=int64) :: atoms_s,bonds_s,pdihs_s


      contains

      function get_index(i,j,N) result(k)
      !
      ! Computes the 1D-index of the 1D-distance matrix. Here i and j
      ! are the i,j indexes of the atoms as they should appear in the
      ! usual 2D distance matrix that has NxN elements, where N is the
      ! total number of atoms.
      !
      ! The following formula is used to get the 1D index:
      !
      ! k=N(i-1)-i(i+1)/2+j
      !
      ! Note that thus derived 1D representation of distance matrix has
      ! N*(N-1)/2 elements.
      !
      integer(kind=int32), intent(in) :: i
      integer(kind=int32), intent(in) :: j
      integer(kind=int32), intent(in) :: N
      integer(kind=int32) :: k

      k=N*(i-1)-i*(i+1)/2+j

      end function get_index

      end module mod_t
