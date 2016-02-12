      module mod_t
      !
      ! Created by Vesselin Kolev <vesso.kolev@gmail.com>
      ! 2016021001
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
         integer(kind=int32) :: kind
      end type atom_t

      type el_atom_t
         real(kind=real32) :: x
         real(kind=real32) :: y
         real(kind=real32) :: z
         integer(kind=int8) :: plane_id
      end type el_atom_t

      type param_t
         real(kind=real32), dimension(3) :: box
         real(kind=real32), dimension(3) :: half_box
         real(kind=real32) :: left_el_charge
         real(kind=real32) :: right_el_charge
         real(kind=real32) :: el_grid_dens_x
         real(kind=real32) :: el_grid_dens_y
         real(kind=real32) :: core_dist_el
         real(kind=real32) :: core_dist_mem
         integer(kind=int64) :: accept_proposals
         integer(kind=int32) :: num_box_img_layers
         real(kind=real32) :: LJ_cutoff
         real(kind=real32) :: LJ_scale_factor
         real(kind=real32) :: left_el_gap_pos
         real(kind=real32) :: right_el_gap_pos
         real(kind=real32) :: left_mem_gap_pos
         real(kind=real32) :: right_mem_gap_pos
         real(kind=real32) :: crit_contact_dist
         real(kind=real32), dimension(3) :: max_displacement
         integer(kind=int32) :: N1
         integer(kind=int32) :: N2
         real(kind=real32) :: f
         real(kind=real32), dimension(2) :: dummy_atom_charge
      end type param_t

      real(kind=real32), parameter :: &
      pi=3.141592653589793116_real32

      real(kind=real32), parameter :: &
      pix2=6.283185307179586232_real32

      real(kind=real32), parameter :: &
      pix3div2=4.712388980384689674_real32

      real(kind=real32), parameter :: &
      f=17.369409810401915_real32 ! kJ/mol*A

      real(kind=real32), parameter :: &
      kT=2.479_real32 ! kJ/mol at 298K

      real(kind=real32), parameter :: &
      f_red=f/kT

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
      ! Computes the 1D-index of the 1D representation of the 2D
      ! symmetric matrix at i<j. Therefore only the elements above the
      ! diagonal of the matrix are taken into account here which is
      ! enough to get all unique elements of the 2D symmetrix matrix.
      !
      ! The following formula is used to get the 1D index based on a
      ! given pair i,j (i<j):
      !
      ! k=N(i-1)-i(i+1)/2+j
      !
      ! Note that the 1D array in this case has N*(N-1)/2 elements
      ! in total. 
      !
      integer(kind=int32), intent(in) :: i
      integer(kind=int32), intent(in) :: j
      integer(kind=int32), intent(in) :: N
      integer(kind=int32) :: k

      k=N*(i-1)-i*(i+1)/2+j

      end function get_index


      function get_index_d(i,j,N) result(k)
      !
      ! This is modification of the function get_index introduced
      ! before. It takes into account both 2D matrix elements above
      ! the diagonal as well as the diagonal elements. Therefore in
      ! this particular case the 1D matrix represents the elements
      ! of 2D matrix when i<=j. The index of the 1D array can be
      ! calculated for any i and j (i<=j) through the formula:
      !
      ! k=N(i-1)-i(i-1)/2+j
      !
      ! Note that the size of thus generated 1D matrix is N*(N+1)/2.
      !
      integer(kind=int32), intent(in) :: i
      integer(kind=int32), intent(in) :: j
      integer(kind=int32), intent(in) :: N
      integer(kind=int32) :: k

      k=N*(i-1)-i*(i-1)/2+j

      end function get_index_d


      function get_rand_int(ll,ul) result(res)
      !
      ! Computes random integer i, so
      !
      ! (1) great or equal to a
      !
      ! (2) less or equal to b
      !
      integer, intent(in) :: ll, ul
      integer :: res
      real :: rand

      call random_number(rand)

      res=ll+int(rand*(ul-ll))+1

      end function get_rand_int


      subroutine rand(seed,rnd)
      !
      ! This is a function that generates an uniformly distributed
      ! floating point random number in [0,1). To call the function one
      ! need to supply as an argument a seed number (that number have
      ! to be of integer type). The result supplied by the function is
      ! based on a linear congruential generator.
      !
      !
      ! The constants l,c, and m bellow are the parameters of the linear
      ! congruential generator. Do not change these values unless you
      ! know exectly what are you doing.
      !
      integer(kind=int32) :: l,c,m
      !
      ! The seed is defined as an varianble that can be changed during
      ! the execution in the body of the function. Hence it is marked
      ! with intent(inout) declaration.
      !
      integer(kind=int32) :: seed
      ! This is the result of the function. It is a floating point
      ! number.
      real(kind=real32) :: rnd,res
      !
      ! This line re-computes value of seed by using the function mod.
      ! The function mod(a,b) computes the remainder of the division
      ! of a by b. For instance: mod(5,9)=5.
      !
      l=1029
      c=221591
      m=1048576
      seed = mod(seed*l+c,m)
      !
      ! Here we compute an uniformly distributed floating point random
      ! number in [0,1) and supply it as a value of res.
      !
      rnd=float(seed)/float(m)
      !
      end


      function sum_(oneDarray,dimen) result(res)
      real(kind=real32), dimension(:), intent(in) :: oneDarray
      integer(kind=real32), intent(in) :: dimen
      real(kind=real32) :: res
      integer(kind=real32) :: i

      res=0.0_real32

      do i=1,dimen
         res=res+oneDarray(i)
      end do

      end function sum_


      subroutine get_random_displacement(seed,amplitude,displacement)
      integer(kind=real32), intent(inout) :: seed
      real(kind=real32), dimension(3), intent(in) :: amplitude
      real(kind=real32), dimension(3), intent(out) :: displacement
      integer(kind=int8) :: i
      real(kind=real32) :: dummy,sign_

      do i=1,3
         call random_number(sign_)
         if (sign_.lt.0.5_real32) then
            sign_= 1.0_real32
         else
            sign_=-1.0_real32
         end if
         call random_number(dummy)
         displacement(i)=sign_*dummy*amplitude(i)
      end do

      end subroutine get_random_displacement


      subroutine extend1DFloatArray(array,newElement)
      real(kind=real32), dimension(:), allocatable, intent(inout) :: &
       array
      real(kind=real32), intent(in) :: newElement
      real(kind=real32), dimension(:), allocatable :: temp
      integer(kind=int64) :: i,array_s

      array_s=size(array,1)

      call move_alloc(array,temp)

      array_s=array_s+1

      allocate(array(array_s))

      array(1:array_s-1)=temp(:)

      deallocate(temp)

      array(array_s)=newElement

      end subroutine extend1DFloatArray


      end module mod_t
