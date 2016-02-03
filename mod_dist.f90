      module mod_dist
      !
      ! Author: Vesselin Kolev <vesso.kolev@gmail.com>
      ! Version: 2016020201
      !
      ! This module handles all distance matrices needed for computing
      ! the electrostatic of Lenard-Jones interactions.
      !
      use mod_t

      implicit none


      contains


      subroutine init_ldist_matrix(atoms,dimen,ldistM,ldistMsq)
      !
      ! Computes the matrix of distances between the atoms (the
      ! electrodes and the membrane atoms are excluded here) placed
      ! inside the simulation box with respect to a certain coordinate
      ! direction only.
      !
      ! atoms - array of atom descriptions. It is created and
      !         inititalized in the main by calling the subroutines in
      !         mod_input.f90. The type of the the elements of the 
      !         array atoms is atom_t, defined in mod_t.f90.
      !
      ! dimen - the dimension (x, y, or z) the distances should be
      !         computed with respect to. One of the following integers
      !         should be assigned to dimen:
      !
      !         1 - to take the distances in x-direction
      !         2 - to take the distances in y-direction
      !         3 - to take the distances in z-direction
      !
      ! ldistM - 1D array containing the elements of the upper triangle
      !          of the respective distance metrix.
      !
      ! ldistMsq - 1D array containing ldistM(k)**2
      !
      ! Because the distance matrix here is symmetric it can be
      ! represented as the 1D array ldistM which has with N*(N-1)/2
      ! elements. To switch from indexes i and j of the ordinary two-
      ! dimension matrix to the index of the 1D array, denoned as k,
      ! bellow, the following conversion formula should be implemented:
      !
      ! k=N(i-1)-i(i+1)/2+j
      !
      ! The implementation of the conversion formula is the function
      ! get_index defined in mod_t.f90
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int8), intent(in) :: dimen
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       ldistM
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       ldistMsq
      integer(kind=int32) :: N
      integer(kind=int32) :: i,j

      N=size(atoms)

      allocate(ldistM(N*(N-1)/2))
      allocate(ldistMsq(N*(N-1)/2))
      !
      ! Note the special way the indexes i and j are changing in the
      ! loop. Only the elements of the upper triangle of the distance
      ! matrix are generating. The diagonal elements are not needed
      ! because they are all zeros.
      ! 
      do i=1,N-1
         do j=i+1,N
            if (dimen.eq.1) then
               ldistM(get_index(i,j,N))=atoms(i)%x-atoms(j)%x
            end if
            if (dimen.eq.2) then
               ldistM(get_index(i,j,N))=atoms(i)%y-atoms(j)%y
            end if
            if (dimen.eq.3) then
               ldistM(get_index(i,j,N))=atoms(i)%z-atoms(j)%z
            end if
         end do
      end do
      !
      ! Get the squares of the distances.
      !
      ldistMsq=ldistM**2

      end subroutine init_ldist_matrix


      subroutine update_ldist_matrix(atoms,dimen,i,N,ldistM,ldistMsq)
      !
      ! When atom i is moved this subroutine updates the elements of
      ! ldistM matrix related to that atom. The goal is to create
      ! highly efficient update routine that speeds up the execution
      ! of the simulation.
      !
      ! Because the function get_index from within the module mod_t.f90
      ! have to be called many times bellow it is good idea to save
      ! computation resources by supplying from outside the total
      ! number of atoms (N), instead of computing it by calling the
      ! intrinsic function size here.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int8), intent(in) :: dimen
      integer(kind=int32), intent(in) :: i
      integer(kind=int32), intent(in) :: N
      real(kind=real32), dimension(:), intent(inout) :: ldistM
      real(kind=real32), dimension(:), intent(inout) :: ldistMsq
      integer(kind=int32) :: q,k
      real(kind=real32) :: coord1,coord2

      do q=1,N
         if (i.ne.q) then
            if (dimen.eq.1) then
               coord1=atoms(q)%x
               coord2=atoms(i)%x
            end if
            if (dimen.eq.2) then
               coord1=atoms(q)%y
               coord2=atoms(i)%y
            end if
            if (dimen.eq.3) then
               coord1=atoms(q)%z
               coord2=atoms(i)%z
            end if
            !
            ! Why we do this if else loop bellow? The function
            ! get_index works only if the first argument is less than
            ! the second one. So if there is a case when it is the
            ! othere way round, the first and second argument in the
            ! argument list of get_index should get their pisitions
            ! replaced.
            if (i.gt.q) then
               k=get_index(q,i,N)
               ldistM(k)=coord1-coord2
               ldistMsq(k)=ldistM(k)**2
            else
               k=get_index(q,i,N)
               ldistM(get_index(i,q,N))=coord2-coord1
               ldistMsq(k)=ldistM(k)**2
            end if
         end if
      end do

      end subroutine update_ldist_matrix


      subroutine init_ie_ldist_matrix(atoms,N1,grid,N2,dimen,ldistM)
      !
      ! Computes the matrix of distances between the ions placed
      ! inside the box and the dummy atoms spread on the electrodes
      ! as grid nodes by the subroutine generate_electrode_grid
      ! supplied by mod_electrode.f90. The result is 2D matrix and it
      ! cannot be reduced to 1D representation (as it is done above for
      ! the box atoms). Crearly it is because this is a distance matrix
      ! between two separate group of atoms. So no reduction in the
      ! dimension is possible.
      !
      ! atoms - array of atom descriptions. It is created and
      !         inititalized in the main by calling the subroutines in
      !         mod_input.f90. The type of the the elements of the 
      !         array atoms is atom_t, defined in mod_t.f90.
      !
      ! N1 - the total number of atoms inside the box. There is no
      !      need this number to be computed every time because N1 is
      !      a parameter of the system and it is constant number during
      !      the simulations.
      !
      ! grid - array of atom descriptions. It is created and
      !        inititalized in generate_electrode_grid subroutine in
      !        the module mod_electrode.f90. The type of the the
      !        elements of the grid array is el_atom_t, as defined in
      !        mod_t.f90.
      !
      ! N2 - the total number of atoms on the electrode. There is no
      !      need this number to be computed every time because N2 is
      !      a parameter of the system and it is constant number during
      !      the simulations.
      !
      ! dimen - the dimension (x, y, or z) the distances should be
      !         computed with respect to. One of the following integers
      !         should be assigned to dimen:
      !
      !         1 - to take the distances in x-direction
      !         2 - to take the distances in y-direction
      !         3 - to take the distances in z-direction
      !
      ! ldistM - 2D matrix contining the distance between atom i in the
      !          box and dummy atom j on an electrode, with respect to
      !          the dimension given by dimen.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      type(el_atom_t), dimension(:), intent(in) :: grid
      integer(kind=int32), intent(in) :: N2
      integer(kind=int8), intent(in) :: dimen
      real(kind=real32), dimension(:,:), allocatable, intent(out) :: &
       ldistM
      integer(kind=int32) :: i,j

      allocate(ldistM(N1,N2))

      do i=1,N1
         do j=1,N2
            if (dimen.eq.1) then
               ldistM(i,j)=atoms(i)%x-grid(j)%x
            end if
            if (dimen.eq.2) then
               ldistM(i,j)=atoms(i)%y-grid(j)%y
            end if
            if (dimen.eq.3) then
               ldistM(i,j)=atoms(i)%z-grid(j)%z
            end if
         end do
      end do

      end subroutine init_ie_ldist_matrix


      subroutine update_ie_ldist_matrix(atoms,grid,N2,i,dimen,ldistM)
      !
      ! This subroutine updates the distance matrix ldistM to take into
      ! account the new distances between i-th atom inside the box and
      ! the dummy atoms on the electrode.
      !
      ! Note than the matrix ldistM is 2D array. Also note that there
      ! is no need to know the number of elelemtns in the array atoms,
      ! because moving the i-th atom in the box changes only one row
      ! in ldistM. This row has N2 elements - the number of the dummy
      ! atoms on the electrode. So only N2 is needed here.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      type(el_atom_t), dimension(:), intent(in) :: grid
      integer(kind=int32), intent(in) :: N2
      integer(kind=int32), intent(in) :: i
      integer(kind=int8), intent(in) :: dimen
      real(kind=real32), dimension(:,:), intent(inout) :: &
       ldistM
      integer(kind=int32) :: j

      do j=1,N2
         if (dimen.eq.1) then
            ldistM(i,j)=atoms(i)%x-grid(j)%x
         end if
         if (dimen.eq.2) then
            ldistM(i,j)=atoms(i)%y-grid(j)%y
         end if
         if (dimen.eq.3) then
            ldistM(i,j)=atoms(i)%z-grid(j)%z
         end if
      end do

      end subroutine update_ie_ldist_matrix


      function get_ii_inn_distance(k,xdistMsq,ydistMsq,zdistMsq) &
                                   result(res)
      !
      ! Computes the distance between the two atoms inside the
      ! simulation boxcentral box - the inner distance. If the indexes
      ! of the atoms are i and j the index k should be computed
      ! by using the function get_index supplied by the module
      ! mod_t.f90.
      !
      ! k - index corresponding to i and j (the indexes of the atoms in
      !     the box through the function get_index (see mod_t.f90).
      !
      ! xdistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in x-direction.
      !
      ! ydistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in y-direction.
      !
      ! zdistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in z-direction.
      !
      integer(kind=int32), intent(in) :: k
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      real(kind=real32) :: res

      res=sqrt(xdistMsq(k)+ydistMsq(k)+zdistMsq(k))

      end function get_ii_inn_distance


      function get_ii_out_distance(i,j,N,params,&
                                   xdistM,ydistM,zdistM,&
                                   xdistMsq,ydistMsq,zdistMsq,&
                                   p,q) result(res)
      !
      ! Computes the distance between two atoms with respect to the
      ! box images - outer distance. The first atom (with index i) is
      ! located inside the simulation box and the second one (j) is
      ! inside any of the box images (the image is identified by the
      ! values of the shift parameters p and q).
      !
      ! k - index corresponding to i and j (the indexes of the atoms in
      !     the box through the function get_index (see mod_t.f90).
      !
      ! xdistM - the 1D distance matrix supplied by either
      !          init_ldist_matrix or update_ldist_matrix
      !          subroutines. It stores the distances between the atoms
      !          in x-direction.
      !
      ! ydistM - the 1D distance matrix supplied by either
      !          init_ldist_matrix or update_ldist_matrix
      !          subroutines. It stores the distances between the atoms
      !          in y-direction.
      !
      ! zdistM - the 1D distance matrix supplied by either
      !          init_ldist_matrix or update_ldist_matrix
      !          subroutines. It stores the distances between the atoms
      !          in z-direction.
      !
      ! xdistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in x-direction.
      !
      ! ydistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in y-direction.
      !
      ! zdistMsq - the 1D square distance matrix supplied by either
      !            init_ldist_matrix or update_ldist_matrix
      !            subroutines. It stores the distances between the atoms
      !            in z-direction.
      !
      ! Computes the distance between ions in box and the images of
      ! the box. Here i and j can be any. No limitations on 
      ! 
      integer(kind=int32), intent(in) :: i,j,N
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), intent(in) :: xdistM
      real(kind=real32), dimension(:), intent(in) :: ydistM
      real(kind=real32), dimension(:), intent(in) :: zdistM
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      integer(kind=int32), intent(in) :: p,q
      real(kind=real32) :: res
      integer(kind=int32) :: k

      if (i.eq.j) then
         res=sqrt((p*params%box(1))**2+(q*params%box(2))**2)
      else
         if (i.gt.j) then
            k=get_index(j,i,N)
            if ((p.ne.0).and.(q.ne.0)) then
               res=sqrt((-xdistM(k)-p*params%box(1))**2+&
                        (-ydistM(k)-q*params%box(2))**2+&
                          zdistMsq(k))
            else
               if (p.eq.0) then
                  res=sqrt(xdistMsq(k)+&
                           (-ydistM(k)-q*params%box(2))**2+&
                             zdistMsq(k))
               end if
               if (q.eq.0) then
                  res=sqrt((-xdistM(k)-p*params%box(1))**2+&
                           ydistMsq(k)+&
                           zdistMsq(k))
               end if
            end if
         else
            k=get_index(i,j,N)
            if ((p.ne.0).and.(q.ne.0)) then
               res=sqrt((xdistM(k)-p*params%box(1))**2+&
                        (ydistM(k)-q*params%box(2))**2+&
                         zdistMsq(k))
            else
               if (p.eq.0) then
                  res=sqrt( xdistMsq(k)+&
                           (ydistM(k)-q*params%box(2))**2+&
                            zdistMsq(k))
               end if
               if (q.eq.0) then
                  res=sqrt((xdistM(k)-p*params%box(1))**2+&
                            ydistMsq(k)+&
                            zdistMsq(k))
               end if
            end if
         end if
      end if

      end function get_ii_out_distance



      function get_ie_distance(i,m,params,ie_xdistM,ie_ydistM,&
                               ie_zdistM,p,q) result(res)
      !
      ! This function computes the distance between i-th atom located
      ! inside the simulation box and m-th dummy atom placed on the
      ! electrode plane.
      !
      ! ie_xdistM - 2D matrix containing the x-distances between the
      !             atoms inside the simulation box and the dummy
      !             atoms placed on the electrode plane.
      !
      ! ie_ydistM - 2D matrix containing the y-distances between the
      !             atoms inside the simulation box and the dummy
      !             atoms placed on the electrode plane.
      !
      ! ie_zdistM - 2D matrix containing the z-distances between the
      !             atoms inside the simulation box and the dummy
      !             atoms placed on the electrode plane.
      !
      ! p - the shift of the simulaiton box in x
      !
      ! q - the shift of the simulaiton box in y
      !
      integer(kind=int32), intent(in) :: i,m
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:,:), intent(in) :: ie_xdistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_ydistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_zdistM
      integer(kind=int32), intent(in) :: p,q
      real(kind=real32) :: res

      res=sqrt((ie_xdistM(i,m)-p*params%box(1))**2+&
               (ie_ydistM(i,m)-q*params%box(2))**2+&
                ie_zdistM(i,m)**2)

      end function get_ie_distance


      end module 

