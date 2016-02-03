      module mod_energy
      !
      ! Author: Veselin Kolev <vesso.kolev@gmail.com>
      ! Version: 2016020301
      !
      use mod_t
      use mod_dist

      implicit none

      contains

      subroutine init_ii_inn_energyM(atoms,N1,xdistMsq,ydistMsq,&
                                     zdistMsq,chargeM,params,&
                                     sigmaM,epsM,ii_inn_energyM)
      !
      ! This subroutine computes the interaction energy between the
      ! ions inside a box - the inner potential of interaction. It is
      ! a sum of the respective Coulomb electrostatic and
      ! Lennard-Jones potentials.
      !
      !
      ! atoms - 1D array of atomic descriptions (one for each atom)
      !         if type atom_t. That type is defined in the module
      !         mod_t.f90.
      !
      ! N1    - the total number of atoms (it is size(atoms)).
      !
      ! xdistMsq - 1D distance matrix containing the difference
      !            (atoms(i)%x-atoms(i)%y)**2.
      !
      ! ydistMsq - 1D distance matrix containing the difference
      !            (atoms(i)%y-atoms(i)%y)**2.
      !
      ! zdistMsq - 1D distance matrix containing the difference
      !            (atoms(i)%z-atoms(i)%z)**2.
      !
      ! chargeM - the product of the atomic charges
      !           atoms(i)%charge*atoms(j)%charge, represented as
      !           1D array.
      !
      ! params - params%LJ_cutoff supplies the subroutine with the
      !          Lennard-Jones the cut-off that helps to get
      !          faster the LJ potential. All ions inside the box
      !          are subject of cut-off test. If two of them are
      !          separated by a distance lower than params%LJ_cutoff
      !          their LJ-potential is computed and taken into account.
      !              
      ! sigmaM - the sigma matrix for van der Waals interactions.
      !
      ! The matrices xdistM,ydistM,zdistM are defined and handled in
      ! the module mod_dist.f90. Charge product matrix chargeM is
      ! defined and handled in mod_charge.f90.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      real(kind=real32), dimension(:), intent(in) :: chargeM
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), intent(in) :: sigmaM
      real(kind=real32), dimension(:), intent(in) :: epsM
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       ii_inn_energyM
      integer(kind=int32) :: i,j,k
      real(kind=real32) :: d

      allocate(ii_inn_energyM(N1*(N1-1)/2))

      do i=1,N1-1
         do j=i+1,N1
            k=get_index(i,j,N1)
            d=get_ii_inn_distance(k,xdistMsq,ydistMsq,zdistMsq)
            ii_inn_energyM(k)=chargeM(k)/d
            if (d.le.params%LJ_cutoff) then
               ii_inn_energyM(k)=ii_inn_energyM(k)+&
                                 get_LJ_inter(atoms,N1,i,j,d,sigmaM,&
                                 epsM)
            end if
         end do
      end do

      end subroutine init_ii_inn_energyM


      subroutine update_ii_inn_energyM(atoms,N1,i,&
                                       xdistMsq,ydistMsq,zdistMsq,&
                                       chargeM,params,sigmaM,epsM,&
                                       ii_inn_energyM)
      !
      ! Updates the inner electrostatic energy after changing the
      ! coordinates of the i-th atom located inside the simulation box.
      ! The matrix ii_inn_energyM is defined and initialized by the
      ! subroutine init_ii_inn_energyM introduced above.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      integer(kind=int32), intent(in) :: i
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      real(kind=real32), dimension(:), intent(in) :: chargeM
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), intent(in) :: sigmaM
      real(kind=real32), dimension(:), intent(in) :: epsM
      real(kind=real32), dimension(:), intent(inout) :: ii_inn_energyM
      integer(kind=int32) :: q,k
      real(kind=real32) :: d

      do q=1,N1
         if (i.ne.q) then
            if (i.gt.q) then
               k=get_index(q,i,N1)
            else
               k=get_index(i,q,N1)
            end if
            d=get_ii_inn_distance(k,xdistMsq,ydistMsq,zdistMsq)
            ii_inn_energyM(k)=chargeM(k)/d
            if (d.le.params%LJ_cutoff) then
               ii_inn_energyM(k)=ii_inn_energyM(k)+&
                                 get_LJ_inter(atoms,N1,i,q,d,sigmaM,&
                                 epsM)
            end if
         end if
      end do

      end subroutine update_ii_inn_energyM


      subroutine init_ii_out_energyM(atoms,N1,params,&
                                     xdistM,ydistM,zdistM,&
                                     xdistMsq,ydistMsq,zdistMsq,&
                                     chargeM,sigmaM,&
                                     epsM,ii_out_energyM)
      !
      ! Computes initially the interaction potential between
      ! the atoms inside the simulation box and those drawn in the
      ! image boxes - the outher interaction potential. It is a sum
      ! of the respective Coulomb electrostatic and Lennard-Jones
      ! potentials.
      !
      ! atoms - 1D array of atomic descriptions (one for each atom)
      !         if type atom_t. That type is defined in the module
      !         mod_t.f90.
      !
      ! N1 - the total number of the atoms located inside the
      !      simulation box. Should be computed by calling the function
      !      size outside the scope of this subroutine.
      !
      ! params - params%num_box_img_layers the number of box layers when
      !          create the images of the simulation box.
      !
      !          params%LJ_cutoff supplies the subroutine with the
      !          Lennard-Jones the cut-off that helps to get
      !          faster the LJ potential. All ions inside the box
      !          are subject of cut-off test. If two of them are
      !          separated by a distance lower than params%LJ_cutoff
      !          their LJ-potential is computed and taken into account.
      !
      !          See the type param_t in the module mod_t.f90 for the
      !          definitions of the properties and their types.
      !
      ! xdistM - The x-distance matrix between the atoms inside the
      !          simulation box. It is 1D array defined and initialized
      !          by the subroutines init_ldist_matrix from and
      !          update_ldist_matrix supplied by the module
      !          mod_dist.f90.
      !
      ! ydistM - The y-distance matrix between the atoms inside the
      !          simulation box. It is 1D array defined and initialized
      !          by the subroutines init_ldist_matrix from and
      !          update_ldist_matrix supplied by the module
      !          mod_dist.f90.
      !
      ! zdistM - The z-distance matrix between the atoms inside the
      !          simulation box. It is 1D array defined and initialized
      !          by the subroutines init_ldist_matrix from and
      !          update_ldist_matrix supplied by the module
      !          mod_dist.f90.
      !
      ! xdistMsq - the square of xdistM.
      !
      ! ydistMsq - the square of ydistM.
      !
      ! zdistMsq - the square of zdistM.
      !
      ! chargeM - The charge matrix computed for the atoms inside the
      !           simulaton box by the subroutine init_chargeM in the
      !           module mod_charge.f90.
      !
      ! sigmaM - the sigma matrix defined and init by the subroutine
      !          init_sigmaM introduced in this module. It is used to
      !          compute LJ-interactions.
      !
      ! epsM - the epsilon matrix defined and init by the subroutine
      !        init_epsM introduced in this module. It is used to
      !        compute LJ-interactions.
      !
      ! ii_out_energyM - The 1D matrix containing the outher
      !                  electrostatic potential of each atom located
      !                  inside the box.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), intent(in) :: xdistM
      real(kind=real32), dimension(:), intent(in) :: ydistM
      real(kind=real32), dimension(:), intent(in) :: zdistM
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      real(kind=real32), dimension(:), intent(in) :: chargeM
      real(kind=real32), dimension(:), intent(in) :: sigmaM
      real(kind=real32), dimension(:), intent(in) :: epsM
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       ii_out_energyM
      integer(kind=int32) :: i,j,k,l,m
      real(kind=real32) :: d,charge,mult,dummy,lj_dummy

      allocate(ii_out_energyM(N1))

      do i=1,N1
         ii_out_energyM(i)=0.0_real32
         do j=1,N1
            if (i.eq.j) then
               charge=atoms(i)%charge**2
               dummy=0.0_real32
               lj_dummy=0.0_real32
               do l=-params%num_box_img_layers,&
                     params%num_box_img_layers
                  do m=-params%num_box_img_layers,&
                        params%num_box_img_layers
                     if ((l.ne.0).and.(m.ne.0)) then
                        d=sqrt((l*params%box(1))**2+&
                               (m*params%box(2))**2)
                        dummy=dummy+1.0_real32/d
                        if (d.le.params%LJ_cutoff) then
                           lj_dummy=lj_dummy+&
                                   get_LJ_inter(atoms,N1,i,j,d,sigmaM,&
                                                epsM)
                        end if
                     end if
                  end do
               end do
            else
               dummy=0.0_real32
               lj_dummy=0.0_real32
               if (i.gt.j) then
                  k=get_index(j,i,N1)
                  charge=chargeM(k)
               else
                  k=get_index(i,j,N1)
                  charge=chargeM(k)
               end if
               do l=-params%num_box_img_layers,&
                     params%num_box_img_layers
                  do m=-params%num_box_img_layers,&
                        params%num_box_img_layers
                     if ((l.ne.0).and.(m.ne.0)) then
                        d=get_ii_out_distance(i,j,N1,params,xdistM,&
                                              ydistM,zdistM,xdistMsq,&
                                              ydistMsq,zdistMsq,l,m)
                        dummy=dummy+1.0_real32/d
                        if (d.le.params%LJ_cutoff) then
                           lj_dummy=lj_dummy+&
                                   get_LJ_inter(atoms,N1,i,j,d,sigmaM,&
                                                epsM)
                        end if
                     end if
                  end do
               end do
            end if
            ii_out_energyM(i)=ii_out_energyM(i)+dummy*charge+lj_dummy
         end do
      end do

      end subroutine init_ii_out_energyM


      subroutine update_ii_out_energyM(atoms,N1,i,params,&
                                       xdistM,ydistM,zdistM,&
                                       xdistMsq,ydistMsq,zdistMsq,&
                                       chargeM,sigmaM,epsM,&
                                       ii_out_energyM)
      !
      ! Updates the outher electrostatic energy after changing the
      ! coordinates of the i-th atom located inside the simulation box.
      ! The matrix ii_out_energyM is defined and initialized by the
      ! subroutine init_ii_out_energyM introduced above.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      integer(kind=int32), intent(in) :: i
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:), intent(in) :: xdistM
      real(kind=real32), dimension(:), intent(in) :: ydistM
      real(kind=real32), dimension(:), intent(in) :: zdistM
      real(kind=real32), dimension(:), intent(in) :: xdistMsq
      real(kind=real32), dimension(:), intent(in) :: ydistMsq
      real(kind=real32), dimension(:), intent(in) :: zdistMsq
      real(kind=real32), dimension(:), intent(in) :: chargeM
      real(kind=real32), dimension(:), intent(in) :: sigmaM
      real(kind=real32), dimension(:), intent(in) :: epsM
      real(kind=real32), dimension(:), intent(inout) :: ii_out_energyM
      integer(kind=int32) :: j,k,l,m
      real(kind=real32) :: d,charge,mult,dummy,lj_dummy


      ii_out_energyM(i)=0.0_real32
      do j=1,N1
         if (i.eq.j) then
            charge=atoms(i)%charge**2
            dummy=0.0_real32
            lj_dummy=0.0_real32
            do l=-params%num_box_img_layers,&
                  params%num_box_img_layers
               do m=-params%num_box_img_layers,&
                     params%num_box_img_layers
                  if (.not.((l.eq.0).and.(m.eq.0))) then
                     d=sqrt((l*params%box(1))**2+&
                            (m*params%box(2))**2)
                     dummy=dummy+1.0_real32/d
                     if (d.le.params%LJ_cutoff) then
                        lj_dummy=lj_dummy+&
                                 get_LJ_inter(atoms,N1,i,j,d,sigmaM,&
                                              epsM)
                     end if
                  end if
               end do
            end do
         else
            dummy=0.0_real32
            lj_dummy=0.0_real32
            if (i.gt.j) then
               k=get_index(j,i,N1)
               charge=chargeM(k)
            else
               k=get_index(i,j,N1)
               charge=chargeM(k)
            end if
            do l=-params%num_box_img_layers,&
                  params%num_box_img_layers
               do m=-params%num_box_img_layers,&
                     params%num_box_img_layers
                  if (.not.((l.eq.0).and.(m.eq.0))) then
                     d=get_ii_out_distance(i,j,N1,params,xdistM,&
                                           ydistM,zdistM,xdistMsq,&
                                           ydistMsq,zdistMsq,l,m)
                     dummy=dummy+1.0_real32/d
                     if (d.le.params%LJ_cutoff) then
                        lj_dummy=lj_dummy+&
                                 get_LJ_inter(atoms,N1,i,j,d,sigmaM,&
                                              epsM)
                     end if
                  end if
               end do
            end do
         end if
         ii_out_energyM(i)=ii_out_energyM(i)+dummy*charge+lj_dummy
      end do

      end subroutine update_ii_out_energyM


      subroutine init_ie_energyM(N1,N2,params,ie_xdistM,&
                                 ie_ydistM,ie_zdistM,chargeM_el,&
                                 ie_energyM)
      !
      ! Computes initially the Coulomb electrostatic potential between
      ! the atoms inside the simulation box and the dummy atoms placed
      ! on the electrodes.
      !
      ! N1 - the total number of the atoms located inside the
      !      simulation box. Should be computed by calling the function
      !      size outside the scope of this subroutine.
      !
      ! N2 - the total number of dummy atoms placed on both electrodes
      !      (sum of the dymmy atoms on left and right electrodes).
      !      Should be computed by calling the function size outside
      !      the scope of this subroutine.
      !
      ! params - params%num_box_img_layers the number of box layers when
      !          create the images of the simulation box.
      !
      !          See the type param_t in the module mod_t.f90 for the
      !          definitions of the properties and their types.
      !
      ! ie_xdistM - The x-distance matrix between the atoms inside the
      !             simulation box and the dummy atoms placed on the
      !             electrodes. This array is defined and initialized
      !             by the subroutine init_ie_ldist_matrix from the
      !             module mod_dist.f90. After each Monte Carlo
      !             proposal the ie_xdistM elements are updated by
      !             the subroutine update_ie_ldist_matrix in the module
      !             mod_dist.f90. ie_xdistM id 2D array.
      !
      ! ie_ydistM - The y-distance matrix between the atoms inside the
      !             simulation box and the dummy atoms placed on the
      !             electrodes. This array is defined and initialized
      !             by the subroutine init_ie_ldist_matrix from the
      !             module mod_dist.f90. After each Monte Carlo
      !             proposal the ie_xdistM elements are updated by
      !             the subroutine update_ie_ldist_matrix in the module
      !             mod_dist.f90. ie_xdistM id 2D array.
      !
      ! ie_zdistM - The z-distance matrix between the atoms inside the
      !             simulation box and the dummy atoms placed on the
      !             electrodes. This array is defined and initialized
      !             by the subroutine init_ie_ldist_matrix from the
      !             module mod_dist.f90. After each Monte Carlo
      !             proposal the ie_xdistM elements are updated by
      !             the subroutine update_ie_ldist_matrix in the module
      !             mod_dist.f90. ie_xdistM id 2D array.
      !
      ! chargeM_el - The charge matrix. Its elements are the product
      !              charge(i)*charge(j) where i is the index of the
      !              corresponding atom inside the simulation box and
      !              j is the undex of the respective dummy atom on
      !              the electrode. The matrix chargeM_el is defined
      !              and initialized by the subroutine
      !              init_charge_matrix_el in the module
      !              mod_energy.f90.
      !
      ! nbl - the number of box layers when create the images of the
      !       simulation box.
      !
      ! ie_energyM - The 2D matrix containg the Coulomb electrostatic
      !              potential between i-th atom located inside the
      !              simulation box and the m-th dummy atom placed on
      !              some of the electrodes.
      !
      ! Note that no LJ-potential is introduced in this particular case
      ! because the atoms in the simulation box and the dummy atoms on
      ! the electrodes are separated by gaps. It is enough to prevent
      ! the charged atoms to go critically close to the electrodes.
      !
      integer(kind=int32), intent(in) :: N1
      integer(kind=int32), intent(in) :: N2
      type(param_t), intent(in) :: params
      real(kind=real32), dimension(:,:), intent(in) :: ie_xdistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_ydistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_zdistM
      real(kind=real32), dimension(:,:), intent(in) :: chargeM_el
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       ie_energyM
      integer(kind=int32) :: i,j,l,m
      real(kind=real32) :: d

      allocate(ie_energyM(N1))

      do i=1,N1
         ie_energyM(i)=0.0_real32
         do j=1,N2
            do l=-params%num_box_img_layers,&
                  params%num_box_img_layers
               do m=-params%num_box_img_layers,&
                     params%num_box_img_layers
                  d=get_ie_distance(i,j,params,ie_xdistM,ie_ydistM,&
                                    ie_zdistM,l,m)
                  ie_energyM(i)=ie_energyM(i)+chargeM_el(i,j)/d
               end do
            end do
         end do
      end do

      end subroutine init_ie_energyM


      subroutine update_ie_energyM(N2,params,i,ie_xdistM,&
                                   ie_ydistM,ie_zdistM,chargeM_el,&
                                   ie_energyM)
      !
      ! Updates the ie_energyM matrix initialized by executing the
      ! subroutine init_ie_energyM (defined above) if the coordinates
      ! of the atom i-th atom inside the sumulation box are changed.
      !
      integer(kind=int32), intent(in) :: N2
      type(param_t), intent(in) :: params
      integer(kind=int32), intent(in) :: i
      real(kind=real32), dimension(:,:), intent(in) :: ie_xdistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_ydistM
      real(kind=real32), dimension(:,:), intent(in) :: ie_zdistM
      real(kind=real32), dimension(:,:), intent(in) :: chargeM_el
      real(kind=real32), dimension(:), intent(inout) ::ie_energyM
      integer(kind=int32) :: j,l,m
      real(kind=real32) :: d

      ie_energyM(i)=0.0_real32
      do j=1,N2
         do l=-params%num_box_img_layers,&
               params%num_box_img_layers
            do m=-params%num_box_img_layers,&
                  params%num_box_img_layers
               d=get_ie_distance(i,j,params,ie_xdistM,ie_ydistM,&
                                 ie_zdistM,l,m)
               ie_energyM(i)=ie_energyM(i)+chargeM_el(i,j)/d
            end do
         end do
      end do

      end subroutine update_ie_energyM


      subroutine init_sigmaM(atoms,N1,sigmaM)
      !
      ! The epsilon matrix which is used for computing Lennard-Jones
      ! interactions has the elements:
      !
      ! sigmaM(i,j)=0.5(sigma(i)+sigma(j))
      !
      ! Then the interactions can be computed according to the formula:
      !
      ! U(i,j)=4*epsM(i,j)*((sigmaM(i,j)/r)**12+(sigmaM(i,j)/r)**6)
      !
      ! To reduce the memory usage the symmetric matrix sigmaM is
      ! represented only with its upper triangular form, stored in
      ! 1D array with dimension N1*(N1-1)/2. The connection between
      ! the index of 1D array (k) and the indexes i and j (i>j) is
      ! provided by the formula:
      !
      ! k=N1(i-1)-i(i+1)/2+j
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       sigmaM
      integer(kind=int32) :: i,j

      allocate(sigmaM(N1*(N1-1)/2))

      do i=1,N1-1
         do j=i+1,N1
            sigmaM(get_index(i,j,N1))=(atoms(i)%sigma+&
                                       atoms(j)%sigma)/2
         end do
      end do
      
      end subroutine init_sigmaM


      subroutine init_epsM(atoms,N1,f,epsM)
      !
      ! The epsilon matrix which is used for computing van der Waals
      ! interactions has the elements:
      !
      ! epsM(i,j)=4*sqrt(eps(i)*eps(j))/f
      !
      ! Here f is a conversion constant that helps to compute the
      ! energy in 1/Angstrom units. Later, then the LJ interation
      ! potential
      !
      ! U(i,j)=epsM(i,j)*((sigmaM(i,j)/r)**12+(sigmaM(i,j)/r)**6)
      !
      ! is adding to the electrostatic one, the result have
      ! to be multiplied by f to get the energy in the proper units.
      ! The value of f is defined in the module mod_t.f90. The
      ! elements sigmaM(i,j) are supplied by the subroutine
      ! init_sigmaM above.
      !
      ! To reduce the memory usage the symmetric matrix epsM is
      ! represented only with its upper triangular form, stored in
      ! 1D array with dimension N1*(N1-1)/2. The connection between
      ! the index of 1D array (k) and the indexes i and j (i>j) is
      ! provided by the formula:
      !
      ! k=N1(i-1)-i(i+1)/2+j
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1
      real(kind=real32), intent(in) :: f
      real(kind=real32), dimension(:), allocatable, intent(out) :: &
       epsM
      integer(kind=int32) :: i,j

      allocate(epsM(N1*(N1-1)/2))

      do i=1,N1-1
         do j=i+1,N1
            epsM(get_index(i,j,N1))=sqrt(atoms(i)%eps*atoms(j)%eps)
         end do
      end do
      !
      ! Note that the constant f is defined in module mod_t.f90
      !
      epsM(:)=epsM(:)*4.0_real32/f

      end subroutine init_epsM


      function get_LJ_inter(atoms,N1,i,j,r,sigmaM,epsM) result(res)
      !
      ! The function computes the Lennard-Jones non-bonded
      ! interactions potential appearing between two atoms
      ! (identified by their indexes - i and j), separated from each
      ! other by distance r.
      !
      ! sigmaM - the sigma matrix. See the subroutine init_sigmaM
      !          above for details on the derivation of the elements.
      !
      ! epsM - the epsilon matrix. See the subroutine init_epsM
      !        above for details on the derivation of the elements.
      !
      type(atom_t), dimension(:), intent(in) :: atoms
      integer(kind=int32), intent(in) :: N1,i,j
      real(kind=real32), intent(in) :: r
      real(kind=real32), dimension(:), intent(in) :: sigmaM
      real(kind=real32), dimension(:), intent(in) :: epsM
      real(kind=real32) :: res,dummy

      if (i.eq.j) then
         dummy=(atoms(i)%sigma/r)**6
         res=4.0_real32*atoms(i)%eps*(dummy**2-dummy)
      else
         if (i.gt.j) then
            dummy=(sigmaM(get_index(j,i,N1))/r)**6
            res=epsM(get_index(j,i,N1))*(dummy**2-dummy)
         else
            dummy=(sigmaM(get_index(i,j,N1))/r)**6
            res=epsM(get_index(i,j,N1))*(dummy**2-dummy)
         end if
      end if

      end function get_LJ_inter


      end module mod_energy
