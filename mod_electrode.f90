      module mod_electrode
      !
      ! Author: Vesselin Kolev <vesso.kolev@gmail.com>
      ! Version: 2016020101
      !
      ! This module helps to draw the grid dummy atoms on each of the
      ! electrodes that are defined on both ends of the box in z.
      ! 
      use mod_t

      implicit none

      contains

      subroutine generate_electrode_grid(params,grid)
      !
      ! Generates the grid of electrode dummy atoms. Each such atom
      ! is represented by its coordinates and charge. The
      ! corresponding type of dummy atoms is el_atom_t and it can be
      ! found defined in mod_t.f90.
      !
      ! params - the parameters of the system. The type param_t is
      !          defined in mod_t.f90. Here params is needed because
      !          it contains the dimensions of the simulation box.
      !
      ! grid - an 1D array of records of type el_atom_t (defined also
      !        in mod_t.f90) that describe the generated grid of dummy
      !        atoms.
      !
      type(param_t), intent(in) :: params
      type(el_atom_t), dimension(:), allocatable, intent(out) :: grid
      real(kind=real32) :: step_x,step_y,shift_x,shift_y,charge
      integer(kind=int32) :: N,Nx,Ny,i,j,counter
      !
      ! Nx - number of nodes (dummy atoms) in x-direction on the
      !      electrode.
      !
      ! Ny - number of nodes (dummy atoms) in y-direction on the
      !      electrode.
      !
      ! To compute Nx and Ny one should have access to the variables
      ! of params that contains the node density in each directions.
      ! These variables are el_grid_dens_x and el_grid_dens_y. They
      ! are read from the input file and their dimenstion is 1/A. That
      ! shows how many nodes are defined on the electrode on each
      ! direction
      !
      Nx=int(params%el_grid_dens_x*params%box(1))
      Ny=int(params%el_grid_dens_y*params%box(2))
      !
      ! N - the total number of nodes in the system. Because we have
      !     two electrodes the product of Nx*Ny which gives the number
      !     of nodes (dummy atoms) on each electrode, have to be
      !     multiplied by two.
      ! 
      N=Nx*Ny*2
      
      allocate(grid(N))
      !
      ! step_x - defines the distance between two nodes in x-direction.
      !
      ! step_y - defines the distance between two nodes in y-direction.
      !
      step_x=params%box(1)/Nx
      step_y=params%box(2)/Ny
      !
      ! Usually when the grid nodes (the position of the dummy atoms)
      ! are generating the first grid node in each direction is always
      ! at x = 0 and y = 0. Which means that in most of the cases the
      ! last node in x or y will not lie on x = box(1) or y=box(2). To
      ! center the dummy atoms the fit have to be fit to be equidistant
      ! with respect to the box borders in x and y. So some shift have
      ! to be implemented. The shift in x and y can be computed as:
      !
      shift_x=step_x/2
      shift_y=step_y/2
      !
      ! Charge is assigned separately to the dummy atoms on the left
      ! and right electrodes (because the planes may have different
      ! total charge and sign of the charge).
      !
      ! Assigning the charge and coordinates to the dummy atoms defined
      ! left electrode (left is the electrode at z = 0):
      ! 
      ! Note: the left_el_charge of params reads the total charge
      !       assigned to the left electrode in electron charge units.
      !       Which means that the charge is dimensionless and it is
      !       equivalent of the number of electrones placed on the
      !       plane of the electrode. left_el_charge also has a sign.
      !
      charge=params%left_el_charge/Nx/Ny

      counter=1

      do i=0,Nx-1
         do j=0,Ny-1
            grid(counter)%x=shift_x+i*step_x
            grid(counter)%y=shift_y+j*step_y
            grid(counter)%z=0.0_real32
            grid(counter)%charge=charge
            counter=counter+1
         end do
      end do
      !
      ! Assigning the charge and coordinates to the dummy atoms defined
      ! right electrode (right is the electrode at z = box(3)):
      ! 
      ! Note: the right_el_charge of params reads the total charge
      !       assigned to the right electrode in electron charge units.
      !       Which means that the charge is dimensionless and it is
      !       equivalent of the number of electrones placed on the
      !       plane of the electrode. right_el_charge also has a sign.
      !
      charge=params%right_el_charge/Nx/Ny

      grid(counter:N)%x=grid(1:counter-1)%x
      grid(counter:N)%y=grid(1:counter-1)%y
      grid(counter:N)%z=params%box(3)
      grid(counter:N)%charge=charge

      end subroutine generate_electrode_grid

      end module mod_electrode
