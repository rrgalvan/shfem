!> Fortran module containing the definition of all types of mesh are defined
!! @author J. Rafael Rodríguez Galván
!! @author Estefanía Alonso Alonso

module mesh
  use nrtypes
  implicit none
  private
  public :: mesh1d, init, ncells, nvertices

  !> Data type for 1-d meshes
  type, public mesh1d
     !> @brief Array which stores the indexes of vertices, for each cell
     !!
     !! Each row, i, stores the global indices of the vertices
     !! contained in the cell i.  I.e. vertices(i,j) = Global index of
     !! the j-th vertex of cell i
     integer(long), allocatable :: vertices(:,:)

     !> @brief Array which stores the coordinates of each vertex
     !!
     !! Each row, i, stores the coordinates of the vertex whose global
     !! index is equal to i
     real(dp), allocatable :: coordinates(:,:)
  end type mesh1d

contains

  !> Builds a 1D mesh in a given domain (in an interval)
  subroutine init(Th, x1, x2, ncells)
    type(mesh1d), intent(out) :: Th     !< Mesh to be built
    real(dp), intent(in) :: x1          !< Left extreme of the domain (interval [x1,x2])
    real(dp), intent(in) :: x2          !< Right extreme of the domain (interval [x1,x2])
    integer(long), intent(in) :: ncells !< Number of cells (subintervals)

    integer(long), parameter :: dimension_of_affine_space = 1
    integer(long), parameter :: nvertices_in_each_element = dimension_of_affine_space+1 ! Number of vertices / elements
    integer(long) :: nvertices, i
    real(dp) :: h

    ! Define vertices for each cell
    allocate( Th%vertices(ncells,nvertices_in_each_element) )
    forall i=1, ncells
       Th%vertices(i,:) = [i,i+1]  ! Define vertices in cell i
    end do

    ! Define coordinates (1d) for each vertex
    nvertices = ncells+1 ! Total number of vertices in this mesh
    h = (x2-x1)/ncells
    allocate( Th%coordinates(nvertices,dimension_of_affine_space) )
    Th%coordinates(:,dimension_of_affine_space) = [ (x1 + (i-1)*h,  i=1,nvertices) ]
  end subroutine init

  !> Return the number of cells in Th
  integer(long) function ncells(Th)
    type(mesh1d), intent(in) :: Th !< 1D mesh
    ncells = size(Th%vertices, 1)
  end function ncells

  !> Return the number of vertices in Th
  integer(long) function nvertices(Th)
    type(mesh1d), intent(in) :: Th !< 1D mesh
    nvertices = size(Th%coordinates, 1)
  end function nvertices

end module mesh
