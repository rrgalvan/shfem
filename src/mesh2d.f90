!> Fortran module containing the definition of all types of mesh are defined
!! @author J. Rafael Rodríguez Galván
!! @author Estefanía Alonso Alonso

module mesh2d_module
  use nrtypes
  implicit none
  private
  public :: triangle_mesh2d, triangle_mesh2d_init_from_file, &
       triangle_mesh2d_ncells, triangle_mesh2d_nvertices, &
       triangle_mesh2d_get_vertices, triangle_mesh2d_get_coordinates

  !> Data type for 2-d meshes which are composed of triangles
  type triangle_mesh2d

     !> @brief Array which stores the indexes of vertices, for each cell
     !!
     !! Each column, j, stores the global indices of the vertices
     !! contained in the cell j.  I.e. vertices(i,j) = Global index of
     !! the i-th vertex of cell j
     integer(long), allocatable :: vertices(:,:)

     !> @brief Array which stores the coordinates of each vertex
     !!
     !! Each column, j, stores the coordinates of the vertex whose global
     !! index is equal to j
     real(dp), allocatable :: coordinates(:,:)

    integer(long) :: ncells !< Number of cells
    integer(long) :: nvertices !< Number of vertices in the mesh

 end type triangle_mesh2d

contains

  !> Builds a triangle 2d mesh which is stored in a .msh (see FreeFem++ doc)
  subroutine triangle_mesh2d_init_from_file(Th, file_name)
    type(triangle_mesh2d), intent(out) :: Th  !< Mesh to be built
    character (len=*), intent(in) :: file_name      !< File to read mesh from (.msh format)

    integer(long), parameter :: dimension_of_affine_space = 2
    integer(long), parameter :: nvertices_in_each_element = 3 ! Number of vertices / elements
    integer(long) :: ncells, nvertices, i

    open(unit=10, status='old', file=file_name)
    read (10,*) nvertices, ncells
    ! print *, "DBG: ncells, nvertices =", ncells, nvertices

    ! Define coordinates for each vertex
    Th%nvertices = nvertices
    allocate( Th%coordinates(dimension_of_affine_space, nvertices) )
    do i=1, nvertices
       read (10,*) Th%coordinates(:,i)
    end do
!    print *, "DBG: coordinates =", Th%coordinates

    ! Define vertices for each cell
    Th%ncells = ncells
    allocate( Th%vertices(nvertices_in_each_element, ncells) )
    do i=1, ncells
       read (10,*) Th%vertices(:,i)
    end do
!    print *, "DBG: vertices =", Th%vertices

  end subroutine triangle_mesh2d_init_from_file

  !> Return the number of cells in Th
  integer(long) function triangle_mesh2d_ncells(Th)
    type(triangle_mesh2d), intent(in) :: Th !< 2d mesh composed of triangles
    triangle_mesh2d_ncells = Th%ncells !size(Th%vertices, 2)
  end function triangle_mesh2d_ncells

  !> Return the number of vertices in Th
  integer(long) function triangle_mesh2d_nvertices(Th)
    type(triangle_mesh2d), intent(in) :: Th !< 2d mesh composed of triangles
    triangle_mesh2d_nvertices = Th%nvertices !size(Th%coordinates, 2)
  end function triangle_mesh2d_nvertices

  !> Return a copy of the vertices contained in the i-th cell of Th
  function triangle_mesh2d_get_vertices(Th,i)
    type(triangle_mesh2d), intent(in) :: Th !< 2d mesh composed of triangles
    integer(long), intent(in) :: i
    integer(long), parameter :: nvertices_in_each_element = 3 ! Number of vertices / elements
    integer(long), dimension(nvertices_in_each_element) :: triangle_mesh2d_get_vertices
    triangle_mesh2d_get_vertices = Th%vertices(:,i)
  end function triangle_mesh2d_get_vertices

  !> Return a copy of the coordinates of the i-th vertex of Th
  function triangle_mesh2d_get_coordinates(Th,i)
    type(triangle_mesh2d), intent(in) :: Th !< 2d mesh composed of triangles
    integer(long), intent(in) :: i
    integer(long), parameter :: dimension_of_affine_space = 2
    real(dp), dimension(dimension_of_affine_space) :: triangle_mesh2d_get_coordinates
    triangle_mesh2d_get_coordinates = Th%coordinates(:,i)
  end function triangle_mesh2d_get_coordinates

end module mesh2d_module
