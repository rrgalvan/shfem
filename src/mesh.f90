!> Fortran module containing the definition of all types of mesh are defined
!! @author J. Rafael Rodríguez Galván
!! @author Estefanía Alonso Alonso

module mesh

  use mesh1d_module
  use mesh2d_module

  interface init
     module procedure mesh1d_init
     module procedure triangle_mesh2d_init_from_file
  end interface init

  interface ncells
     module procedure mesh1d_ncells
     module procedure triangle_mesh2d_ncells
  end interface ncells

  interface nvertices
     module procedure mesh1d_nvertices
     module procedure triangle_mesh2d_nvertices
  end interface nvertices

  interface get_vertices
     module procedure mesh1d_get_vertices
     module procedure triangle_mesh2d_get_vertices
  end interface get_vertices

  interface get_coordinates
     module procedure mesh1d_get_coordinates
     module procedure triangle_mesh2d_get_coordinates
  end interface get_coordinates

end module mesh
