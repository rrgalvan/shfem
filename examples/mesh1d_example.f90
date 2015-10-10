program mesh1d_example

  use nrtypes
  use mesh, only : mesh1d, init, ncells, nvertices, get_vertices, get_coordinates
  implicit none

  integer(long) :: i, nc, nv
  type(mesh1d) :: Th

  call init(Th, x1=0.0_dp, x2=1.0_dp, ncells=10)

  nc = ncells(Th)
  print *, "Number of cells in Th:", nc
  do i=1,nc
     print *, "Cell", i, ": vertices...", get_vertices(Th,i)
  end do

  nv = nvertices(Th)
  print *, "Number of vertices in Th:", nv
  do i=1,nv
     print *, "Vertex", i, ": coordinates...", get_coordinates(Th,i)
  end do

end program mesh1d_example
