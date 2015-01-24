fe.set_mesh(m);

//   // Get the stored mesh and assert it matches the original mesh
//   const Mesh& m2 = fe.get_mesh();
//   BOOST_CHECK_EQUAL(m.get_ncel(), m2.get_ncel());
//   BOOST_CHECK_EQUAL(m.get_nver(), m2.get_nver());
// }

BOOST_AUTO_TEST_CASE(QuadRuleInOrder)
{
  int nb_nodes = 3;

  // Define a quadrature rule
  VerticesQuadRule qr;
  BOOST_CHECK_EQUAL(qr.size(), nb_nodes);

  // Define a function on quadrature rule nodes
  shfem::QuadFunction f(nb_nodes);
  for(int i=0; i<nb_nodes; i++) f[i]=1.0;  // Set f={1,1,1}

  // Assure integral on reference triangle of previous function is correct
  real_t integral = qr.integrate_on_ref_element(f);
  real_t area_of_ref_triangle = 0.5; // Area of reference triangle
  BOOST_CHECK_CLOSE(integral, area_of_ref_triangle, 1.e-20);
}

BOOST_AUTO_TEST_CASE(FiniteElementLoopInOrder)
{
  // Define mesh
  Mesh mesh;
  mesh.read_file_msh("squared-mesh-2x2.msh");

  // Define quadrature rule
  VerticesQuadRule qr;

  // // Define a Finite Element space and attach previous mesh and QuadRule to it
  // CG_FESpace fe_space;
  // BOOST_CHECK_EQUAL(fe_space.get_nelt(), 0);

  // fe_space.set_mesh(mesh);
  // BOOST_CHECK_EQUAL(fe_space.get_nelt(), mesh.get_ncel());
  // fe_space.set_default_quadrature_rule(qr);

  shfem::FiniteElement fe;
  fe.set_quadrature_rule(qr);
  for(size_t idx_cel=0; idx_cel<mesh.get_ncel(); ++idx_cel) {
    fe.reinit(mesh, idx_cel); // Compute element-specific data
  }
}
