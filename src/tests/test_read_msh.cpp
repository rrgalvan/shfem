#include <hpfem.hpp>

int main() {
  hpfem::Mesh m;
  m.read_file_msh("squared-mesh-2x2.msh");
  m.print();
  
  for(unsigned i=0; i<m.get_nt(); ++i) {
    std::cout << "Area of triangle " << i << ": "
	      << m.area(i) << std::endl; // ASSERT m.area(i)==0.5
  }
}
