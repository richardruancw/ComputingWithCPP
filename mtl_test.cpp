/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix {

	IdentityMatrix(size_t n): n_(n) {}
	/** Compute the product ot a vector with this identity matrix
	 */

	// Helper function to perform multiplication.
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn &v, VectorOut &w, Assign) const {
		assert(int(size(v)) == n_);
		for(size_t i = 0; i < n_; ++i) {
			Assign::apply(w[i], v[i]);
		}
	}

	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector> operator*(const Vector &v) const {
		return {*this, v};
	}
	size_t n_;
};

inline std::size_t size(const IdentityMatrix &A) {return A.n_ * A.n_;}
inline std::size_t num_rows(const IdentityMatrix &A) {return A.n_;}
inline std::size_t num_cols(const IdentityMatrix &A) {return A.n_;}

namespace mtl {
	namespace ashape {
		// Define IdentityMatrix to be a non-scalar type.
		template<>
		struct ashape_aux<IdentityMatrix> {
			typedef nonscal type;
		};
	} // end namespace ashape

	// IdentityMatrix inplements the Collection concept with value_type and size_type
	template<>
	struct Collection<IdentityMatrix> {
		typedef double value_type;
		typedef unsigned size_type;
	};
} // end namespace mtl


int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
	const size_t size = 6;

	IdentityMatrix I(size);
	mtl::dense_vector<double> b(size), x(size);
	std::iota(b.begin(), b.end(), 0);
	
	itl::pc::identity<IdentityMatrix> P(I);
	itl::noisy_iteration<double> iter(b, 500, 1.e-6);
	itl::cg(I, x, b, P, iter);
	std::for_each(x.begin(), x.end(), [] (double t) {std::cout<< t << "; ";});

  return 0;
}
