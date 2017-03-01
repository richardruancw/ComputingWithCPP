/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */

#include <fstream>


#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include "Graph.hpp"


// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
using GraphType = Graph<bool,char>;  
using NodeType  = typename GraphType::node_type;

struct poisson_matrix {

	/** Consturct the matrix for a graph.
	 * @param[in] graph The graph used interneally 
	 * @pre For all edges in graph, edge.node1() != edge.node2()
	 */
	poisson_matrix(GraphType* graph):graph_(graph) {}

	// Helper function to perform multiplication.
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn &v, VectorOut &w, Assign) const {
		assert(int(size(v)) == graph_->num_nodes());
		for (auto ni = graph_->node_begin(); ni != graph_->node_end(); ++ni) {
			double sum = 0;
			NodeType n = *ni;
			if (n.value()) {
				sum += v[n.index()];
			} else {
				sum -= n.degree() * v[n.index()];
			}
			for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
				NodeType other_node = (*ii).node2();
				assert(n != other_node);
				if (!n.value() && !other_node.value()) {
					sum += v[other_node.index()];
				}
			}
			Assign::apply(w[n.index()], sum);
		}
		return;
	}

	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<poisson_matrix, Vector> operator*(const Vector &v) const {
		return {*this, v};
	}

	size_t nrow() const{
		return graph_->num_nodes();
	}

private:
	GraphType* graph_;
};

inline std::size_t size(const poisson_matrix &A) {return A.nrow() * A.nrow();}
inline std::size_t num_rows(const poisson_matrix &A) {return A.nrow();}
inline std::size_t num_cols(const poisson_matrix &A) {return A.nrow();}

namespace mtl {
	namespace ashape {
		// Define IdentityMatrix to be a non-scalar type.
		template<>
		struct ashape_aux<poisson_matrix> {
			typedef nonscal type;
		};
	} // end namespace ashape

	// IdentityMatrix inplements the Collection concept with value_type and size_type
	template<>
	struct Collection<poisson_matrix> {
		typedef double value_type;
		typedef unsigned size_type;
	};
} // end namespace mtl



/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  auto ni = g.node_begin();
  while (ni != g.node_end()) {
  	if (bb.contains((*ni).position())) {
  		ni = g.remove_node(ni);
  	} else {
  		++ni;
  	}
  }
  return;
}

/** Set the node value of the boundary nodes to be true. 
 * @param[in] g The Graph whose nodes' value to be set
 *
 */
void set_boundary_flag(GraphType &g) {
	auto bb = Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2,1));
	for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
		Point p = (*ni).position();
		if (norm_inf(p) == 1 || 
			norm_inf(p - Point(0.6, 0.6, 0)) < 0.2  || norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2  ||
			norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 || norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2 ||
			bb.contains(p)) {
			(*ni).value() = true;
		}
	}
}


/** The function used in modeling poission equation*/
struct f_function  {
	double operator()(NodeType n) const{
		return (5 * cos(norm_1(n.position())));
	}
};
/** The function used in modeling poission equation
 * @pre: the n.value() = true;
 */
struct g_function {
	double operator()(NodeType n) const{
		assert(n.value() == true);
		Point p = n.position();
		if (norm_inf(p) == 1) { 
			return 0;
		} else if (norm_inf(p - Point(0.6, 0.6, 0)) < 0.2  || norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2  ||
					norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 || norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2) {
			return -0.2;
		} else {
			return 1;
		}
	}
};

/** Return updated b vector for possion equation given the matrix(graph) and forces*/
template <typename B, typename F, typename G>
void construct_b(GraphType & graph, B & b, const F &f, const G &g, size_t h) {
  	for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
	  	NodeType n =  *ni;
	  	if (n.value()) {
	  		b[n.index()] += g(n);
	  	} else {
	  		b[n.index()] += h * h * f(n);
	  		for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
	  			NodeType other_node = (*ii).node2();
	  			if (other_node.value()) {
	  				b[n.index()] -= g(other_node);
	  			}
	  		}
	  	}
  	}
}

/** Return the associated color for each node given its value and the z-normalizing coefficient*/
template<typename vector_type, typename data_value_type>
  struct MyNodeColor {
  	MyNodeColor(const vector_type &x, const data_value_type &max_value, const data_value_type &min_value):
  	x_(x), max_value_(max_value), min_value_(min_value) {}

   CME212::Color operator()(const NodeType& node) const {
   		if (min_value_ == max_value_) {
   			return CME212::Color::make_heat(0);
   		} else {
   			return (CME212::Color::make_heat((x_[node.index()] - min_value_) / (max_value_ - min_value_)));
   		}
    }
	private:
		const vector_type &x_;
		const data_value_type &max_value_;
		const data_value_type &min_value_;
};

/** Return Point(x, y f(x, y)) given (x, y, 0)*/
template<typename vector_type>
struct NodePosition {
	NodePosition(const vector_type &x): x_(x) {}
	Point operator()(const NodeType &n) const {
		return Point(n.position().x, n.position().y, x_[n.index()]);
	}
private:
	const vector_type x_;
};

/** Plot the current solution and print the residual during the solving process*/
template<typename vector_type, typename viewer_type, typename node_map_type>
class visual_iteration: public itl::cyclic_iteration<double> {
	typedef itl::cyclic_iteration<double> super;
	typedef vector_type Vector;
public:

	visual_iteration(const GraphType &graph, const vector_type &x, viewer_type &viewer, node_map_type &node_map, 
					 const Vector &r0, int max_iter, double tol, double atol = 0, int cycle = 100)
	: super(r0, max_iter, tol, atol, cycle), graph_(graph), x_(x), viewer_(viewer), node_map_(node_map) {}

	bool finished() {
		return super::finished();
	}

	template<typename T>
	bool finished(const T& r) {
		update_current_solution_view();
		return super::finished(r);
	}

private:
	void update_current_solution_view() {
		double max_value, min_value;
	    max_value = *(std::max_element(x_.begin(), x_.end()));
	    min_value = *(std::min_element(x_.begin(), x_.end()));
	    // Center the view and enter the event loop for interactivity
	    MyNodeColor<vector_type, double> headColor(x_, max_value, min_value);
	    NodePosition<vector_type> resultPosition(x_);
		//viewer_.clear();
		//node_map_.clear();
		viewer_.add_nodes(graph_.node_begin(), graph_.node_end(), headColor, resultPosition, node_map_);
		viewer_.add_edges(graph_.edge_begin(), graph_.edge_end(), node_map_);
	}
	const GraphType &graph_;
	const vector_type &x_;
	viewer_type &viewer_;
	node_map_type & node_map_;
};


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

	// Create a nodes_file from the first input argument
	std::ifstream nodes_file(argv[1]);
	// Interpret each line of the nodes_file as a 3D Point and add to the Graph
	std::vector<NodeType> node_vec;
	Point p;
	while (CME212::getline_parsed(nodes_file, p)) {
	node_vec.push_back(graph.add_node(2*p - Point(1,1,0), false));
	}


    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }



  // Get the edge length, should be the same for each edge

  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());
  // Make holes in our Graph

  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));


  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  set_boundary_flag(graph);

  size_t dim = graph.num_nodes();
  mtl::dense_vector<double> b(dim), x(dim);
  b = 0;
  f_function f;
  g_function g;
  construct_b(graph, b, f, g, h);

  poisson_matrix A(&graph);
  itl::pc::identity<poisson_matrix> P(A);
  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  visual_iteration<decltype(x), CME212::SFML_Viewer, decltype(node_map)> 
  visual_iter(graph, x, viewer, node_map, b, 10000, 1.e-10, 0, 50);
  
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]() {
  	itl::cg(A, x, b, P, visual_iter);
  });
  
  viewer.event_loop();
  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
