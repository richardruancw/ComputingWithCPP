/**
 * @file Combine.cpp
 * Test script for viewing a subgraph from our Graph combing the plots in @file shortest_path.cpp
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include <queue>

#include "Graph.hpp"

// Define our types
using GraphType = Graph<float>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>>
{
 public:
  // Get all of the iterator traits and make them our own
  using value_type        = typename std::iterator_traits<It>::value_type;
  using pointer           = typename std::iterator_traits<It>::pointer;
  using reference         = typename std::iterator_traits<It>::reference;
  using difference_type   = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last) {
    // HW1 #4: YOUR CODE HERE
      	it_ = first;
      	end_ = last;
      	p_ = p;
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // value_type operator*() const;
  // filter_iterator& operator++();
  // bool operator==(const self_type&) const;

    /** Return the dereference of the iterator
   */
  value_type operator*() const {
  	return *it_;
  }

  /** Increment the filter_iterator forward
   * @post if iter == end() then iter++ == end()
   */
  filter_iterator& operator++() {
  	while(true) {
  		++it_;
  		if (p_(*it_) || it_ == end_) break;
  	}
  	return *this;
  }

  // Return true if two filter_iterator equal
  bool operator==(const filter_iterator& other_iter) const {
  	return it_ == other_iter.it_ && end_ == other_iter.end_;
  }

  // Return the begin of a filter_iterator
  filter_iterator begin() {
  	return *this;
  }

  // Return the end of a filter_iterator
  filter_iterator end() {
  	return filter_iterator(p_, end_, end_);
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/** Return true when the node.value() > than 0.222. This predicate functor
 *  will discard the nodes near the root points
 */
template <typename NODE>
struct SlicePredicate {
  bool operator()(const NODE& n) {
    return n.value() > 0.222;
  }
};

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */


NodeIter nearest_node(const GraphType& g, const Point& point)
{
  // HW1 #3: YOUR CODE HERE
  double small_distance = norm((*g.node_begin()).position() - point);
  auto res = g.node_begin();
  for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
    double temp = norm((*iter).position() - point);
    if (temp < small_distance) {
      res = iter;
      small_distance = temp;
    }
  }
  return res;
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */


int shortest_path_lengths(GraphType& g, NodeType& root)
{
  	// HW1 #3: YOUR CODE HERE
  	// initialize all node's value to 0
  	for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
    	(*ni).value() = -1;
  	}
  	float maximum = 0;
  	// Run BFS
  	std::queue<NodeType> nodeQueue;
  	root.value() = 0;
  	nodeQueue.push(root);
  
  	while(!nodeQueue.empty()) {
    NodeType& curr = nodeQueue.front();
    nodeQueue.pop();
    for (auto ii = curr.edge_begin(); ii != curr.edge_end(); ++ii) {
      NodeType temp;
      if ((*ii).node1() == curr) {
        temp = (*ii).node2();
      } else {
        temp = (*ii).node1();
      }
      // Don't revisit node
      if (temp.value() < 0) {
        temp.value() = curr.value() + 1;
        maximum = std::max(temp.value(), maximum);
        nodeQueue.push(temp);
      }
    }
  } 
  	return (int)(maximum); 
}

/** Normalize a graph's node's value 
 * @param[in,out] g     Input graph
 * @param[in,out] max The normalizing factor
 * @return The maximum path length found.
 *
 * @post root.value() == (old root.value() / max) if old root.value() > 0
 * @post Graph has modified node values by dividing the normalizing factor to non-negative values.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 */
void normalize_node_value(GraphType& g, int max) {
	float maximum = max + 0.0;
  	for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
  		if ((*ni).value() > 0) {
    		(*ni).value() = (*ni).value() / maximum;
		}
  	}  
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph<float>;
  using NodeType  = typename GraphType::node_type;

  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  
  auto node_map = viewer.empty_node_map(graph); 

  Point pp(10, 10, 10);
  NodeType nearest = *(nearest_node(graph, pp));
  // Update the path length
  int maximum = shortest_path_lengths(graph, nearest);
  // Normalize the value in nodes
  normalize_node_value(graph, maximum);
  // Center the view and enter the event loop for interactivity
  struct MyNodeColor {
   CME212::Color operator()(const NodeType& node) const {
      return (CME212::Color::make_heat(node.value()));
    }
  };
  MyNodeColor headColor;

  SlicePredicate<NodeType> slicePre;
  auto Nodeiter = make_filtered(graph.node_begin(), graph.node_end(), slicePre);

  viewer.add_nodes(Nodeiter.begin(), Nodeiter.end(), node_map); 
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.add_nodes(Nodeiter.begin(), Nodeiter.end(), headColor, node_map);

  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
