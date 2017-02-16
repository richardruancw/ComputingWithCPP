/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include <queue>

#include "Graph.hpp"

// Define our types
using GraphType = Graph<float>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;


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


//--comment
//--use min_element & lambda function
//--START
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
//--END

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


void shortest_path_lengths(GraphType& g, NodeType& root)
{
  // HW1 #3: YOUR CODE HERE
  //--comment
  //--default vaule is -1 for unreachable nodes
  //--START
  // initialize all node's value to 0
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
    (*ni).value() = 0;
  }
  //--END

  // Run BFS
  std::queue<NodeType> nodeQueue;
  root.value() = 0.001;
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
      if (temp.value() == 0) {
        temp.value() = curr.value() + 1;
        nodeQueue.push(temp);
      }
    }
  }
  float maximum = 0;
   for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
   maximum = std::max((*ni).value(), maximum);
  } 
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
    (*ni).value() = (*ni).value() / maximum;
  }   
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

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

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer

  Point pp(-1, 0, 1);
  NodeType nearest = *(nearest_node(graph, pp));

  shortest_path_lengths(graph, nearest);
  // Center the view and enter the event loop for interactivity

  struct MyNodeColor {
   CME212::Color operator()(const NodeType& node) const {
      return (CME212::Color::make_heat(node.value()));
    }
  };

  MyNodeColor headColor;

  auto node_map = viewer.empty_node_map(graph); 
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map); 
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), headColor, node_map);

  viewer.center_view();
  viewer.event_loop();

  return 0;
}
