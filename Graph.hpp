#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->node_list[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

//--comment
//--This implementation doesn't satisfy the specs! Two nodes from different graphs
//--that have the same point in their respective graphs would here be returned
//--as equal nodes, when they're not even part of the same graph 
//--Also the specs specifically say that equal nodes have the same *index*, no necessarily
//--the same point. In fact, our specs overall allo for multiple (distinct) nodes that have
//--the same point. These should be considered different nodes.
//--START
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return norm_1(graph_->node_list[uid_]) == norm_1(n.graph_->node_list[n.uid_]);
    }
//--END
    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
//--comment
//--your global ordering here is *technically* correct, but only because your == op
//--is incorrect. Once you fix op==, this < implementation will no longer obey
//--trichotomy (two nodes x and y with same index but different graphs will not return
//--true for *any* of  x==y, x<y, y<x). Also, you'll have some weirdness of having 
//--op== be implemented iwth index and < with points. Rethink this function before
//--the next assignment! Come to OH if you want to chat in more detail!
//--START 
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      return norm_1(graph_->node_list[uid_]) < norm_1(n.graph_->node_list[n.uid_]);
    }
//--END
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_list.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    node_list.push_back(position);
    return Node(this, node_list.size() - 1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type pos = graph_->edge_list[uid_].first;
      return Node(graph_, pos);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type pos = graph_->edge_list[uid_].second;
      return Node(graph_, pos);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.graph_ != graph_) {
        return false;
      } else  {
//--comment
//--In op< you do a good job of testing all posibilities of first/second uid being min/max,
//--but here, you'll return that edge (a ,b) !== (b,a), which violates the spec.
//--START
        return graph_->edge_list[uid_].first == e.graph_->edge_list[e.uid_].first && 
        graph_->edge_list[uid_].second == e.graph_->edge_list[e.uid_].second;
      }
    }
//--END
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
	size_type min_a = std::min(graph_->edge_list[uid_].first, graph_->edge_list[uid_].first);
	size_type max_a = std::max(graph_->edge_list[uid_].first, graph_->edge_list[uid_].first);
	size_type min_b = std::min(e.graph_->edge_list[e.uid_].first, e.graph_->edge_list[e.uid_].second);
	size_type max_b = std::max(e.graph_->edge_list[e.uid_].first, e.graph_->edge_list[e.uid_].second);
	if (min_a < min_b) {
		return true;
	} else if (min_a == min_b && max_a < max_b) {
		return true;
	} else {
		return false;
	}			
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    Edge(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type i = 0; i < edge_list.size(); ++i) {
        const auto &temp = edge_list[i];
        if ((a.uid_ == temp.first && b.uid_ == temp.second)
         || (a.uid_ == temp.second && b.uid_ == temp.first)) {
          return true;
        }
    }   
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
//--comment
//--Duplicate code. Why not just use the has_edge() function?
//--START
    for (size_type i = 0; i < edge_list.size(); ++i) {
        const auto &temp = edge_list[i];
        if ((a.uid_ == temp.first && b.uid_ == temp.second)
         || (a.uid_ == temp.second && b.uid_ == temp.first)) {
//--END
//--comment
//--By having only an edge id in your Edge class and not two node ids, you're unable to 
//--meet the specification in this function that edge.node1()==a and edge.node2()==b
//--In order to meet all the specs, you should consider another implementation! 
//--(You'll probably have to end up storing two node ids. But if you think of a more
//--creative solution, let me know! I'd be really curious to hear it!)
//--START
          return edge(i);
//--END
        }
    } 
    edge_list.push_back(std::make_pair(a.uid_, b.uid_));
    return edge(edge_list.size() - 1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_list.clear();
    node_list.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> node_list;
  std::vector<std::pair<size_type, size_type>> edge_list;

};

#endif // CME212_GRAPH_HPP

//--grade8
