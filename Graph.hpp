#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Type of the node value */
  using node_value_type = V;
  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
  	num_of_edges = 0;
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
  class Node : private totally_ordered<Node>{

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

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return the value of this node
     *  @post The value of this node might be changed by user
     */
    node_value_type& value() { 
      return graph_->node_value_list[uid_];
    }

    /** Return a const reference of the value of this node
     */
    const node_value_type& value() const{ 
      return graph_->node_value_list[uid_];
    }

    // Return the degree of this node;
    // Complexity O(1)
    size_type degree() const {
    	return graph_->adj_list[uid_].size();
    }

    // Return the begin of the incident_iterator
    incident_iterator edge_begin() const {
    	return incident_iterator(graph_, uid_, 0);
    }

    // Return the end of the incident_iterator
    incident_iterator edge_end() const {
    	return incident_iterator(graph_, uid_, graph_->adj_list[uid_].size());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (n.graph_ != graph_) {
        return false;
      }
      return n.uid_ == uid_;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      return uid_ < n.uid_;
    }

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
   * Complexity: O(1).
   */
  size_type size() const {
    return node_list.size();
  }

  /** Synonym for size(). 
   * Complexity: O(1)
   */
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE

    node_value_list.push_back(node_value);
    node_list.push_back(position);
    std::vector<size_type> temp;
    adj_list.push_back(temp);
    return Node(this, node_list.size() - 1);        
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, nid_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, nid_b_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.graph_ != graph_) {
        return false;
      } else  {
        return (nid_a_ == e.nid_a_ && 
        nid_b_ == e.nid_b_) ||
        (nid_a_ == e.nid_b_ && 
        nid_b_ == e.nid_a_) ;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
  size_type min_a = std::min(nid_a_, nid_b_);
  size_type max_a = std::max(nid_a_, nid_b_);
  size_type min_b = std::min(e.nid_a_, e.nid_b_);
  size_type max_b = std::max(e.nid_a_, e.nid_b_);
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
    size_type nid_a_;
    size_type nid_b_;
    /*
    Edge(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    */

    Edge(const Graph* graph, size_type nid_a, size_type nid_b) {
      graph_ = const_cast<Graph*>(graph);
      nid_a_ = std::min(nid_a, nid_b);
      nid_b_ = std::max(nid_a, nid_b);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_of_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type edge_count = 0;
    for (auto ei = this->edge_begin(); ei != this>edge_end(); ++ei) {
    	if (edge_count == i) {
    		return *ei;
    	} else {
    		edge_count++;
    	}
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(degree)
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Use adjacent list to check the exisitence of edge

    if (a.graph_ != b.graph_ || a.uid_ >= adj_list.size() || b.uid_ >= adj_list.size()) {
      return false;
    } else {
      for (const auto &id : adj_list[a.uid_]) {
          if (id == b.uid_) {
            return true;
          }
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

    if (has_edge(a, b)) {
        size_type left = std::min(a.uid_, b.uid_);
        size_type right = std::max(a.uid_, b.uid_);
        return Edge(this, left, right);
    } else {
      num_of_edges++;
      adj_list[a.uid_].push_back(b.uid_);
      if (a.uid_ != b.uid_) {
      	adj_list[b.uid_].push_back(a.uid_);
      }
      return Edge(this, a.uid_, b.uid_); 
    }
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    adj_list.clear();
    node_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    NodeIterator(const Graph* graph, size_type pos) {
    	graph_ = const_cast<Graph*> (graph);
    	position_ = pos;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const


    /** Dereference operator for this iterator
     *  @ pre The nodes of the graph should not be changed after the creation of
     *  this iterator.
     */
    Node operator*() const {
    	return Node(graph_, position_);
    }

    /** Forward operator for this iterator
     *  @pre The iterator != g.node_end()
     *  @ post @a position increase by one
     */
    NodeIterator& operator++() {
    	//std::cout << "NodeIterator++ called!!" << std::endl;
    	position_++;
    	return *this;
    }

    // Return true if two NodeIterator equalsS
    bool operator==(const NodeIterator& inputIter) const {
    	return inputIter.graph_ == graph_ && inputIter.position_ == position_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    Graph* graph_;
    size_type position_;
  };

  // HW1 #2: YOUR CODE HERE
    // Return the begin of node iterator
     node_iterator node_begin() {
     	return NodeIterator(this, size_type(0));
     }
    // Return the begin of node iterator
     node_iterator node_begin() const{
     	return NodeIterator(this, size_type(0));
     }

    // Return the end of node iterator
     node_iterator node_end() {
     	return NodeIterator(this, size_type(node_list.size()));
     }

    // Return the end of node iterator
     node_iterator node_end() const{
     	return NodeIterator(this, size_type(node_list.size()));
     }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>{
    public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    IncidentIterator(const Graph* graph, size_type root_count, size_type end_count) {
    	graph_ = const_cast<Graph*>(graph);
    	root_count_ = root_count;
    	end_count_ = end_count;
    }

    // HW1 #3: YOUR CODE HERE

    /** Return the dereferenced edge object
     *  @pre the incident iterator is in the range of edge_begin() and edge_end()
     */
    Edge operator*() {
    	return Edge(graph_, root_count_, graph_->adj_list[root_count_][end_count_]);
    }

    /** Increment the incident iterator 
     *  @post new iter = old iter++
     */
    IncidentIterator& operator++() {
    	end_count_++;
    	return *this;
    }

    /** Return true if two incident iterator equal
     */
    bool operator==(const IncidentIterator other_iter) const{
    	return graph_ == other_iter.graph_ && root_count_ == other_iter.root_count_ &&
    		   end_count_ == other_iter.end_count_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type root_count_;
    size_type end_count_;
  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    EdgeIterator(const Graph* graph, size_type position,
    			size_type root_count, size_type end_count) {

    	graph_ = const_cast<Graph*> (graph);
    	root_count_ = root_count;
    	end_count_ = end_count;
    	position_ = position;

    }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Derefrence operator for the iterator
     * @pre The graph strucutre should not be modified after the creation of iterator and
     * before using the iteraot.
     * @post The graph does not change
     */
    Edge operator*() const{
    	return Edge(graph_, root_count_, graph_->adj_list[root_count_][end_count_]);
    }


    /** Increment EdgeIterator forward
     *  @pre This iterator is in the range of graph.edge_begin() and graph.edge_end()
     *  @post new edgeIterator = old ++edgeIterator;
     */
    EdgeIterator& operator++() { 
    	bool jobDone = false;
    	// Have not yet reached the end of the adj_list in the next step
    	while (root_count_ < graph_->adj_list.size()) {
    		if (!jobDone) {
		    	if (graph_->adj_list[root_count_].size() == 0 ) {
		    		root_count_++;
		    		end_count_ = 0;
		    	} else if(end_count_ < graph_->adj_list[root_count_].size() - 1) {
		    		end_count_++;
		    	} else {
		    		root_count_++;
		    		end_count_ = 0;
		    	}
		    	jobDone = true;
	    	} else {
	    		if (graph_->adj_list[root_count_].size() == 0) {
	    			root_count_++;
	    			end_count_ = 0;
	    			continue;
	    		}
	    		if (root_count_ <= graph_->adj_list[root_count_][end_count_]) {
	    			break;
	    		}
	    		if (end_count_ < graph_->adj_list[root_count_].size() - 1) {
	    			end_count_++;
	    		} else {
	    			root_count_++;
	    			end_count_ = 0;
	    		}
	    	}
    	}
    	position_++;
    	return *this;
    }

    /** Test whether two EdgeIterator equal
     *  @pre the graph structure have not been changed
     *  @post the graph does not change
     */
    bool operator==(const EdgeIterator& other_iter) const {
    	return graph_ == other_iter.graph_ && position_ == other_iter.position_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type position_;
    size_type root_count_;
    size_type end_count_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Return the begin of edge_iterator
   
  edge_iterator edge_begin() const {
  	size_type i = 0;
  	while (i < adj_list.size()) {
  		if (adj_list[i].size() != 0) {
  			break;
  		}
  		i++;
  	}
  	return EdgeIterator(this, 0, i, 0);
  }
  
  // Return the end of the edge_iterator
  edge_iterator edge_end() const {
  	return EdgeIterator(this, num_of_edges, adj_list.size(), adj_list[adj_list.size() - 1].size());
  }
  private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  // The list of node
  std::vector<Point> node_list;
  // The list of node value
  std::vector<node_value_type> node_value_list;
  // adjacent list
  std::vector<std::vector<size_type>> adj_list;

  size_type num_of_edges;
};

#endif // CME212_GRAPH_HPP
