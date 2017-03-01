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
template <typename V, typename E>
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Type of the node value */
  using node_value_type = V;
  using edge_value_type = E;
  /** Type of this graph. */
  using graph_type = Graph<V, E>;
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
  /** Type of the uid */
  using uid_type = size_type;

	struct edge_info {
		uid_type uid_other_;
		edge_value_type edge_value_;
		// Default constructor
		edge_info(): uid_other_(), edge_value_() {}
		//
		edge_info(const uid_type &uid_other, const edge_value_type &edge_value = edge_value_type()):
		uid_other_(uid_other), edge_value_(edge_value) {} 
	};

  	struct node_info {
	  	Point p_;
	  	node_value_type node_value_;

	  	size_type index_;
	  	std::vector<edge_info> adj_;
	  	// Default constructor
	  	node_info():p_(), node_value_(), adj_() {}
	  	// Creat new node
	  	node_info(const Point &position, size_type index, const node_value_type &node_value = node_value_type()):
	  	p_(position), index_(index), node_value_(node_value), adj_() {}

  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  /** Construct an empty graph. */
  Graph() {
  	num_of_edges_ = 0;
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
    Point& position() {
    	assert(valid());
    	return graph_->node_list_[uid_].p_;
    }

    const Point& position() const {
    	assert(valid());
    	return graph_->node_list_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type& index() const {
    	assert(valid());
    	return graph_->node_list_[uid_].index_;
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
    	assert(valid());
    	return graph_->node_list_[uid_].node_value_;
    }

    /** Return a const reference of the value of this node
     */
    const node_value_type& value() const{
    	assert(valid()); 
    	return graph_->node_list_[uid_].node_value_;
    }

    // Return the degree of this node;
    // Complexity O(1)
    size_type degree() const {
    	return graph_->node_list_[uid_].adj_.size();
    }

    // Return the begin of the incident_iterator
    incident_iterator edge_begin() const {
    	return incident_iterator(graph_, uid_, 0);
    }

    // Return the end of the incident_iterator
    incident_iterator edge_end() const {
    	return incident_iterator(graph_, uid_, graph_->node_list_[uid_].adj_.size());
    }


    /** Test whether this node and @a n are equal.
     * @return True if two nodes have the same graph and the same index.
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

    bool valid() const {
    	return 0 <= uid_ && uid_ < graph_->node_list_.size() && 
    	graph_->node_list_[uid_].index_ < graph_->index_to_uid_.size() &&
    	graph_->index_to_uid_[graph_->node_list_[uid_].index_] == uid_;
    }
    /** Return a node object given the access to the graph and position
     * @pre 0 <= @a uid < node_list.size()
     * @post the returned node might change the value of the node in graph.
     * Complexity: O(1)
     */
    Node(const Graph* graph, uid_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

    std::vector<edge_info> & get_adj_() {
    	return graph_->node_list_[uid_].adj_;
    }
  };

  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    return index_to_uid_.size();
  }

  /** Synonym for size(). 
   * Complexity: O(1)
   */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @pre The graph did not use the delete node method before.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    node_list_.push_back(node_info(position, num_nodes(), node_value));

    index_to_uid_.push_back(node_list_.size() - 1);


    return Node(this, index_to_uid_[index_to_uid_.size() - 1]);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this && n.valid();
  }
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type index) const {
    // HW0: YOUR CODE HERE
    return Node(this, index_to_uid_[index]);
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

    /** Compute the length of the edge in term of euclidean distance
     * 
     */
    double length() const {
      return norm(graph_->node_list_[nid_a_].p_ - graph_->node_list_[nid_b_].p_);
    }

    /** Return the reference of the value associated with this edge
     *  @pre The edge is valid
     *  @post if edge_a == edge_b then edge_a.value() == edge_b.value()
     */
    edge_value_type& value() {
      size_type left = std::min(nid_a_, nid_b_);
      size_type right = std::max(nid_a_, nid_b_);
      for (size_type i = 0; i < graph_->node_list_[left].adj_.size(); ++i) {
        if (graph_->node_list_[left].adj_[i].uid_other_ == right) {
          return graph_->node_list_[left].adj_[i].edge_value_;
        }
      }
    }

    /** Return the reference of the value associated with this edge
     *  @post if edge_a == edge_b then edge_a.value() == edge_b.value()
     */
    const edge_value_type& value() const {
      return value();
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
	    return graph_ < e.graph_;
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
    uid_type nid_a_;
    uid_type nid_b_;

    Edge(const Graph* graph, uid_type nid_a, uid_type nid_b):
    graph_(const_cast<Graph*>(graph)), nid_a_(nid_a), nid_b_(nid_b) { }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_of_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type edge_count = 0;
    for (auto ei = this->edge_begin(); ei != this->edge_end(); ++ei) {
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

    if (a.graph_ != b.graph_ || a.uid_ >= node_list_.size() || b.uid_ >= node_list_.size()) {
      return false;
    } else {
      for (const auto &id : node_list_[a.uid_].adj_) {
          if (id.uid_other_ == b.uid_) {
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
      return Edge(this, a.uid_, b.uid_);
    } else {
      // update edge
      num_of_edges_++;
      node_list_[a.uid_].adj_.push_back(edge_info(b.uid_));

      if (a.uid_ != b.uid_) {
		node_list_[b.uid_].adj_.push_back(edge_info(a.uid_));
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
    node_list_.clear();
    index_to_uid_.clear();
    num_of_edges_ = 0;
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
    	return Node(graph_, graph_->index_to_uid_[position_]);
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
    uid_type uid_;
  };

  // HW1 #2: YOUR CODE HERE
    // Return the begin of node iterator
     node_iterator node_begin() {
     	return NodeIterator(this, uid_type(0));
     }
    // Return the begin of node iterator
     node_iterator node_begin() const{
     	return NodeIterator(this, uid_type(0));
     }

    // Return the end of node iterator
     node_iterator node_end() {
     	return NodeIterator(this, uid_type(index_to_uid_.size()));
     }

    // Return the end of node iterator
     node_iterator node_end() const{
     	return NodeIterator(this, uid_type(index_to_uid_.size()));
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

    IncidentIterator(const Graph* graph, uid_type root_count, size_type end_count):
    graph_(const_cast<Graph*>(graph)), root_count_(root_count), end_count_(end_count) { }

    // HW1 #3: YOUR CODE HERE

    /** Return the dereferenced edge object
     *  @pre The incident iterator is in the range of edge_begin() and edge_end().
     *  @return The returned edge.node1() is the root node, and edge.node2() is the end node.
     */
    Edge operator*() const{
    	return Edge(graph_, root_count_, graph_->node_list_[root_count_].adj_[end_count_].uid_other_);

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
    // The uid of the root node
    uid_type root_count_;
    // The index of the end node in the adjacent list of the root node
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

    /** Return an EdgeIterator given the access to the graph, the position, the position
     *  of the root node and the position of the end node
     * @pre 0 <= @a position < num_edges()
     * @pre 0 <= @a root_count < num_nodes()
     * @pre 0<= @a end_count < number of nodes root links to
     * @post return an EdgeIterator that might change the value of the graph.
     */
    EdgeIterator(const Graph* graph, node_iterator root_it, incident_iterator inci_it) {
    	graph_ = const_cast<Graph*> (graph);
    	root_it_ = root_it;
    	inci_it_ = inci_it;
    }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Derefrence operator for the iterator
     * @pre The graph strucutre should not be modified after the creation of iterator and
     * before using the iteraot.
     * @post The graph does not change
     * Complexity: O(1)
     */
    Edge operator*() const{
    	return (*inci_it_);
    }

    /** Increment EdgeIterator forward
     *  @pre This iterator is in the range of graph.edge_begin() and graph.edge_end()
     *  @post new edgeIterator = old ++edgeIterator;
     *  @post if int count = 0, for(auto e = g.edge_begin(); e != e.edge_end(); count++, ++e)
     *        then, count + 1 == num_edges().
     */
    EdgeIterator& operator++() { 

    	// Have not yet reached the end of the adj_list in the next step
    	while(root_it_ != graph_->node_end()) {
    		auto n = *root_it_;
    		while(inci_it_ != n.edge_end()) {			
    			++inci_it_;
    			if (inci_it_ == n.edge_end()) {
    				break;
    			} else if (n.uid_ <= (*inci_it_).node2().uid_) {
    				return *this;
    			}
    		}
    		++root_it_;
    		if (root_it_ != graph_->node_end()) {inci_it_ = (*root_it_).edge_begin();}
    	}
    	return *this;
    }
    /** Test whether two EdgeIterator equal
     *  @pre the graph structure have not been changed
     *  @post the graph does not change
     */
    bool operator==(const EdgeIterator& other_iter) const {
    	return graph_ == other_iter.graph_ && root_it_ == other_iter.root_it_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Acess to the graph object
    Graph* graph_;
    // Inner uid for the root node
    node_iterator root_it_;
    // Inner index for the end nodes that in the adj_ of root node.
    incident_iterator inci_it_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Return the begin of edge_iterator
   
  edge_iterator edge_begin() const {
  	node_iterator node_it = node_begin();
  	auto inci_it = (*node_it).edge_begin();
  	return EdgeIterator(this, node_it, inci_it);
  }
  
  /** Return the end of the edge_iterator
   *  @pre The num_of_edges_ is valid
   */
  edge_iterator edge_end() const {
  	node_iterator node_it = node_end();
  	incident_iterator inci_it = (*node_begin()).edge_begin();
  	return EdgeIterator(this, node_it, inci_it);
  }



  /** Remove a edge linking to given two nodes
   *  @return Return 1 if the edge was in the graph and has been removed. Return 0 if the edge was
   *  not in the graph.
   *  @post If return 1 then: 
   *        1) old num_edges() = new num_edges() + 1.
   *        2) has_edge(@a n1, @a n2) == false;
   *        3) All other Edge objects but the deleted one remain valid.
   *        4) All outstanding iterator are invalidated
   *  Complexity: O(d), d is the maxium degree of nodes in the graph.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {

  	if (!has_edge(n1, n2)) {
  		return 0;
  	}
  	num_of_edges_--;
  	for (size_type i = 0; i < node_list_[n1.uid_].adj_.size(); ++i) {
  		if (node_list_[n1.uid_].adj_[i].uid_other_ == n2.uid_) {
  			node_list_[n1.uid_].adj_.erase(node_list_[n1.uid_].adj_.begin() + i);
  			break;
  		}
  	}
  	for (size_type i = 0; i < node_list_[n2.uid_].adj_.size(); ++i) {
  		if (node_list_[n2.uid_].adj_[i].uid_other_ == n1.uid_) {
  			node_list_[n2.uid_].adj_.erase(node_list_[n2.uid_].adj_.begin() + i);
  			break;
  		}
  	}
  	return 1;
  }

  /** Remove a given edge. This method calls the remove_edge(n1, n2). The return value and post
   *  conditions are identical to it.
   *  Complexity: O(d), d is the maximum degree of a node in the graph
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  } 

  /** Remove the edge that the edge iterator points to.
   *  @pre The @a e_it is valid
   *  @return An edge_iterator object new_it such that new_it = old @a e_it++) if @a e_it != edge_end()
   *  @post 1) old num_edges() = new num_edges() + 1.
   *        2) has_edge(old *e_it) == false;
   *        3) All other Edge objects but the deleted one remain valid.
   *        4) All outstanding iterator are invalidated
   *  Complexity: O(d), d is the maximum degree of a node in the graph
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    node_iterator root = e_it.root_it_;
    incident_iterator inci = e_it.inci_it_;
    remove_edge(*e_it);
    if (inci == (*root).edge_end()) {++e_it;}
    return e_it;
   }
 
  /** Remove the given node and the edges that link to it
   *  @return If has_node(@a n) return 1 else return 0
   *  @post 
   *		1) old num_nodes() = new num_nodes() + 1.
   *        2) old num_edges() = new num_edges() + old @a n.degree()
   *        3) has_node(@a n) == false;
   *        4) All other node objects but the deleted one remain valid.
   *        5) All other edge objects but those link to node @a n remian valid.
   *        4) All outstanding iterator of the graph are invalidated
   *  Complexity: O(num_nodes()) given a sparse graph
   */
  size_type remove_node(const Node& n) {
  	if (!has_node(n)) {
  		return 0;
  	}
    size_type index = n.index();
    // Adjust the index of "node"
    for (size_type i = index + 1; i < num_nodes(); ++i) {
    	node_list_[index_to_uid_[i]].index_--;
    }
    // Remove edges that link to this node
    std::vector<edge_type> edges_to_remove;
    std::copy(n.edge_begin(), n.edge_end(), std::back_inserter(edges_to_remove));
    std::for_each(edges_to_remove.begin(), edges_to_remove.end(), [this] (edge_type &e) { remove_edge(e); });
    // Remove the nodes 
    index_to_uid_.erase(index_to_uid_.begin() + index);
    return 1;
  } 


  /** Remove the given node and the edges that link to it given the position of the node
   *  @pre @a n_it is a valid interator and @a n_it != node_end().
   *  @return A node_iterator new n_it such that if has_node(* n_it) then new n_it == ++n_it.
   *  @post 
   *		1) old num_nodes() = new num_nodes() + 1.
   *        2) old num_edges() = new num_edges() + old (*n_it).degree()
   *        3) has_node(* old n_it) == false;
   *        4) All other node objects but the deleted one remain valid.
   *        5) All other edge objects but those link to node @a n remian valid.
   *        4) All outstanding iterator of the graph are invalidated
   *  Complexity: O(num_nodes()) given a sparse graph.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  } 



  private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // The list of node
  std::vector<node_info> node_list_;
  //
  size_type num_of_edges_;
  //
  std::vector<size_type> index_to_uid_;

};

#endif // CME212_GRAPH_HPP
