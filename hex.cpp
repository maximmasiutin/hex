/*

ASCII Hexagonal Grid HEX Board Game v1.4 (with Monte-Carlo AI)

Copyright 2020-2025 Maxim Masiutin. All rights reserved.

2025-11-23

This program is distributed under the GNU GPL v2.0

HEX is a two-person abstract strategy board game in which players try to
connect opposite sides of a hexagonal board.

The game rules were first described by the mathematician and poet Piet Hein in
1942.

This program implements Monte-Carlo simulation to get a good move for a
computer player.

Each legal move will be evaluated using a million (1,000,000) of trials. See
the "monte_carlo_iterations" variable to configure this value.

Each trial winds the game forward by randomly selecting successive moves until
there is a winner, and the trial is counted as a win or loss. The ratio:
wins/trials are the AI's metric for picking which next move to make.

The Monte-Carlo simulation is implemented very efficiently, so it takes just
about a second to make a move from a million trials on an average notebook on
a 7x7 field, or about 5 seconds on an 11x11 field.

This is a pure Monte-Carlo implementation, without the min-max algorithm
or alpha-beta pruning.

The computer can play HEX intelligently against a human on an 11 by 11 board,
and it only takes a few seconds for the computer to make the intelligent move.
You can decrease the "monte_carlo_iterations" variable to make the computer
make each move faster.

The program uses a special algorithm to determine who wins (as opposed to the
Dijkstra's), which is based on a high-speed function AiPathExistsBlue() to
check if the blue player has a connection, the winning position. This function
uses a precomputed table of node connections and is only called from
Monte-Carlo simulation when there is no empty area on the board, that is, when
the entire board is completed randomly by red and blue stones. Since there is
no "draw" in the HEX game, once the whole board is full and the "blue" player
is not connected, it means that the "red" player already has a connection.
That's why the function is only implemented for the "blue" player.

The program uses the new C++11 <random> library just to fill the seed data for
the fast xorshift random number generator, which is much quicker than one
supplied in the STL, and is better suitable for the Monte-Carlo. It uses the
xorshift implementation published by George Marsaglia from The Florida State
University, in 2003 in the Journal of Statistical Software, DOI:
10.18637/jss.v008.i14

This implementation:

- draws the board using ASCII symbols, and it has three various ways to draw the
board:
    (1) conventional hexagonal grid in the parallelogram board, inclined right
    (2) non-standard rectangular board
    (3) compact board with just 4x2 ASCII characters per cell

- inputs a move and determines if a move is legal, it has the following options
  by the number of players:
    (0) computer vs computer (demo mode)
    (1) player vs. computer (the computer makes a Monte-Carlo simulation without
        any tree search to get a better move)
    (2) player vs player (two human players)

- determines who won.

To determine who won, the program uses graph (nodes and edges) representation
and treats the game as a path-finding problem. It implements and uses
Dijkstra's algorithm to find the shortest path within a conventional graph.
Each internal node (hexagon) has six neighbors, so each has six edges. The
corners and edges are special. A corner has two or three neighbors, and a
non-corner edge has four neighbors.

The program has a simple text-mode interface for entering a move for the first
player, displaying the board, and then, unless playing with the computer,
entering a move for the second player, and so on, until one player wins.

Therefore, the program can determine when the game is over and then announce
the winner, and it does so immediately.

The program uses a straightforward method to input players' moves. That is, it
lets the player enter the (x, y) coordinates corresponding to the currently
empty cell column (x) and row (y). After the user entered the move, the
program checks whether this is a legal move, and if not, the program asks the
user again to make a move.

*/

#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

// Holds whether a cell is empty (no move done yet) or it is occupied by either
// Red or Blue player
enum class CellValue : uint8_t { Blank, Blue, Red };

// The Field class holds a 2-dimensional array of cells that store players'
// moves
class Field {
public:
  const unsigned int width;
  const unsigned int height;
  Field(const unsigned int a_width, const unsigned int a_height)
      : width(a_width), height(a_height) {}
  using FieldData = std::vector<std::vector<CellValue>>;
  FieldData cells;

  // clear the field where all the players moves are recorded
  void clear(void) {
    cells.clear();
    cells.resize(
        static_cast<FieldData::size_type>(width)); // to suppress a warning
    for (unsigned int x = 0; x < width; ++x) {
      cells[x].resize(
          static_cast<FieldData::size_type>(height)); // to suppress a warning
      for (unsigned int y = 0; y < height; ++y) {
        cells[x][y] = CellValue::Blank;
      }
    }
  }
  bool coord_out_of_range_xy(const unsigned int x, const unsigned int y) const {
    return (x >= width) || (y >= height);
  }

  bool coord_out_of_range_x(const unsigned int x) const { return (x >= width); }

  bool coord_out_of_range_y(const unsigned int y) const {
    return (y >= height);
  }
};

using ConsoleBuffer = std::vector<std::vector<char>>; // used to draw the board

/*********************************************************************

Console draw helper classes

**********************************************************************/

inline bool
is_odd_uint(const unsigned int a) // just returns whether an integer number is
                                  // odd, i.e. is 1, 3, etc...
{
  return (a & 1) != 0u;
}

/* The "incline" only affects how the adjacent hexagons are connected and how
the board is painted. Still, the winning path is always between 0-index and
max-index rows (or columns) regardless of how they are connected or painted,
so if the player connects a 0-index hexagon in a path towards the max-index
hexagon, that player wins. Specifically, if the start and stop and edges of
the board are horizontal, i.e., there is a path from up to down from y=0 to
y=max_height, then the winner is a red player, and vice versa: a way from left
to right, or x=0 to x=max_width from one vertical edge of the board to another
means the blue player wins; there can be no "draw". */

enum class BoardIncline : uint8_t { left, right, rect };

class ConsoleAbstract {
protected:
  const char blank_console_background =
      ' '; // an ASCII character that fills blank empty background space in the
           // console buffer
  const char cell_red = 'O';  // the cell is occupied by the red player
  const char cell_blue = 'X'; // the cell is occupied by the blue player
  const char cell_available_by_default =
      ' '; // the cell is not yet occupied by a player, this is the default
           // characters, but a different one may be used depending on the way
           // to draw the field

  const Field &field;
  ConsoleBuffer console;
  unsigned int console_width;
  unsigned int console_height;

public:
  ConsoleAbstract(const Field &a_field)
      : field(a_field), console_width(0),
        console_height(0) {}; // default constructor

  // returns the character to paint when there was no player's move to the cell
  // in the field, so the cell is available
  virtual char cell_available(void) = 0;

  // convert cell value to an ASCII character to be printed
  char cellvalue2char(CellValue value) {
    switch (value) {
    case CellValue::Red:
      return cell_red;
    case CellValue::Blue:
      return cell_blue;
    default: // case CellValue::Blank:
      return cell_available();
    }
  }

  // clear the ASCII screen console buffer before printing the board
  void clear(void) {
    console.resize(static_cast<ConsoleBuffer::size_type>(
        console_width)); // to suppress a warning
    for (unsigned int x = 0; x < console_width; ++x) {
      console[x].resize(static_cast<ConsoleBuffer::size_type>(
          console_height)); // to suppress a warning
      for (unsigned int y = 0; y < console_height; ++y) {
        console[x][y] = blank_console_background;
      }
    }
  }

  // output the console buffer to std::cout (STDOUT)
  void print(void) const {
    for (unsigned int y = 0; y < console_height; ++y) {
      for (unsigned int x = 0; x < console_width; ++x) {
        std::cout << console[x][y];
      }
      std::cout << std::endl;
    }
  }

  // draw the board foreground: grid and the cell contents
  virtual void draw_field(void) = 0;

  // return whether there is the non-conventional, rectangular field with
  // alternative hexagon adjacency rules
  virtual BoardIncline get_board_incline(void) = 0;
  virtual ~ConsoleAbstract() {};
};

class ConsoleCondensed : public ConsoleAbstract {
protected:
public:
  ConsoleCondensed(const Field &a_field) : ConsoleAbstract(a_field) {
    unsigned int double_height = field.height << 1;  // *(2^1)
    unsigned int quadruple_width = field.width << 2; // *(2^2)
    console_width =
        quadruple_width + double_height - 2; // do not paint the last
    console_height = double_height - 1;      // do not paint the last line
  };

  virtual char cell_available(void) { return '.'; }
  virtual BoardIncline get_board_incline(void) { return BoardIncline::left; }
  virtual void draw_field(void) {
    unsigned int y = 0;
    for (unsigned int row = 0; row < field.height; ++row) {
      // paint first line, with cell values
      unsigned int x =
          row *
          2; // each row starts shifted right comparing to the previous row
      for (unsigned int column = 0; column < field.width; ++column) {
        console[x][y] = cellvalue2char(field.cells[column][row]);
        x += 2;
        if (column + 1 < field.width) // do not paint the last one
        {
          console[x][y] = '-';
          x += 2;
        }
      }

      // paint second line, that is in between the cell values
      if (row + 1 < field.height) // do not paint last line
      {
        unsigned int x = row * 2 + 1;
        y++;
        for (unsigned int column = 0; column < field.width; ++column) {
          console[x][y] = '\\';
          x += 2;
          if (column + 1 < field.width) {
            console[x][y] = '/'; // do not paint the last one
            x += 2;
          }
        }
        y++;
      }
    }
  }
};

class ConsoleConventionalHexagons : public ConsoleAbstract {
protected:
  /*

  scale = 1
   __
  /  \__
  \__/  \
  /  \__/
  \__/  \
     \__/

  scale = 2
    _____
   /     \
  /       \_____
  \       /     \
   \_____/       \
   /     \       /
  /       \_____/
  \       /     \
   \_____/       \
         \       /
          \_____/


  scale = 3
     ________
    /        \
   /          \
  /            \________
  \            /        \
   \          /          \
    \________/            \
    /        \            /
   /          \          /
  /            \________/
  \            /        \
   \          /          \
    \________/            \
             \            /
              \          /
               \________/

  */

  // internal private variables needed for the draw_cell method
  const unsigned int cell_border_overlap = 1;
  unsigned int cell_halfwidth;
  unsigned int cell_height;
  unsigned int cell_halfheight;
  unsigned int edge_width;

  // draw the hexagonal grid with the contents of the field
  // each hexagonal cell is drawn one by one in a nested loop
  virtual void draw_field(void) {
    for (unsigned int y = 0; y < field.height; ++y) {
      for (unsigned int x = 0; x < field.width; ++x) {
        draw_cell(x, y);
      }
    }
  }

  virtual char cell_available(void) { return cell_available_by_default; }

  virtual void calculate_cell_left_top(const unsigned int ax,
                                       const unsigned int ay,
                                       unsigned int &left,
                                       unsigned int &top) = 0;

  // Draw a hexagonal cell in the rectangular console buffer, the cell size is
  // configured in the "cell_scale" variable
  void draw_cell(const unsigned int ax, const unsigned int ay) {
    unsigned int x = 0;
    unsigned int y = 0;

    calculate_cell_left_top(ax, ay, x, y);

    x += cell_halfheight;

    // draw hexagonal cell's edge #1 (top)
    for (unsigned int i = 0; i < edge_width; ++i) {
      console[x][y] = '_';
      x++;
    }
    y++;

    // draw hexagonal cell's edge #2 (right upper)
    for (unsigned int i = 0; i < cell_halfheight; ++i) {
      console[x][y] = '\\';
      x++;
      y++;
    }

    // draw field contents (e.g. a space or a player's move)
    {
      unsigned int vx = x - edge_width;
      unsigned int vy = y - 1;
      console[vx][vy] = cellvalue2char(field.cells[ax][ay]);
    }
    x--;

    // draw hexagonal cell's edge #3 (right lower)
    for (unsigned int i = 0; i < cell_halfheight; ++i) {
      console[x][y] = '/';
      x--;
      y++;
    }
    y--;

    // draw hexagonal cell's edge #4 (bottom)
    for (unsigned int i = 0; i < edge_width; ++i) {
      console[x][y] = '_';
      x--;
    }

    // draw hexagonal cell's edge #5 (lower left)
    for (unsigned int i = 0; i < cell_halfheight; ++i) {
      console[x][y] = '\\';
      x--;
      y--;
    }
    x++;

    // draw hexagonal cell's edge #6 (upper left)
    for (unsigned int i = 0; i < cell_halfheight; ++i) {
      console[x][y] = '/';
      x++;
      y--;
    }
  }

public:
  // default constructor
  ConsoleConventionalHexagons(const Field &a_field,
                              const unsigned int a_cell_scale)
      : ConsoleAbstract(a_field) {
    const unsigned int ascii_hexagon_scale_factor = 4;
    unsigned int half_size = a_cell_scale * ascii_hexagon_scale_factor;
    cell_halfwidth = (half_size - cell_border_overlap);
    cell_height = (half_size >> 1);
    cell_halfheight = (cell_height >> 1);
    edge_width = cell_halfwidth - cell_halfheight;
  };
};

class ConsoleRectangularBoard : public ConsoleConventionalHexagons {
protected:
  // calculate initial coordinates of a cell's left upper corner
  virtual void calculate_cell_left_top(const unsigned int ax,
                                       const unsigned int ay,
                                       unsigned int &left, unsigned int &top) {
    left = ax * cell_halfwidth;
    top = ay * cell_height;

    if (is_odd_uint(ax)) {
      top += cell_halfheight;
    }
  }

public:
  ConsoleRectangularBoard(const Field &a_field, const unsigned int a_cell_scale)
      : ConsoleConventionalHexagons(a_field, a_cell_scale) {
    console_width =
        cell_halfwidth * field.width + cell_border_overlap * a_cell_scale;
    console_height = cell_height * field.height + cell_halfheight + 1;
  };
  virtual BoardIncline get_board_incline(void) { return BoardIncline::rect; }
};

class ConsoleParallelogramBoard : public ConsoleConventionalHexagons {
protected:
  // calculate initial coordinates of a cell's left upper corner
  virtual void calculate_cell_left_top(const unsigned int ax,
                                       const unsigned int ay,
                                       unsigned int &left, unsigned int &top) {
    left = ax * cell_halfwidth;
    top = ay * cell_height + ((field.width - ax - 1)) * cell_halfheight;
  }

public:
  virtual BoardIncline get_board_incline(void) { return BoardIncline::right; }
  ConsoleParallelogramBoard(const Field &a_field,
                            const unsigned int a_cell_scale)
      : ConsoleConventionalHexagons(a_field, a_cell_scale) {
    console_width =
        cell_halfwidth * field.width + cell_border_overlap * a_cell_scale;
    console_height = (cell_height * field.height) +
                     (cell_halfheight * (field.height - 1)) + 1;
  };
};

/*********************************************************************

The Graph classes: the abstract one and one that is familiar with the
HEX cells can can build connected nodes from the HEX board of moves

**********************************************************************/

// An enum type that holds the result of the calculation of the shortest path
// algorithm. If the result is "No", than the graph is not connected and there
// is no path between the "start" and "finish" nodes
enum class PathExists : uint8_t { DidNotCalculate, No, Yes };

// The ShortestPath class holds the results of the shortest path calculation,
// and it is separate from the Graph
template <typename Node, typename Distance> class ShortestPath {
public:
  using NodeSet = std::set<Node>;
  using NodeList = std::list<Node>;
  using WeightedNodeMap = std::map<Distance, NodeSet>;
  using NodeWeights = std::map<Node, Distance>;
  PathExists path_exists;
  NodeList path;

  // here we use initialize constructor to fill the class members with default
  // values
  ShortestPath() : path_exists(PathExists::DidNotCalculate) {};
};

// the main class - the "Graph" - implemented as a set of edges ("the edge
// lists") rather than an array of nodes ("connectivity matrices")
template <typename Node, typename Distance> class BasicGraph {
public:
  using Edge = std::pair<Node, Node>;

private:
  using NodeSet = std::set<Node>;
  using EdgeLayoutMap = std::map<Node, NodeSet>;
  using EdgeLengthMap = std::map<Edge, Distance>;

  // We have a map of sets to keep the edges - keys of the map are "a" and the
  // values of the set are "b"
  EdgeLayoutMap Edges;

  // Here is the map of edge lengths: key is the edge in form "std::pair<a,b>"
  // (b is never lower than a) and the value is the length per se of that node
  EdgeLengthMap EdgeLengths;

  inline Edge make_edge(const Node &a, const Node &b) {
    if (b < a) // make the correct order of the nodes in the edge-pair, so the
               // node with lower value will never be the second
      return std::make_pair(b, a);
    else
      return std::make_pair(a, b);
  }

public:
  bool node_exists(const Node &a) {
    const typename EdgeLayoutMap::const_iterator map_iter = Edges.find(a);
    return (map_iter != Edges.end());
  }

  // tests whether the specified edge exists in the graph
  bool adjacent(const Node &a, const Node &b) {
    const typename EdgeLayoutMap::const_iterator map_iter = Edges.find(a);
    if (map_iter == Edges.end())
      return false;
    const typename NodeSet::const_iterator set_iter = map_iter->second.find(b);
    if (set_iter == map_iter->second.end())
      return false;
    return true;
  };

  // lists all nodes "b" such that there is an edge from "a" to "b"
  NodeSet &neighbors(const Node &a) { return Edges[a]; }

  // adds the edge from a to b, if it is not there.
  void add_edge(const Node &a, const Node &b) {
    // Since we have undirected (bidirectional) edges, add for both directions,
    // to be able to find the edge faster in future - we have a map of sets
    Edges[a].insert(
        b); // Insert does nothing if a value already exists in a set
    Edges[b].insert(a);
  }

  // removes the edge from a to b, if it is there
  void delete_edge(const Node &a, const Node &b) {
    typename EdgeLayoutMap::iterator map_iter;

    // delete the edge from a to b
    map_iter = Edges.find(a);
    if (map_iter != Edges.end()) {
      map_iter->second.erase(b);
      if (map_iter->second.empty())
        Edges.erase(map_iter); // erase the empty set from the map by the
                               // iterator which points to the empty set
    }

    // delete the edge from b to a
    map_iter = Edges.find(b);
    if (map_iter != Edges.end()) {
      map_iter->second.erase(a);
      if (map_iter->second.empty())
        Edges.erase(map_iter);
    }
  }

  // set the distance value of an existing edge from a to b
  void set_edge_length(const Node &a, const Node &b, const Distance &v) {
    EdgeLengths[make_edge(a, b)] = v;
  }

  // get the distance value of an existing edge from a to b
  Distance get_edge_length(const Node &a, const Node &b) {
    return EdgeLengths[make_edge(a, b)];
  }

  // add the edge from a to b and set the distance value of the newly created
  // edge
  void add_edge_with_length(const Node &a, const Node &b, const Distance &v) {
    add_edge(a, b);
    set_edge_length(a, b, v);
  }

  // print all the edges in a graph, in human-readable form
  void print_edges(void) {
    typename EdgeLayoutMap::const_iterator map_iterator;
    for (map_iterator = Edges.begin(); map_iterator != Edges.end();
         ++map_iterator) {
      typename NodeSet::const_iterator set_iterator;
      for (set_iterator = map_iterator->second.begin();
           set_iterator != map_iterator->second.end(); ++set_iterator) {
        if (map_iterator->first < *set_iterator) {
          std::cout << map_iterator->first << " <-> " << (*set_iterator)
                    << " ( "
                    << get_edge_length(map_iterator->first, *set_iterator)
                    << " )" << std::endl;
        }
      }
    }
  }

  /*

  The "calculate_shortest_path" method calculates the shortest path and fills
  the "p" class with the resulting track, if any.

  This is the implementation of the Dijkstra's algorithm called the
  "Uniform-cost search (UCS)". To use the algorithm, fill the graph(an instance
  of the "Graph" class) with weighted edges and call the
  "calculate_shortest_path()" method of the graph. It will fill the instance of
  the "ShortestPath" with the results.

  Note: The shortest path algorithms, for convenience, are all, collectively,
  called the "Dijkstra's algorithm", although the implementation can differ
  significantly. In common presentations of Dijkstra's algorithm, initially, all
  nodes are entered into the priority queue. This is, however, not necessary.
  This is the implementation of the algorithm, as explained by Ira Pohl, called
  the "Uniform-cost search (UCS)". The UCS implementation starts with a priority
  queue ("open set" or "frontier") that contains only one item - the "start"
  node, and inserts new items as they are discovered. This variant has the same
  worst-case bounds as the common variant, but maintains a smaller priority
  queue in practice, speeding up the queue operations. In common presentations
  of Dijkstra's algorithm, initially, all nodes are entered into the priority
  queue. This is, however, not necessary: the algorithm can start with a
  priority queue that contains only one item, and insert new items as they are
  discovered.

  */

  void calculate_shortest_path_on_general_graph(
      const Node &start_node, const Node &finish_node,
      const Distance &initial_weight, ShortestPath<Node, Distance> &p) {
    using Map = typename ShortestPath<Node, Distance>::WeightedNodeMap;
    using Set = typename ShortestPath<Node, Distance>::NodeSet;
    using Weights = typename ShortestPath<Node, Distance>::NodeWeights;
    using SetSizeType = typename Set::size_type;

    Weights node_weights;
    SetSizeType total_nodes_walked;

    // first step of the algorithm - walk from the "start" node throughout the
    // whole graph and assign weights to all the reachable nodes
    {
      Map open_set;   // the moving frontier, i.e. the priority queue, initially
                      // containing the "start" node only
      Set closed_set; // set of explored nodes

      // add first node to the "open set", assigning the weight of the first
      // node to the initial_weight value (most probably, zero)
      open_set[initial_weight].insert(start_node);
      node_weights[start_node] = initial_weight;

      bool found_finish_node = false;
      Distance finish_node_weight;

      while (!open_set.empty()) {
        const typename Map::iterator &iter =
            open_set
                .begin(); // get the pointer to the first element in the "open
                          // set"; since the set is sorted by distance, we get
                          // the node with the lowest distance

        // get the node and its weight, from the "open set"
        const Distance current_weight =
            iter->first; // copy the node weight by value since the node will be
                         // deleted from the "open set"
        Set &nodeset = iter->second;
        typename Set::iterator node_iterator = nodeset.begin();
        const Node n =
            *node_iterator; // copy the node id by value since the node will be
                            // deleted from the "open set"

        // delete the node from the "open set"
        nodeset.erase(node_iterator);

        // if there are no more nodes with the same priority, then remove the
        // entry for the priority
        if (nodeset.empty())
          open_set.erase(iter);

        // add the node to the closed set
        closed_set.insert(n);

        // add node neighbours to the "open set"
        const NodeSet &neighbours = neighbors(n);
        for (typename NodeSet::const_iterator neighbors_iter =
                 neighbours.begin();
             neighbors_iter != neighbours.end(); ++neighbors_iter) {
          const Node &neighbour = *neighbors_iter;

          const Distance new_weight =
              current_weight + get_edge_length(n, neighbour);

          // update node weight
          typename Weights::iterator node_weight_iter =
              node_weights.find(neighbour);
          if (node_weight_iter == node_weights.end()) {
            node_weights[neighbour] = new_weight;
          } else {
            if (new_weight < node_weight_iter->second) {
              node_weight_iter->second = new_weight;
            }
          }

          if ((finish_node == neighbour) &&
              ((!found_finish_node) || (new_weight < finish_node_weight))) {
            found_finish_node = true;
            finish_node_weight = new_weight;
          }

          // not inserting all nodes allows us to find the shortest path from a
          // single source to the closest of a set of target nodes on a very
          // large graph without iterating through all the nodes - we can even
          // work with infinite graphs this way
          if (found_finish_node && finish_node_weight < new_weight) {
            // do not check the nodes that are already weighed more than the
            // finish node
            continue;
          }

          if (closed_set.find(neighbour) == closed_set.end()) {
            // add node to the open set if it is not in the closed set already,
            // i.e. if the node is not yet explored, add it to the frontier
            open_set[new_weight].insert(neighbour);
          }
        }
      }
      if (!found_finish_node) {
        p.path_exists = PathExists::No;
        return;
      }

      // we have walked through all reachable nodes
      total_nodes_walked = closed_set.size();
    }

    // second step of the algorithm: walk from the "finish" node back to the
    // "start" node, building the shortest path along the way
    {
      Node current_node = finish_node;
      for (SetSizeType i = 0; i < total_nodes_walked; ++i) {
        p.path.push_front(current_node);
        if (current_node == start_node) {
          p.path_exists = PathExists::Yes;
          return;
        }
        const NodeSet &neighbours = neighbors(current_node);
        Node least_node;
        Distance least_weight;
        bool found_neighbour = false;
        const Distance current_weight = node_weights[current_node];
        for (typename NodeSet::const_iterator iter = neighbours.begin();
             iter != neighbours.end(); ++iter) {
          const Distance this_weight = node_weights[*iter];
          const Distance edge_value = get_edge_length(*iter, current_node);
          if (this_weight + edge_value != current_weight)
            continue;

          if ((!found_neighbour) || (this_weight < least_weight)) {
            found_neighbour = true;
            least_weight = this_weight;
            least_node = *iter;
          }
        }
        if (!found_neighbour) {
          // we didn't find a new neighbour - path does not exist
          break;
        }
        current_node = least_node;
      }
      p.path_exists = PathExists::No;
    }
  }

  // the helper function that returns just the boolean result if we don't need
  // the full path to follow
  bool path_exists(const Node &start_node, const Node &finish_node) {
    const Distance initial_weight = 1;
    typedef ShortestPath<Node, Distance> PathResult;
    PathResult *p = new PathResult;
    calculate_shortest_path_on_general_graph(start_node, finish_node,
                                             initial_weight, *p);
    bool path_exists = (p->path_exists == PathExists::Yes);
    delete p;
    p = nullptr;
    return path_exists;
  }
};

// the main class - the "Graph" - implemented as a set of edges ("the edge
// lists") rather than an array of nodes ("connectivity matrices")
template <typename Node, typename Distance>
class FieldGraph : public BasicGraph<Node, Distance> {
public:
  const Field &field;
  const BoardIncline board_incline;

  int coord2nodeid(const unsigned int x, const unsigned int y) const {
    return (y * field.width) + x;
  }

  void connect_node_if_marked(const Node &a, const CellValue mark,
                              const unsigned int x, const unsigned int y) {
    if (!field.coord_out_of_range_xy(x, y)) {
      if (field.cells[x][y] == mark) {
        const Node b = coord2nodeid(x, y);
        const Distance default_distance = 1;
        BasicGraph<Node, Distance>::add_edge_with_length(a, b,
                                                         default_distance);
      }
    }
  }

  // initialize the graph from the 2-dimensional field
  FieldGraph(const Field &a_field, const CellValue mark,
             const BoardIncline a_board_incline)
      : field(a_field), board_incline(a_board_incline) {
    for (unsigned int y = 0; y < field.height; ++y) {
      for (unsigned int x = 0; x < field.width; ++x) {
        if (field.cells[x][y] == mark) {
          const Node a = coord2nodeid(x, y);
          if (y > 0)
            connect_node_if_marked(a, mark, x, y - 1);
          connect_node_if_marked(a, mark, x, y + 1);
          if (x > 0)
            connect_node_if_marked(a, mark, x - 1, y);
          connect_node_if_marked(a, mark, x + 1, y);
          switch (board_incline) {
          case BoardIncline::rect:
            // use non-conventionally-calculated coordinates for the 5th and 6th
            // nodes
            {
              unsigned int y2 = 0;
              bool have_y2 = false;
              if (is_odd_uint(x)) {
                if (y + 1 < field.height) {
                  y2 = y + 1;
                  have_y2 = true;
                }
              } else {
                if (y > 0) {
                  y2 = y - 1;
                  have_y2 = true;
                }
              }
              if (have_y2) {
                if (x > 0)
                  connect_node_if_marked(a, mark, x - 1, y2);
                connect_node_if_marked(a, mark, x + 1, y2);
              }
            }
            break;

          case BoardIncline::right:
            if (x > 0 && y > 0)
              connect_node_if_marked(a, mark, x - 1, y - 1);
            connect_node_if_marked(a, mark, x + 1, y + 1);
            break;

          default: // case BoardInlcline::left:
            if (x > 0)
              connect_node_if_marked(a, mark, x - 1, y + 1);
            if (y > 0)
              connect_node_if_marked(a, mark, x + 1, y - 1);
            break;
          }
        }
      }
    }
  }

  void add_terminal_node_if_path_exists(const Node &a, const unsigned int x,
                                        const unsigned int y) {
    const Node n = coord2nodeid(x, y);
    if (BasicGraph<Node, Distance>::node_exists(n)) {
      const Distance default_distance = 1;
      BasicGraph<Node, Distance>::add_edge_with_length(a, n, default_distance);
    }
  }

  void add_start_stop_edge_horizontal(const Node &start, const Node &stop) {
    const unsigned int ystart = 0;
    const unsigned int ystop = field.height - 1;
    for (unsigned int x = 0; x < field.width; ++x) {
      add_terminal_node_if_path_exists(start, x, ystart);
      add_terminal_node_if_path_exists(stop, x, ystop);
    }
  }

  void add_start_stop_edge_vertical(const Node &start, const Node &stop) {
    const unsigned int xstart = 0;
    const unsigned int xstop = field.width - 1;
    for (unsigned int y = 0; y < field.height; ++y) {
      add_terminal_node_if_path_exists(start, xstart, y);
      add_terminal_node_if_path_exists(stop, xstop, y);
    }
  }
};

// Check whether a winning path exists according to the rules of the HEX game,
// using a conventional graph and Dijkstra's algorithm
bool is_game_won(const Field &field, const CellValue mark,
                 const bool start_stop_edges_vertical,
                 const BoardIncline board_incline) {
  using MyNodeType = uint8_t;
  using MyDistanceType = uint8_t;
  using MyGraph = FieldGraph<MyNodeType, MyDistanceType>;

  // assign the IDs for the virtual nodes
  // we will use virtual, "start" and "stop" nodes connected to all side lines of
  // the board, either horizontal or vertical

  // create a graph from the field of all the moves, marked by a "mark"
  // character, made by the appropriate single player
  MyGraph g(field, mark, board_incline);
  const MyNodeType start_id =
      g.coord2nodeid(field.width - 1, field.height - 1) + 1;
  const MyNodeType stop_id = start_id + 1;

  if (start_stop_edges_vertical)
    g.add_start_stop_edge_vertical(start_id, stop_id);
  else
    g.add_start_stop_edge_horizontal(start_id, stop_id);

  return g.path_exists(start_id, stop_id);
}

/******************************************************************************
The classes that implement the computer's moves.
(1) The basic abstract class
(2) The ComputerMoveRandom that just uses a random generator to make an
arbitrary move (3) The ComputerMoveAI that uses Monte-Carlo simulation to make a
better move
*******************************************************************************/

// the row and column are 1-based with the move coordinates
const unsigned int first_row_base = 1;
const unsigned int first_column_base = 1;

class ComputerMoveAbstract {
public:
  virtual void make_move(const CellValue mark, const BoardIncline board_incline,
                         unsigned int &row, unsigned int &column) = 0;
  virtual ~ComputerMoveAbstract() {};
};

class ComputerMoveRandom : public ComputerMoveAbstract {
private:
  std::default_random_engine random_engine;
  std::random_device random_device;

  using MyDistribution = std::uniform_int_distribution<int>;
  MyDistribution *random_distribution_row, *random_distribution_column;

  Field &field;

public:
  ComputerMoveRandom(Field &a_field) : field(a_field) {
    // create a random number generator using C++ STL
    random_engine.seed(random_device());
    auto board_width = field.width;
    auto board_height = field.height;
    random_distribution_row = new MyDistribution(first_row_base, board_height);
    random_distribution_column =
        new MyDistribution(first_column_base, board_width);
  }

  ~ComputerMoveRandom() {
    delete random_distribution_row;
    random_distribution_row = nullptr;
    delete random_distribution_column;
    random_distribution_column = nullptr;
  }

  virtual void make_move(const CellValue player_mark,
                         const BoardIncline /* board_incline */,
                         unsigned int &row, unsigned int &column) {
    for (unsigned int x, y, r, c;;) {
      do {
        r = (*random_distribution_row)(random_engine);
        y = r - first_row_base;
      } while (field.coord_out_of_range_y(y));

      do {
        c = (*random_distribution_column)(random_engine);
        x = c - first_column_base;
      } while (field.coord_out_of_range_x(x));

      CellValue &cv = field.cells[x][y]; // initialize the reference to further
                                         // modify the field value
      if (cv == CellValue::Blank) {
        cv = player_mark;
        row = r;
        column = c;
        break;
      } else {
        continue;
      }
    }
  }
};

/****************************************************************************
The "XorShift" Pseudorandom Number Generator

For more information, see https://en.wikipedia.org/wiki/Xorshift or the paper
"Xorshift RNGs" by George Marsaglia from The Florida State University,
published in 2003 in the Journal of Statistical Software,
DOI: 10.18637/jss.v008.i14
or just search for "Xorshift RNGs" on Semantic Scholar
*****************************************************************************/

const uint64_t initial_seed =
    88172645463325252LL; // the intial seed value from the paper above mentioned

uint64_t current_seed = initial_seed;

inline uint64_t xorshift64(uint64_t &seed) {
  seed ^= (seed << 13);
  seed ^= (seed >> 7);
  return (seed ^= (seed << 17));
}

// use one xorshift call to return one pseudorandom bytes each in the given
// range from 0 but less to the given limit, i.e. [0..limit)

inline uint8_t rand_byte_lim(uint64_t &seed, const uint8_t limit) {
  uint64_t x = xorshift64(seed);
  x &= 0xffffffffffffff; // 7 bytes
  x *= limit;
  x >>= 56; // 7 bytes * 8 bits = 56 bits
  return x & 0xff;
}

// use one xorshift call to return 8 pseudorandom bytes each in own range
// [0..limN)
inline void rand_8_bytes(uint64_t &seed, const uint8_t lim1, uint8_t &res1,
                         const uint8_t lim2, uint8_t &res2, const uint8_t lim3,
                         uint8_t &res3, const uint8_t lim4, uint8_t &res4,
                         const uint8_t lim5, uint8_t &res5, const uint8_t lim6,
                         uint8_t &res6, const uint8_t lim7, uint8_t &res7,
                         const uint8_t lim8, uint8_t &res8) {
  // unrolled for CPU instruction-level parallelism to allow out-of-order
  // execution for better performance

  uint64_t x = xorshift64(seed);
  uint32_t t1 = (x >> (0 * 8)) & 0xff;
  uint32_t t2 = (x >> (1 * 8)) & 0xff;
  uint32_t t3 = (x >> (2 * 8)) & 0xff;
  uint32_t t4 = (x >> (3 * 8)) & 0xff;
  uint32_t t5 = (x >> (4 * 8)) & 0xff;
  uint32_t t6 = (x >> (5 * 8)) & 0xff;
  uint32_t t7 = (x >> (6 * 8)) & 0xff;
  uint32_t t8 = (x >> (7 * 8)) & 0xff;
  t1 *= lim1;
  t2 *= lim2;
  t3 *= lim3;
  t4 *= lim4;
  t5 *= lim5;
  t6 *= lim6;
  t7 *= lim7;
  t8 *= lim8;
  res1 = (t1 >> 8) & 0xff;
  res2 = (t2 >> 8) & 0xff;
  res3 = (t3 >> 8) & 0xff;
  res4 = (t4 >> 8) & 0xff;
  res5 = (t5 >> 8) & 0xff;
  res6 = (t6 >> 8) & 0xff;
  res7 = (t7 >> 8) & 0xff;
  res8 = (t8 >> 8) & 0xff;
}

/*********************************************************************************************************************

START OF ARTIFICIAL INTELLIGENCE CODE / THE MONTE-CARLO SIMULATION

*********************************************************************************************************************/

using NodeCoord = uint8_t; // the node coordinate on one axe, i.e. from 0 to 10
                           // on an 11x11 board
using NodeIndex = uint8_t; // the node index from 0 to the total number of cells
                           // on the board, i.e. from 0 to 120 on an 11x11 board
using AiMonteCarloAttemptCount =
    unsigned int; // an integer type used to keep the maximum number of
                  // Monte-Carlo attempts that we are doing for each turn

/*
The class to implement the Monte-Carlo simulation. It does not use graphs or
STL in order to be very fast. It has fixed-size arrays and is thus built via a
template for each of the board sizes, so the board size is compiled to the
machine code as a constant rather than as a variable, thus achieving faster
calculations for Monte-Carlo.
*/

template <const NodeCoord AiFieldSideLength>
class ComputerMoveAI : public ComputerMoveAbstract {
public:
  static const NodeIndex AiFieldSize =
      AiFieldSideLength *
      AiFieldSideLength; // the total number of cells on the board
  static const NodeIndex AiFieldSizeInQWords =
      (AiFieldSize + 63) >> 6; // used to caclulate the number of bits to
                               // allocate one bit for every cell
  static const NodeIndex AiBlueFinishVerticalRightNodeIndex = AiFieldSize;

  // calculated once for a given board size, e.g. 11x11
  // this table is needed in order for each of the nodes has indexes of all its
  // adjacent nodes for faster path calculation
  uint64_t AiPermanentBoardCellsMeshForBlue[AiFieldSize];

  /*
  To check whether the path exists for the Monte-Carlo simulation, we do not
  necessarily need the "shortest" one, so we use a simplified algorithm that
  does not calculate the distance. Thus we only keep track of the tentative
  nodes ("open set") and don't track visited nodes ("closed set") separately.
  */

  using AiTentativeNodeIndexesArrayType = NodeIndex;
  static const AiTentativeNodeIndexesArrayType AiTentativeNodeIndexesArraySize =
      AiFieldSize; // each node may add 6 more nodes
  NodeIndex AiTentativeNodeIndexes[AiTentativeNodeIndexesArraySize];

  CellValue AiCurrentTurnFieldCellValuesSaved[AiFieldSize];
  CellValue AiFilled[AiFieldSize];

  AiMonteCarloAttemptCount AiProposedMovesPlayer[AiFieldSize];
  AiMonteCarloAttemptCount AiProposedMovesOpponent[AiFieldSize];

  NodeIndex EmptyCellsSavedIndices[AiFieldSize];
  NodeIndex EmptyCellsSavedCount;

  NodeIndex RedCellsSavedIndices[AiFieldSize];
  NodeIndex RedCellsSavedCount;

  NodeIndex BlueCellsSavedIndices[AiFieldSize];
  NodeIndex BlueCellsSavedCount;

  NodeIndex EmptyFieldsShuffled[AiFieldSize];

  void ShuffleEmptyCells(void) {
    uint64_t seed =
        current_seed; // load the xorshif seed from memory to the local register
                      // variable for faster processing

    NodeIndex i = 0;

    /* Shuffle 8 cells in one block to combine the call to the pseudorandom
    function for CPU instruction-level parallelism to allow out-of-order
    execution for better performance. */

    while (i + 8 < EmptyCellsSavedCount) {
      NodeIndex ja, jb, jc, jd, je, jf, jg, jh;

      // Call one xorshift to get 8 pseudorandom values in the given range
      rand_8_bytes(seed, i + 1, ja, i + 2, jb, i + 3, jc, i + 4, jd, i + 5, je,
                   i + 6, jf, i + 7, jg, i + 8, jh);

      EmptyFieldsShuffled[i] = EmptyFieldsShuffled[ja];
      EmptyFieldsShuffled[ja] = EmptyCellsSavedIndices[i];

      EmptyFieldsShuffled[i + 1] = EmptyFieldsShuffled[jb];
      EmptyFieldsShuffled[jb] = EmptyCellsSavedIndices[i + 1];

      EmptyFieldsShuffled[i + 2] = EmptyFieldsShuffled[jc];
      EmptyFieldsShuffled[jc] = EmptyCellsSavedIndices[i + 2];

      EmptyFieldsShuffled[i + 3] = EmptyFieldsShuffled[jd];
      EmptyFieldsShuffled[jd] = EmptyCellsSavedIndices[i + 3];

      EmptyFieldsShuffled[i + 4] = EmptyFieldsShuffled[je];
      EmptyFieldsShuffled[je] = EmptyCellsSavedIndices[i + 4];

      EmptyFieldsShuffled[i + 5] = EmptyFieldsShuffled[jf];
      EmptyFieldsShuffled[jf] = EmptyCellsSavedIndices[i + 5];

      EmptyFieldsShuffled[i + 6] = EmptyFieldsShuffled[jg];
      EmptyFieldsShuffled[jg] = EmptyCellsSavedIndices[i + 6];

      EmptyFieldsShuffled[i + 7] = EmptyFieldsShuffled[jh];
      EmptyFieldsShuffled[jh] = EmptyCellsSavedIndices[i + 7];

      i += 8;
    }

    // shuffle remaining cells which are less than 8
    while (i < EmptyCellsSavedCount) {
      uint8_t j = rand_byte_lim(seed, i + 1);
      EmptyFieldsShuffled[i] = EmptyFieldsShuffled[j];
      EmptyFieldsShuffled[j] = EmptyCellsSavedIndices[i];
      i++;
    }

    current_seed =
        seed; // save back current seed from the local register back to memory
  }

  void FillBoardNonEmpty(void) {
    for (NodeIndex i = 0; i < BlueCellsSavedCount; ++i) {
      AiFilled[BlueCellsSavedIndices[i]] = CellValue::Blue;
    }

    for (NodeIndex i = 0; i < RedCellsSavedCount; ++i) {
      AiFilled[RedCellsSavedIndices[i]] = CellValue::Red;
    }
  }

  inline NodeIndex FillBoardEmptyOdd(const CellValue player_mark) {
    NodeIndex i = 0;
    NodeIndex half = EmptyCellsSavedCount >> 1;
    for (NodeIndex k = 0; k < half; ++k) {
      AiFilled[EmptyFieldsShuffled[i++]] = CellValue::Blue;
      AiFilled[EmptyFieldsShuffled[i++]] = CellValue::Red;
    }
    NodeIndex proposed_move = EmptyFieldsShuffled[i];
    AiFilled[proposed_move] = player_mark;
    return proposed_move;
  }

  inline NodeIndex FillBoardEmptyEven(const CellValue player_mark) {
    NodeIndex ec0 = EmptyFieldsShuffled[0];
    NodeIndex ec1 = EmptyFieldsShuffled[1];
    NodeIndex proposed_move;

    // set node idex of the player mark
    if (player_mark == CellValue::Blue)
      proposed_move = ec0;
    else
      proposed_move = ec1;

    // fill first 2 nodes
    AiFilled[ec0] = CellValue::Blue;
    AiFilled[ec1] = CellValue::Red;

    // fill remaining nodes is more than 2
    NodeIndex i = 2;
    NodeIndex half = EmptyCellsSavedCount >> 1;
    for (NodeIndex k = 1; k < half; ++k) {
      AiFilled[EmptyFieldsShuffled[i++]] = CellValue::Blue;
      AiFilled[EmptyFieldsShuffled[i++]] = CellValue::Red;
    }

    return proposed_move;
  }

  inline void AiAddNeighbour(uint64_t &pack, const NodeCoord x,
                             const NodeCoord y, uint8_t &bytes) {
    const NodeIndex idx = (y * AiFieldSideLength) + x;
    pack = ((pack << 8) | idx);
    bytes++;
  }

  inline bool
  is_odd_NodeIndex(const NodeIndex i) // just returns whether an node index is
                                      // odd, i.e. is 1, 3, etc...
  {
    return i & 1;
  }

  inline bool
  is_odd_NodeCoord(const NodeCoord c) // just returns whether an node coordinate
                                      // is odd, i.e. is 1, 3, etc...
  {
    return c & 1;
  }

  void FillAiFieldMeshForBlue(const BoardIncline board_incline) {
    for (NodeCoord y = 0; y < AiFieldSideLength; ++y) {
      for (NodeCoord x = 0; x < AiFieldSideLength; ++x) {
        uint64_t n = 0;
        uint8_t b = 0;
        const NodeIndex a = (y * AiFieldSideLength) + x;
        if (y > 0)
          AiAddNeighbour(n, x, y - 1, b);
        if (y + 1 < AiFieldSideLength)
          AiAddNeighbour(n, x, y + 1, b);
        if (x > 0)
          AiAddNeighbour(n, x - 1, y, b);
        if (x + 1 < AiFieldSideLength)
          AiAddNeighbour(n, x + 1, y, b);
        switch (board_incline) {
        case BoardIncline::rect:
          // use non-conventionally-calculated coordinates for the 5th and 6th
          // nodes
          if (is_odd_NodeCoord(x)) // is odd
          {
            if (y + 1 < AiFieldSideLength) {
              NodeIndex y2 = y + 1;
              if (x > 0)
                AiAddNeighbour(n, x - 1, y2, b);
              if (x + 1 < AiFieldSideLength)
                AiAddNeighbour(n, x + 1, y2, b);
            }
          } else {
            if (y > 0) {
              NodeIndex y2 = y - 1;
              if (x > 0)
                AiAddNeighbour(n, x - 1, y2, b);
              if (x + 1 < AiFieldSideLength)
                AiAddNeighbour(n, x + 1, y2, b);
            }
          }
          break;

        case BoardIncline::right:
          if ((x > 0) && (y > 0))
            AiAddNeighbour(n, x - 1, y - 1, b);
          if ((x + 1 < AiFieldSideLength) && (y + 1 < AiFieldSideLength))
            AiAddNeighbour(n, x + 1, y + 1, b);
          break;

        default: // case BoardIncline::left:
          if ((x > 0) && (y + 1 < AiFieldSideLength))
            AiAddNeighbour(n, x - 1, y + 1, b);
          if ((y > 0) && (x + 1 < AiFieldSideLength))
            AiAddNeighbour(n, x + 1, y - 1, b);
          break;
        }
        if (b < 6) {
          // special handling for edge nodes, having less than 6 neighbours
          NodeIndex idx;
          if (x + 1 == AiFieldSideLength) {
            // right vertical edge for the blue
            idx = AiBlueFinishVerticalRightNodeIndex;
          } else {
            idx = a; // the node itself
          }
          for (; b < 6; ++b) {
            n = ((n << 8) | idx);
          }
        }

        AiPermanentBoardCellsMeshForBlue[a] = n;
      }
    }
  }

  inline uint64_t bit_test_and_set_tentative(const uint8_t idx,
                                             uint64_t *AiTentativeNodesBits) {
    const uint64_t one = 1;
    if (AiFieldSizeInQWords > 1) {
      const uint8_t ofs = idx >> 6;
      const uint8_t bit = idx & 0x3f;
      uint64_t &addr = AiTentativeNodesBits[ofs];
      uint64_t mask = one << bit;
      uint64_t ret = addr & mask;
      addr |= mask;
      return ret;
    } else {
      // for the boards of less than 64 cells (e.g. 7x7) we just need a single
      // 64-bit variable to keep all the bits of the tentative nodes, one bit
      // per node
      uint64_t mask = one << idx;
      uint64_t ret = AiTentativeNodesBits[0] & mask;
      AiTentativeNodesBits[0] |= mask;
      return ret;
    }
  }

  inline void
  AddTentativeNode(const NodeIndex a_node,
                   AiTentativeNodeIndexesArrayType &AiTentativeNodeEndIndex,
                   uint64_t *AiTentativeNodesBits) {
    if (!bit_test_and_set_tentative(a_node, AiTentativeNodesBits)) {
      // the node has not yet been visited
      AiTentativeNodeIndexes[AiTentativeNodeEndIndex++] = a_node;
    }
  }

  inline void AddTentativeNodeIfBlue(
      const NodeIndex a_node,
      AiTentativeNodeIndexesArrayType &AiTentativeNodeEndIndex,
      uint64_t *AiTentativeNodesBits) {
    if (AiFilled[a_node] == CellValue::Blue) {
      AddTentativeNode(a_node, AiTentativeNodeEndIndex, AiTentativeNodesBits);
    }
  }

  /*
  The very fast routine to check whether the blue player has a connection (e.g.,
  the winning position). It is only called from the Monte-Carlo simulation when
  there are no empty fields, i.e., the whole board is filled randomly - what was
  empty space is shuffled randomly with blue and red marks. Since there is no
  "draw" in the hex game, once the whole board is filled and there is no
  connection for the "blue" player, this means that there is the connection for
  the "red" player, thus lack of connection for "blue" indicates the connection
  for "red".
  */

  bool AiPathExistsBlue(void) {

    uint64_t ba[AiFieldSizeInQWords]; // the bitmap to track all

    for (uint8_t i = 0; i < AiFieldSizeInQWords; ++i) {
      ba[i] = 0;
    }

    AiTentativeNodeIndexesArrayType ub;
    AiTentativeNodeIndexesArrayType ue;

    ub = 0;
    ue = 0;

    // Add leftmost vertical column to the set of unvisited nodes,
    // since the blue player moves from right to left (or vice versa),
    // i.e. needs to make a horizontal connection.

    for (NodeIndex y = 0; y < AiFieldSideLength; ++y) {
      const NodeIndex x = 0;
      const NodeIndex EdgeNode = (y * AiFieldSideLength) + x;
      AddTentativeNodeIfBlue(EdgeNode, ue, ba);
    }

    while (ub < ue) {
      const NodeIndex current_node = AiTentativeNodeIndexes[ub++];

      // process 6 nodes independently, for CPU instruction-level parallelism
      // to allow out-of-order execution for better performance

      const uint64_t neighbours =
          AiPermanentBoardCellsMeshForBlue[current_node];
      const NodeIndex n1 = (neighbours >> (8 * 0)) & 0xff;
      const NodeIndex n2 = (neighbours >> (8 * 1)) & 0xff;
      const NodeIndex n3 = (neighbours >> (8 * 2)) & 0xff;
      const NodeIndex n4 = (neighbours >> (8 * 3)) & 0xff;
      const NodeIndex n5 = (neighbours >> (8 * 4)) & 0xff;
      const NodeIndex n6 = (neighbours >> (8 * 5)) & 0xff;
      if ((n1 == AiBlueFinishVerticalRightNodeIndex) ||
          (n2 == AiBlueFinishVerticalRightNodeIndex) ||
          (n3 == AiBlueFinishVerticalRightNodeIndex) ||
          (n4 == AiBlueFinishVerticalRightNodeIndex) ||
          (n5 == AiBlueFinishVerticalRightNodeIndex) ||
          (n6 == AiBlueFinishVerticalRightNodeIndex)) {
        return true; // path exists for blue
      }
      AddTentativeNodeIfBlue(n1, ue, ba);
      AddTentativeNodeIfBlue(n2, ue, ba);
      AddTentativeNodeIfBlue(n3, ue, ba);
      AddTentativeNodeIfBlue(n4, ue, ba);
      AddTentativeNodeIfBlue(n5, ue, ba);
      AddTentativeNodeIfBlue(n6, ue, ba);
      assert(ub < AiTentativeNodeIndexesArraySize);
      assert(ue < AiTentativeNodeIndexesArraySize);
    }

    return false; // path not found for "blue" - this means that there is the
                  // path for "red"
  }

  Field &field;
  bool MeshFilledForBlue;

  void SaveFieldSnapshotForFurtherAiCalculations(void) {
    EmptyCellsSavedCount = 0;
    RedCellsSavedCount = 0;
    BlueCellsSavedCount = 0;

    for (NodeCoord x = 0; x < AiFieldSideLength; ++x) {
      for (NodeCoord y = 0; y < AiFieldSideLength; ++y) {
        const NodeIndex i = (y * AiFieldSideLength) + x;
        const CellValue v = field.cells[x][y];
        AiCurrentTurnFieldCellValuesSaved[i] = v;
        switch (v) {
        case CellValue::Red:
          RedCellsSavedIndices[RedCellsSavedCount++] = i;
          break;
        case CellValue::Blue:
          BlueCellsSavedIndices[BlueCellsSavedCount++] = i;
          break;
        default: //  case CellValue::Blank:
          EmptyCellsSavedIndices[EmptyCellsSavedCount++] = i;
          break;
        }
      }
    }
    for (NodeIndex i = RedCellsSavedCount; i < AiFieldSize; ++i) {
      RedCellsSavedIndices[i] =
          AiFieldSize; // fill the rest of the array with invalid values
    }

    for (NodeIndex i = BlueCellsSavedCount; i < AiFieldSize; ++i) {
      BlueCellsSavedIndices[i] =
          AiFieldSize; // fill the rest of the array with invalid values
    }

    for (NodeIndex i = EmptyCellsSavedCount; i < AiFieldSize; ++i) {
      EmptyCellsSavedIndices[i] =
          AiFieldSize; // fill the rest of the array with invalid values
    }
  }

public:
  ComputerMoveAI(Field &a_field)
      : EmptyCellsSavedCount(0), RedCellsSavedCount(0), BlueCellsSavedCount(0),
        field(a_field), MeshFilledForBlue(false) {
    assert((AiFieldSideLength == field.width) &&
           (AiFieldSideLength == field.height));
  }

  ~ComputerMoveAI() {}

  virtual void make_move(const CellValue player_mark,
                         const BoardIncline board_incline, unsigned int &row,
                         unsigned int &column) {

    if (!MeshFilledForBlue) {
      MeshFilledForBlue = true;
      FillAiFieldMeshForBlue(board_incline);
    }
    SaveFieldSnapshotForFurtherAiCalculations();
    switch (player_mark) {
    case CellValue::Red:
      assert(RedCellsSavedCount <=
             BlueCellsSavedCount); // we can only generate a move for a Red
                                   // player if there is no more red cells
                                   // already than blue cells on the board
      assert(BlueCellsSavedCount - RedCellsSavedCount <=
             1); // the difference between the number of blue and red cell
                 // should be no more than one
      break;
    case CellValue::Blue:
      assert(BlueCellsSavedCount <=
             RedCellsSavedCount); // we can only generate a move for a Blue
                                  // player if there is no more blue cells
                                  // already than red cells on the board
      assert(RedCellsSavedCount - BlueCellsSavedCount <=
             1); // the difference between the number of blue and red cell
                 // should be no more than one
      break;
    default:
      assert(false); // The player_mark should be either Red or Blue
      break;
    }

    assert(AiFieldSize ==
           EmptyCellsSavedCount + RedCellsSavedCount + BlueCellsSavedCount);

    assert(EmptyCellsSavedCount >
           0); // there should be at least one empty cell to fill

    for (NodeIndex n = 0; n < AiFieldSize; ++n) {
      AiProposedMovesPlayer[n] = 0;
      AiProposedMovesOpponent[n] = 0;
    }

    /* The total number of iterations for all possible moves. When less possible
    moves are remaining closer to the end of the game, more iterations are made
    for each remaining move, so the total number of Monte-Carlo evaluations is
    always constant. */

    const AiMonteCarloAttemptCount monte_carlo_iterations = 1000000;

    for (AiMonteCarloAttemptCount iter = 0; iter < monte_carlo_iterations;
         ++iter) {
      ShuffleEmptyCells();

      FillBoardNonEmpty();

      NodeIndex proposed_move;
      if (is_odd_NodeIndex(EmptyCellsSavedCount)) {
        proposed_move = FillBoardEmptyOdd(player_mark);
      } else {
        proposed_move = FillBoardEmptyEven(player_mark);
      }

      bool path_exists_for_blue = AiPathExistsBlue();

      bool is_red_player = (player_mark == CellValue::Red);

      // we count winning moves for both red and blue players for debugging only
      // to check the sum later since only counting for one player would have
      // been enough
      if (is_red_player ^ path_exists_for_blue) {
        ++(AiProposedMovesPlayer[proposed_move]);
      } else {
        ++(AiProposedMovesOpponent[proposed_move]);
      }
    }

    AiMonteCarloAttemptCount total_proposed_moves_for_player = 0;
    AiMonteCarloAttemptCount total_proposed_moves_for_opponent = 0;
    for (NodeIndex n = 0; n < AiFieldSize; ++n) {
      total_proposed_moves_for_player += AiProposedMovesPlayer[n];
      total_proposed_moves_for_opponent += AiProposedMovesOpponent[n];
    }

    assert(total_proposed_moves_for_player +
               total_proposed_moves_for_opponent ==
           monte_carlo_iterations);

    NodeIndex BestMoveIdx = 0;
    AiMonteCarloAttemptCount BestMoveWinningCount = 0;
    for (NodeIndex n = 0; n < AiFieldSize; n++) {
      const AiMonteCarloAttemptCount &current_move_winning_count =
          AiProposedMovesPlayer[n];
      if (current_move_winning_count > BestMoveWinningCount) {
        BestMoveIdx = n;
        BestMoveWinningCount = current_move_winning_count;
      }
    }
    NodeCoord x = (BestMoveIdx % AiFieldSideLength);
    NodeCoord y = (BestMoveIdx / AiFieldSideLength);

    CellValue &cv = field.cells[x][y]; // initialize the reference to further
                                       // modify the field value
    assert(cv == CellValue::Blank);
    cv = player_mark;
    column = x + first_column_base;
    row = y + first_row_base;
  }
};

/* The "factory" to get the custom Monte-Carlo AI computer player class
depending on the board size custom classes have constant-size buffers and the
constants built-in into the machine code that work faster than variables */

ComputerMoveAbstract *get_computer_move_ai_factory(const NodeCoord a_len,
                                                   Field &field) {
  using Ai3 = ComputerMoveAI<3>;
  using Ai4 = ComputerMoveAI<4>;
  using Ai5 = ComputerMoveAI<5>;
  using Ai6 = ComputerMoveAI<6>;
  using Ai7 = ComputerMoveAI<7>;
  using Ai8 = ComputerMoveAI<8>;
  using Ai9 = ComputerMoveAI<9>;
  using Ai10 = ComputerMoveAI<10>;
  using Ai11 = ComputerMoveAI<11>;
  using Ai12 = ComputerMoveAI<12>;

  switch (a_len) {
  case 3:
    return new Ai3(field);
  case 4:
    return new Ai4(field);
  case 5:
    return new Ai5(field);
  case 6:
    return new Ai6(field);
  case 7:
    return new Ai7(field);
  case 8:
    return new Ai8(field);
  case 9:
    return new Ai9(field);
  case 10:
    return new Ai10(field);
  case 11:
    return new Ai11(field);
  case 12:
    return new Ai12(field);
  default:
    return nullptr;
  }
}

/**************************************

END OF THE ARTIFICIAL INTELLIGENCE CODE

***************************************/

// Default setup variables

const int input_board_type_parallelogram = 1;
const int input_board_type_rectangle = 2;
const int input_board_type_compact_grid = 3;
const int input_board_type_min = 1;
const int input_board_type_max = 3;

const int input_board_type_default = input_board_type_parallelogram;

void print_field(const Field &field, BoardIncline &board_incline,
                 const int cell_scale,
                 const int board_type) // draws a filed to the console buffer
                                       // and ouputs the buffer to STDOUT
{
  ConsoleAbstract *console = nullptr;
  switch (board_type) {
  case input_board_type_rectangle:
    console = new ConsoleRectangularBoard(field, cell_scale);
    break;
  case input_board_type_compact_grid:
    console = new ConsoleCondensed(field);
    break;
  default: // input_board_type_parallelogram
    console = new ConsoleParallelogramBoard(field, cell_scale);
    break;
  }
  board_incline = console->get_board_incline();
  console->clear();
  console->draw_field();
  console->print();
  delete console;
  console = nullptr;
}

void ignore_line(void) {
  std::cin.clear();
  std::string str;
  std::getline(std::cin, str);
  std::cout << "\"" << str << "\" is not a valid integer" << std::endl;
}

int read_integer(const int a_default) {
  int retval = a_default;
  std::string str;
  std::getline(std::cin, str);
  if (!str.empty())
    try {
      retval = std::stoi(str);
    } catch (std::exception const &e) {
      std::cerr << e.what() << std::endl;
      retval = a_default;
    }
  return retval;
}

const int min_board_size = 3;
const int max_board_size = 12;
const int default_board_size = 4;

using Turn = unsigned int;

inline bool is_odd_turn(const Turn t) // just returns whether the turn number is
                                      // odd, i.e. is 1, 3, etc...
{
  return t & 1;
}

void randomize_xor64(void) {
  std::default_random_engine random_engine;
  std::random_device random_device;
  random_engine.seed(random_device());
  std::uniform_int_distribution<unsigned int> distribution(1000, 2000);
  unsigned int repeat;
  repeat = distribution(random_engine);
  for (unsigned int i = 0; i < repeat; ++i)
    xorshift64(current_seed);
}

int main() {

  std::cout << "Welcome to HEX!" << std::endl;

  randomize_xor64();

  std::cout << "Please select number of human players:" << std::endl
            << " 0 - Computer vs Computer / Auto Demo" << std::endl
            << " 1 - Human vs Computer (Default)" << std::endl
            << " 2 - Human vs Human " << std::endl
            << "Enter your choice: ";

  int number_of_players = read_integer(1);
  if ((number_of_players < 0) || (number_of_players > 2)) {
    number_of_players = 1;
  }

  CellValue human_color_choice = CellValue::Blue;
  if (number_of_players == 1) {
    const int input_color_blue = 1;
    const int input_color_red = 2;
    std::cout << "Please select the color of the human player:" << std::endl
              << " " << input_color_blue
              << " - Blue - always going first (Default)" << std::endl
              << " " << input_color_red << " - Red " << std::endl
              << "Enter your choice: ";
    int player_color = read_integer(input_color_blue);
    if (player_color == input_color_red) {
      human_color_choice = CellValue::Red;
    }
  }

  std::cout << "Please select board size, i.e. the number of rows and columns; "
               "between "
            << min_board_size << " and " << max_board_size
            << " (default=" << default_board_size << "): ";

  // the board size, in the number of hexagonal cells in each row and column
  int board_size = read_integer(default_board_size);
  if ((board_size < min_board_size) || (board_size > max_board_size)) {
    board_size = default_board_size;
  }

  Field field(board_size,
              board_size); // the field keeps all the players' moves, in a
                           // 2-dimensional array of characters
  field.clear();

  int board_type;

  if (board_size <= 11) // for boards up to 11x11 you can select, but for larger
                        // boards only compact type is supported
  {
    std::cout << "Please select board type:" << std::endl
              << " " << input_board_type_parallelogram
              << " - Parallelogram comprised of hexagonal cells - conventional "
                 "(Default)"
              << std::endl
              << " " << input_board_type_rectangle
              << " - Rectangle comprised of hexagonal cells - nonstandard"
              << std::endl
              << " " << input_board_type_compact_grid << " - Compact grid"
              << std::endl
              << "Enter your choice: ";

    board_type = read_integer(input_board_type_default);
    if ((board_type < input_board_type_min) ||
        (board_type > input_board_type_max)) {
      board_type = input_board_type_default;
    }
  } else {
    board_type = input_board_type_compact_grid;
  }

  // configurable variables
  int cell_scale = 1; // the "scale factor" of a single cell as drawn in ASCII
                      // characters, e.g. 1 for smallest cell, 2 for a larger
                      // cell, 3 for even larger, etc ...

  // smaller the board, larger the cells
  if (board_size < 7)
    cell_scale++;
  if (board_size < 5)
    cell_scale++;

  randomize_xor64(); // call randomize again, after the user have made the
                     // choices

  ComputerMoveAbstract *computer_move_getter_ai;
  if (number_of_players <= 1) {
    // ask the class factory to supply us with the one tailored for the board
    // size
    computer_move_getter_ai = get_computer_move_ai_factory(board_size, field);
    assert(computer_move_getter_ai != nullptr);
  } else {
    computer_move_getter_ai = nullptr;
  }

  Turn turn = 0; // The counter of the number of moves all players have done;
                 // odd means that it is the first player's turn, even mean the
                 // second player's turn

  // the main loop of the game
  while (true) {
    turn++;
    BoardIncline
        board_incline; // the board incline will be set by the subseqent
                       // print_field call, uninitialized for now during the
                       // first iteration of the "while" loop
    print_field(field, board_incline, cell_scale, board_type);

    std::string player_color, player_name, player_name_and_color;
    CellValue player_mark;
    bool start_stop_edges_vertical;
    if (is_odd_turn(turn)) {
      player_color = "blue";
      player_name = "first";
      player_mark = CellValue::Blue;
      start_stop_edges_vertical = true;
    } else {
      player_color = "red";
      player_name = "second";
      player_mark = CellValue::Red;
      start_stop_edges_vertical = false;
    }

    bool is_computer_move =
        (number_of_players == 0) ||
        ((number_of_players == 1) && (human_color_choice != player_mark));

    if (number_of_players == 1) {
      if (human_color_choice == player_mark)
        player_name = "human";
      else
        player_name = "computer";
    }

    player_name_and_color = player_name + " (" + player_color + ")";

    if (is_computer_move) {
      // Computer move

      unsigned int row;
      unsigned int column;
      computer_move_getter_ai->make_move(player_mark, board_incline, row,
                                         column);

      // verify the computer move
      {
        bool column_valid = column >= first_column_base;
        bool row_valid = row >= first_row_base;
        unsigned int ux = column_valid ? column - first_column_base : 0;
        unsigned int uy = row_valid ? row - first_row_base : 0;
        assert(column_valid && row_valid);
        assert(!(field.coord_out_of_range_xy(ux, uy)));
        CellValue &c = field.cells[ux][uy];
        assert(c == player_mark);
      }
    } else {

      // Human move

      while (true) {
        // Human move

        unsigned int row;
        unsigned int column;

        if (number_of_players == 2) {
          std::cout << "Input column (x coordinate) for the "
                    << player_name_and_color << " player (1-" << field.width
                    << "): ";
        } else {
          std::cout << "Input column (x coordinate) [1.." << field.width
                    << "]: ";
        }
        std::cin >> column;
        if (std::cin.fail()) {
          ignore_line();
          continue;
        }
        if (column < first_column_base)
          continue;

        if (number_of_players == 2) {
          std::cout << "Input row (y coordinate) for the "
                    << player_name_and_color << " player (1-" << field.height
                    << "): ";
        } else {
          std::cout << "Input row (y coordinate) [1.." << field.height << "]: ";
        }
        std::cin >> row;
        if (std::cin.fail()) {
          ignore_line();
          continue;
        }
        if (row < first_row_base)
          continue;

        // column >= first_column_base and row >= first_row_base are guaranteed
        // by the continue statements above
        unsigned int ux = column - first_column_base;
        unsigned int uy = row - first_row_base;
        if (field.coord_out_of_range_xy(ux, uy)) {
          std::cout << "Invalid coordinates entered (" << row << "," << column
                    << ") Row should be between 1 and " << field.width
                    << "; column should be between 1 and " << field.height
                    << "." << std::endl;
          continue;
        }

        CellValue &c = field.cells[ux][uy]; // initialize the reference to further
                                          // modify the field value
        if (c == CellValue::Blank) {
          c = player_mark;
          break;
        } else {
          std::cout << "This is not a valid move. The cell (" << column << ","
                    << row << ") is already occupied" << std::endl;
          continue;
        }
      }

      // End of human move
    }

    if (is_game_won(field, player_mark, start_stop_edges_vertical,
                    board_incline)) {
      // print the board with the winning position
      print_field(field, board_incline, cell_scale, board_type);

      // announce the winner and exit
      std::cout << "The " << player_name_and_color << " player wins!"
                << std::endl;

      break;
    } else {
      // in a computer demo, draw an empty line between the fields unless the
      // game is won and finished
      if (number_of_players < 1)
        std::cout << std::endl;
    }
  }

  if (computer_move_getter_ai) {
    delete computer_move_getter_ai;
    computer_move_getter_ai = nullptr;
  }

  return 0;
}
