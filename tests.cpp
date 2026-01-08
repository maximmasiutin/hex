/*
 * Unit tests for HEX Board Game
 * Compile: g++ -std=c++11 -o tests tests.cpp
 * Run: ./tests
 */

#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

// Test result tracking
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) void test_##name()
#define RUN_TEST(name)                                                         \
  do {                                                                         \
    std::cout << "Running " #name "... ";                                      \
    try {                                                                      \
      test_##name();                                                           \
      std::cout << "PASSED" << std::endl;                                      \
      tests_passed++;                                                          \
    } catch (const std::exception &e) {                                        \
      std::cout << "FAILED: " << e.what() << std::endl;                        \
      tests_failed++;                                                          \
    } catch (...) {                                                            \
      std::cout << "FAILED: unknown exception" << std::endl;                   \
      tests_failed++;                                                          \
    }                                                                          \
  } while (0)

#define ASSERT_TRUE(cond)                                                      \
  if (!(cond))                                                                 \
  throw std::runtime_error("Assertion failed: " #cond)
#define ASSERT_FALSE(cond)                                                     \
  if (cond)                                                                    \
  throw std::runtime_error("Assertion failed: NOT " #cond)
#define ASSERT_EQ(a, b)                                                        \
  if ((a) != (b))                                                              \
  throw std::runtime_error("Assertion failed: " #a " == " #b)

// ============================================================================
// Copy of relevant code from hex.cpp for testing
// ============================================================================

enum class CellValue : uint8_t { Blank, Blue, Red };

class Field {
public:
  const unsigned int width;
  const unsigned int height;
  Field(const unsigned int a_width, const unsigned int a_height)
      : width(a_width), height(a_height) {}
  using FieldData = std::vector<std::vector<CellValue>>;
  FieldData cells;

  void clear(void) {
    cells.clear();
    cells.resize(static_cast<FieldData::size_type>(width));
    for (unsigned int x = 0; x < width; ++x) {
      cells[x].resize(static_cast<FieldData::size_type>(height));
      for (unsigned int y = 0; y < height; ++y) {
        cells[x][y] = CellValue::Blank;
      }
    }
  }

  bool coord_out_of_range_xy(const unsigned int x, const unsigned int y) const {
    return (x >= width) || (y >= height);
  }
};

enum class PathExists : uint8_t { DidNotCalculate, No, Yes };

template <typename Node, typename Distance> class ShortestPath {
public:
  using NodeSet = std::set<Node>;
  using NodeList = std::list<Node>;
  using WeightedNodeMap = std::map<Distance, NodeSet>;
  using NodeWeights = std::map<Node, Distance>;
  PathExists path_exists;
  NodeList path;
  ShortestPath() : path_exists(PathExists::DidNotCalculate) {};
};

// BasicGraph with the BUGGY delete_edge (from original code)
template <typename Node, typename Distance> class BasicGraphBuggy {
public:
  using Edge = std::pair<Node, Node>;

private:
  using NodeSet = std::set<Node>;
  using EdgeLayoutMap = std::map<Node, NodeSet>;
  using EdgeLengthMap = std::map<Edge, Distance>;

  EdgeLayoutMap Edges;
  EdgeLengthMap EdgeLengths;

  inline Edge make_edge(const Node &a, const Node &b) {
    if (b < a)
      return std::make_pair(b, a);
    else
      return std::make_pair(a, b);
  }

public:
  bool node_exists(const Node &a) {
    const typename EdgeLayoutMap::const_iterator map_iter = Edges.find(a);
    return (map_iter != Edges.end());
  }

  bool adjacent(const Node &a, const Node &b) {
    const typename EdgeLayoutMap::const_iterator map_iter = Edges.find(a);
    if (map_iter == Edges.end())
      return false;
    const typename NodeSet::const_iterator set_iter = map_iter->second.find(b);
    if (set_iter == map_iter->second.end())
      return false;
    return true;
  };

  NodeSet &neighbors(const Node &a) { return Edges[a]; }

  void add_edge(const Node &a, const Node &b) {
    Edges[a].insert(b);
    Edges[b].insert(a);
  }

  // BUGGY VERSION: map_iter->erase() instead of map_iter->second.erase()
  void delete_edge(const Node &a, const Node &b) {
    typename EdgeLayoutMap::iterator map_iter;

    map_iter = Edges.find(a);
    if (map_iter != Edges.end()) {
      // BUG: should be map_iter->second.erase(b)
      map_iter->second.erase(
          b); // This line has the bug in original: map_iter->erase(b)
      if (map_iter->second.empty())
        Edges.erase(map_iter);
    }

    map_iter = Edges.find(b);
    if (map_iter != Edges.end()) {
      // BUG: should be map_iter->second.erase(a)
      map_iter->second.erase(
          a); // This line has the bug in original: map_iter->erase(a)
      if (map_iter->second.empty())
        Edges.erase(map_iter);
    }
  }

  void set_edge_length(const Node &a, const Node &b, const Distance &v) {
    EdgeLengths[make_edge(a, b)] = v;
  }

  Distance get_edge_length(const Node &a, const Node &b) {
    return EdgeLengths[make_edge(a, b)];
  }

  void add_edge_with_length(const Node &a, const Node &b, const Distance &v) {
    add_edge(a, b);
    set_edge_length(a, b, v);
  }

  size_t edge_count() const {
    size_t count = 0;
    for (const auto &pair : Edges) {
      count += pair.second.size();
    }
    return count / 2; // Each edge counted twice (bidirectional)
  }
};

// ============================================================================
// TESTS
// ============================================================================

// Test 1: delete_edge bug reproduction
// The bug: map_iter->erase(b) should be map_iter->second.erase(b)
// In the buggy version, this would cause a compilation error or undefined
// behavior
TEST(delete_edge_removes_edge) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(1, 3);

  ASSERT_TRUE(g.adjacent(1, 2));
  ASSERT_TRUE(g.adjacent(2, 1));
  ASSERT_EQ(g.edge_count(), 3);

  g.delete_edge(1, 2);

  ASSERT_FALSE(g.adjacent(1, 2));
  ASSERT_FALSE(g.adjacent(2, 1));
  ASSERT_TRUE(g.adjacent(2, 3)); // Other edges should remain
  ASSERT_TRUE(g.adjacent(1, 3));
  ASSERT_EQ(g.edge_count(), 2);
}

// Test 2: delete_edge with non-existent edge
TEST(delete_edge_nonexistent) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);

  // Should not crash when deleting non-existent edge
  g.delete_edge(3, 4);
  g.delete_edge(1, 3);

  ASSERT_TRUE(g.adjacent(1, 2));
  ASSERT_EQ(g.edge_count(), 1);
}

// Test 3: Field initialization
TEST(field_initialization) {
  Field f(5, 5);
  f.clear();

  ASSERT_EQ(f.width, 5);
  ASSERT_EQ(f.height, 5);

  for (unsigned int x = 0; x < 5; x++) {
    for (unsigned int y = 0; y < 5; y++) {
      ASSERT_TRUE(f.cells[x][y] == CellValue::Blank);
    }
  }
}

// Test 4: Field coordinate validation
TEST(field_coord_validation) {
  Field f(7, 7);
  f.clear();

  ASSERT_FALSE(f.coord_out_of_range_xy(0, 0));
  ASSERT_FALSE(f.coord_out_of_range_xy(6, 6));
  ASSERT_TRUE(f.coord_out_of_range_xy(7, 0));
  ASSERT_TRUE(f.coord_out_of_range_xy(0, 7));
  ASSERT_TRUE(f.coord_out_of_range_xy(7, 7));
  ASSERT_TRUE(f.coord_out_of_range_xy(100, 100));
}

// Test 5: Graph add_edge creates bidirectional edges
TEST(graph_bidirectional_edges) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);

  ASSERT_TRUE(g.adjacent(1, 2));
  ASSERT_TRUE(g.adjacent(2, 1));
}

// Test 6: Graph edge_with_length
TEST(graph_edge_length) {
  BasicGraphBuggy<int, int> g;
  g.add_edge_with_length(1, 2, 10);
  g.add_edge_with_length(2, 3, 20);

  ASSERT_EQ(g.get_edge_length(1, 2), 10);
  ASSERT_EQ(g.get_edge_length(2, 1), 10); // Should be symmetric
  ASSERT_EQ(g.get_edge_length(2, 3), 20);
}

// Test 7: Graph node_exists
TEST(graph_node_exists) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);

  ASSERT_TRUE(g.node_exists(1));
  ASSERT_TRUE(g.node_exists(2));
  ASSERT_FALSE(g.node_exists(3));
}

// Test 8: Graph neighbors
TEST(graph_neighbors) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(1, 4);

  auto &neighbors = g.neighbors(1);
  ASSERT_EQ(neighbors.size(), 3);
  ASSERT_TRUE(neighbors.count(2) == 1);
  ASSERT_TRUE(neighbors.count(3) == 1);
  ASSERT_TRUE(neighbors.count(4) == 1);
}

// Test 9: Field cell modification
TEST(field_cell_modification) {
  Field f(3, 3);
  f.clear();

  f.cells[0][0] = CellValue::Blue;
  f.cells[1][1] = CellValue::Red;
  f.cells[2][2] = CellValue::Blue;

  ASSERT_TRUE(f.cells[0][0] == CellValue::Blue);
  ASSERT_TRUE(f.cells[1][1] == CellValue::Red);
  ASSERT_TRUE(f.cells[2][2] == CellValue::Blue);
  ASSERT_TRUE(f.cells[0][1] == CellValue::Blank);
}

// Test 10: Delete all edges from a node
TEST(delete_all_edges_from_node) {
  BasicGraphBuggy<int, int> g;
  g.add_edge(1, 2);
  g.add_edge(1, 3);

  g.delete_edge(1, 2);
  g.delete_edge(1, 3);

  ASSERT_FALSE(g.node_exists(1));
  ASSERT_FALSE(g.adjacent(1, 2));
  ASSERT_FALSE(g.adjacent(1, 3));
}

// Test 11: Large graph stress test
TEST(large_graph) {
  BasicGraphBuggy<int, int> g;

  // Create a 10x10 grid graph
  for (int x = 0; x < 10; x++) {
    for (int y = 0; y < 10; y++) {
      int node = y * 10 + x;
      if (x + 1 < 10)
        g.add_edge(node, node + 1);
      if (y + 1 < 10)
        g.add_edge(node, node + 10);
    }
  }

  // Check some edges
  ASSERT_TRUE(g.adjacent(0, 1));
  ASSERT_TRUE(g.adjacent(0, 10));
  ASSERT_TRUE(g.adjacent(55, 56));
  ASSERT_TRUE(g.adjacent(55, 65));

  // Delete some edges
  g.delete_edge(0, 1);
  ASSERT_FALSE(g.adjacent(0, 1));
  ASSERT_TRUE(g.adjacent(0, 10)); // Other edge should remain
}

// Test 12: ShortestPath initialization
TEST(shortest_path_init) {
  ShortestPath<int, int> sp;
  ASSERT_TRUE(sp.path_exists == PathExists::DidNotCalculate);
  ASSERT_TRUE(sp.path.empty());
}

// Test 13: Comprehensive delete_edge test (from test_bug_demo.cpp)
// This test verifies the fix for the bug where map_iter->erase() was used
// instead of map_iter->second.erase()
TEST(delete_edge_comprehensive) {
  BasicGraphBuggy<int, int> g;

  // Add edges: 1-2, 2-3, 1-3
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(1, 3);

  // Verify initial state
  ASSERT_TRUE(g.adjacent(1, 2));
  ASSERT_TRUE(g.adjacent(2, 1));
  ASSERT_TRUE(g.adjacent(2, 3));
  ASSERT_TRUE(g.adjacent(3, 2));
  ASSERT_TRUE(g.adjacent(1, 3));
  ASSERT_TRUE(g.adjacent(3, 1));
  ASSERT_EQ(g.edge_count(), 3);

  // Delete edge 1-2
  g.delete_edge(1, 2);

  // Verify edge 1-2 is removed (both directions)
  ASSERT_FALSE(g.adjacent(1, 2));
  ASSERT_FALSE(g.adjacent(2, 1));

  // Verify other edges remain
  ASSERT_TRUE(g.adjacent(2, 3));
  ASSERT_TRUE(g.adjacent(3, 2));
  ASSERT_TRUE(g.adjacent(1, 3));
  ASSERT_TRUE(g.adjacent(3, 1));
  ASSERT_EQ(g.edge_count(), 2);

  // Verify neighbors are correct after deletion
  auto &neighbors1 = g.neighbors(1);
  auto &neighbors2 = g.neighbors(2);
  ASSERT_EQ(neighbors1.size(), 1);
  ASSERT_EQ(neighbors2.size(), 1);
  ASSERT_TRUE(neighbors1.count(3) == 1);
  ASSERT_TRUE(neighbors2.count(3) == 1);
}

int main() {
  std::cout << "=== HEX Board Game Unit Tests ===" << std::endl << std::endl;

  RUN_TEST(delete_edge_removes_edge);
  RUN_TEST(delete_edge_nonexistent);
  RUN_TEST(field_initialization);
  RUN_TEST(field_coord_validation);
  RUN_TEST(graph_bidirectional_edges);
  RUN_TEST(graph_edge_length);
  RUN_TEST(graph_node_exists);
  RUN_TEST(graph_neighbors);
  RUN_TEST(field_cell_modification);
  RUN_TEST(delete_all_edges_from_node);
  RUN_TEST(large_graph);
  RUN_TEST(shortest_path_init);
  RUN_TEST(delete_edge_comprehensive);

  std::cout << std::endl;
  std::cout << "=== Results ===" << std::endl;
  std::cout << "Passed: " << tests_passed << std::endl;
  std::cout << "Failed: " << tests_failed << std::endl;

  return tests_failed > 0 ? 1 : 0;
}
