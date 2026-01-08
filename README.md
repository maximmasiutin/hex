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

Each legal move will be evaluated using a million (1,000,000) trials. See
the "monte_carlo_iterations" variable to configure this value.

Each trial winds the game forward by randomly selecting successive moves until
there is a winner, and the trial is counted as a win or loss. The ratio:
wins/trials is the AI's metric for picking which next move to make.

The Monte-Carlo simulation is implemented very efficiently, so it takes just
about a second to make a move from a million trials on an average notebook on
a 7x7 field, or about 5 seconds on an 11x11 field.

This is a pure Monte-Carlo implementation, without the min-max algorithm
or the alpha-beta pruning.

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
empty cell column (x) and row (y). After the user enters the move, the
program checks whether this is a legal move, and if not, the program asks the
user again to make a move.

## Building and Testing with Docker

Several Dockerfiles are provided for different purposes:

### Testing (syntax check and unit tests)

Build and run unit tests:

    docker build -f Dockerfile.test -t hex-test .
    docker run --rm hex-test

The test Dockerfile compiles with -O1 -Wall -Wextra flags to catch potential
issues.

### Playing the game

Pre-built images are available on DockerHub. Choose the version matching your
CPU:

**For newer Intel processors (Sapphire Rapids and later):**

    docker run -it --rm maximmasiutin/hex:sapphirerapids

This version uses GCC 15 and targets Intel Sapphire Rapids microarchitecture,
taking advantage of newer instruction sets including AVX-512.

**For older Intel processors (Skylake and later):**

    docker run -it --rm maximmasiutin/hex:skylake

This version uses GCC 12 and targets Intel Skylake microarchitecture.

Note: Use -it flags for interactive play (required for keyboard input).

### Building locally

To build the Docker images locally instead of pulling from DockerHub:

    docker build -f Dockerfile.sapphirerapids -t hex-sapphirerapids .
    docker build -f Dockerfile.skylake -t hex-skylake .
