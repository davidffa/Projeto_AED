//
// Algoritmos e Estruturas de Dados --- 2023/2024
//
// Joaquim Madeira, Joao Manuel Rodrigues - June 2021, Nov 2023
//
// Graph - Using a list of adjacency lists representation
//

#include "Graph.h"

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "SortedList.h"
#include "instrumentation.h"

struct _Vertex {
  unsigned int id;
  unsigned int inDegree;
  unsigned int outDegree;
  List* edgesList;
};

struct _Edge {
  unsigned int adjVertex;
  double weight;
};

struct _GraphHeader {
  int isDigraph;
  int isComplete;
  int isWeighted;
  unsigned int numVertices;
  unsigned int numEdges;
  List* verticesList;
};

#define VERTMEM InstrCount[0]
#define EDGEMEM InstrCount[1]
#define SUMS InstrCount[3]

// Used for instrumentation stuff
static bool init = false;

static void GraphInit(void) {
  InstrCalibrate();
  InstrName[0] = "vertmem";
  InstrName[1] = "edgemem";
  InstrName[2] = "qinsert";
  InstrName[3] = "sums";
  InstrName[4] = "listmove";
}

// The comparator for the VERTICES LIST

int graphVerticesComparator(const void* p1, const void* p2) {
  unsigned int v1 = ((struct _Vertex*)p1)->id;
  unsigned int v2 = ((struct _Vertex*)p2)->id;
  int d = v1 - v2;
  return (d > 0) - (d < 0);
}

// The comparator for the EDGES LISTS

int graphEdgesComparator(const void* p1, const void* p2) {
  unsigned int v1 = ((struct _Edge*)p1)->adjVertex;
  unsigned int v2 = ((struct _Edge*)p2)->adjVertex;
  int d = v1 - v2;
  return (d > 0) - (d < 0);
}

Graph* GraphCreate(unsigned int numVertices, int isDigraph, int isWeighted) {
  // Initialize the instrumentation
  if (!init) {
    GraphInit();
    init = true;
  }

  Graph* g = (Graph*)malloc(sizeof(struct _GraphHeader));
  if (g == NULL) abort();

  g->isDigraph = isDigraph;
  g->isComplete = 0;
  g->isWeighted = isWeighted;

  g->numVertices = numVertices;
  g->numEdges = 0;

  g->verticesList = ListCreate(graphVerticesComparator);

  for (unsigned int i = 0; i < numVertices; i++) {
    struct _Vertex* v = (struct _Vertex*)malloc(sizeof(struct _Vertex));
    if (v == NULL) abort();

    v->id = i;
    v->inDegree = 0;
    v->outDegree = 0;

    v->edgesList = ListCreate(graphEdgesComparator);

    ListInsert(g->verticesList, v);
  }

  assert(g->numVertices == ListGetSize(g->verticesList));

  return g;
}

Graph* GraphCreateComplete(unsigned int numVertices, int isDigraph) {
  Graph* g = GraphCreate(numVertices, isDigraph, 0);

  g->isComplete = 1;

  List* vertices = g->verticesList;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++) {
    struct _Vertex* v = ListGetCurrentItem(vertices);
    List* edges = v->edgesList;
    for (unsigned int j = 0; j < g->numVertices; j++) {
      if (i == j) {
        continue;
      }
      struct _Edge* new = (struct _Edge*)malloc(sizeof(struct _Edge));
      if (new == NULL) abort();
      new->adjVertex = j;
      new->weight = 1;

      ListInsert(edges, new);
    }
    if (g->isDigraph) {
      v->inDegree = g->numVertices - 1;
      v->outDegree = g->numVertices - 1;
    } else {
      v->outDegree = g->numVertices - 1;
    }
  }
  if (g->isDigraph) {
    g->numEdges = numVertices * (numVertices - 1);
  } else {
    g->numEdges = numVertices * (numVertices - 1) / 2;
  }

  return g;
}

void GraphDestroy(Graph** p) {
  assert(*p != NULL);
  Graph* g = *p;

  List* vertices = g->verticesList;
  if (ListIsEmpty(vertices) == 0) {
    ListMoveToHead(vertices);
    unsigned int i = 0;
    for (; i < g->numVertices; ListMoveToNext(vertices), i++) {
      struct _Vertex* v = ListGetCurrentItem(vertices);

      List* edges = v->edgesList;
      if (ListIsEmpty(edges) == 0) {
        unsigned int i = 0;
        ListMoveToHead(edges);
        for (; i < ListGetSize(edges); ListMoveToNext(edges), i++) {
          struct _Edge* e = ListGetCurrentItem(edges);
          free(e);
        }
      }
      ListDestroy(&(v->edgesList));
      free(v);
    }
  }

  ListDestroy(&(g->verticesList));
  free(g);

  *p = NULL;
}

/// Makes a Deep Copy of a graph
Graph* GraphCopy(const Graph* g) {
  assert(g != NULL);

  Graph* copy = (Graph*)malloc(sizeof(struct _GraphHeader));
  if (copy == NULL) abort();

  // Copy graph attributes
  copy->isDigraph = g->isDigraph;
  copy->isComplete = g->isComplete;
  copy->isWeighted = g->isWeighted;

  copy->numVertices = g->numVertices;
  copy->numEdges = g->numEdges;

  copy->verticesList = ListCreate(graphVerticesComparator);

  ListMoveToHead(g->verticesList);

  // Copy graph vertices
  for (unsigned int i = 0; i < copy->numVertices; ListMoveToNext(g->verticesList), i++) {
    struct _Vertex* copy_v = (struct _Vertex*)malloc(sizeof(struct _Vertex));
    if (copy_v == NULL) abort();

    // Load current vertex
    VERTMEM += 1;
    struct _Vertex* v = ListGetCurrentItem(g->verticesList);

    // Ensure that is the same vertex
    assert(v != NULL && v->id == i);

    // Copy vertex attributes
    copy_v->id = i;
    copy_v->inDegree = v->inDegree;
    copy_v->outDegree = v->outDegree;
    copy_v->edgesList = ListCreate(graphEdgesComparator);

    // Store copied vertex
    VERTMEM += 1;
    ListInsert(copy->verticesList, copy_v);

    // Copy graph edges
    List* edges = v->edgesList;
    ListMoveToHead(edges);
    for (unsigned int j = 0; j < ListGetSize(edges); ListMoveToNext(edges), j++) {
      struct _Edge* copy_e = (struct _Edge*)malloc(sizeof(struct _Edge));
      if (copy_e == NULL) abort();

      // Load edge
      EDGEMEM += 1;
      struct _Edge* edge = ListGetCurrentItem(edges);
      copy_e->adjVertex = edge->adjVertex;
      copy_e->weight = edge->weight;

      // Store edge
      EDGEMEM += 1;
      ListInsert(copy_v->edgesList, copy_e);
    }

    // Ensure that all edges were inserted successfully
    assert(copy_v->outDegree == ListGetSize(copy_v->edgesList));
  }

  // Ensure that all vertices were inserted successfully
  assert(copy->numVertices == ListGetSize(copy->verticesList));
  assert(GraphCheckInvariants(copy));

  return copy;
}

/// Constructs a Graph structure from a file with graph data
/// Graph file format:
/// 0 / 1 Is it directed ?
/// 0 / 1 Is it a weighted graph ?
/// Number of vertices
/// Number of edges
/// src vertex dest vertex weight (if it exists)
Graph* GraphFromFile(FILE* f) {
  assert(f != NULL);

  int directed, weighted;
  unsigned int numVertices, numEdges;

  // Read first four lines with numbers
  if (fscanf(f, "%d%d%u%u", &directed, &weighted, &numVertices, &numEdges) != 4) {
    fprintf(stderr, "ERROR: Invalid graph file format!\n");
    exit(1);
  }

  if (directed != 0 && directed != 1) {
    fprintf(stderr, "ERROR: Invalid graph file format!\n");
    fprintf(stderr, "Invalid flag isDirected, expected 0 or 1, but found %d\n", directed);
    exit(1);
  }

  if (weighted != 0 && weighted != 1) {
    fprintf(stderr, "ERROR: Invalid graph file format!\n");
    fprintf(stderr, "Invalid flag isWeighted, expected 0 or 1, but found %d\n", weighted);
    exit(1);
  }

  Graph* g = GraphCreate(numVertices, directed, weighted);

  double weight;
  unsigned int v, w;

  // We assume that the file contains at least numEdges+4 lines, so the fscanf should read the data successfully
  for (unsigned int i = 0; i < numEdges; ++i) {
    if (weighted) {
      if (fscanf(f, "%u %u %lf", &v, &w, &weight) != 3) {
        fprintf(stderr, "ERROR: Invalid graph file format! Could not read an weighted edge.\n");
        exit(1);
      }

      // The graph cannot have self-loops
      if (v == w) continue;

      GraphAddWeightedEdge(g, v, w, weight);
    } else {
      if (fscanf(f, "%u %u", &v, &w) != 2) {
        fprintf(stderr, "ERROR: Invalid graph file format! Could not read an edge.\n");
        exit(1);
      }

      // The graph cannot have self-loops
      if (v == w) continue;

      GraphAddEdge(g, v, w);
    }

    // NOTE: This Graph implementation does not allow parallel edges, but we don't need to
    //       check for that case because the ListInsert only allows inserting unique members,
    //       so that case will be handled "automatically" by ListInsert function
  }

  if (g->isDigraph) {
    if (g->numEdges == numVertices * (numVertices - 1))
      g->isComplete = 1;
  } else {
    if (g->numEdges == numVertices * (numVertices - 1) / 2)
      g->isComplete = 1;
  }

  assert(GraphCheckInvariants(g));

  return g;
}

// Graph

int GraphIsDigraph(const Graph* g) { return g->isDigraph; }

int GraphIsComplete(const Graph* g) { return g->isComplete; }

int GraphIsWeighted(const Graph* g) { return g->isWeighted; }

unsigned int GraphGetNumVertices(const Graph* g) { return g->numVertices; }

unsigned int GraphGetNumEdges(const Graph* g) { return g->numEdges; }

//
// For a graph
//
double GraphGetAverageDegree(const Graph* g) {
  assert(g->isDigraph == 0);
  return 2.0 * (double)g->numEdges / (double)g->numVertices;
}

static unsigned int _GetMaxDegree(const Graph* g) {
  List* vertices = g->verticesList;
  if (ListIsEmpty(vertices)) return 0;

  unsigned int maxDegree = 0;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++) {
    struct _Vertex* v = ListGetCurrentItem(vertices);
    if (v->outDegree > maxDegree) {
      maxDegree = v->outDegree;
    }
  }
  return maxDegree;
}

//
// For a graph
//
unsigned int GraphGetMaxDegree(const Graph* g) {
  assert(g->isDigraph == 0);
  return _GetMaxDegree(g);
}

//
// For a digraph
//
unsigned int GraphGetMaxOutDegree(const Graph* g) {
  assert(g->isDigraph == 1);
  return _GetMaxDegree(g);
}

// Vertices

//
// returns an array of size (outDegree + 1)
// element 0, stores the number of adjacent vertices
// and is followed by indices of the adjacent vertices
//
unsigned int* GraphGetAdjacentsTo(const Graph* g, unsigned int v) {
  assert(v < g->numVertices);

  // Node in the list of vertices
  List* vertices = g->verticesList;
  ListMove(vertices, v);

  // Load vertex
  VERTMEM += 1;
  struct _Vertex* vPointer = ListGetCurrentItem(vertices);
  unsigned int numAdjVertices = vPointer->outDegree;

  unsigned int* adjacent =
      (unsigned int*)calloc(1 + numAdjVertices, sizeof(unsigned int));

  if (numAdjVertices > 0) {
    adjacent[0] = numAdjVertices;
    List* adjList = vPointer->edgesList;
    ListMoveToHead(adjList);
    for (unsigned int i = 0; i < numAdjVertices; ListMoveToNext(adjList), i++) {
      // Load edge
      EDGEMEM += 1;
      struct _Edge* ePointer = ListGetCurrentItem(adjList);
      adjacent[i + 1] = ePointer->adjVertex;
    }
  }

  return adjacent;
}

//
// returns an array of size (outDegree + 1)
// element 0, stores the number of adjacent vertices
// and is followed by the distances to the adjacent vertices
//
double* GraphGetDistancesToAdjacents(const Graph* g, unsigned int v) {
  assert(v < g->numVertices);

  // Node in the list of vertices
  List* vertices = g->verticesList;
  ListMove(vertices, v);
  struct _Vertex* vPointer = ListGetCurrentItem(vertices);
  unsigned int numAdjVertices = vPointer->outDegree;

  double* distance = (double*)calloc(1 + numAdjVertices, sizeof(double));

  if (numAdjVertices > 0) {
    distance[0] = numAdjVertices;
    List* adjList = vPointer->edgesList;
    ListMoveToHead(adjList);
    for (unsigned int i = 0; i < numAdjVertices; ListMoveToNext(adjList), i++) {
      struct _Edge* ePointer = ListGetCurrentItem(adjList);
      distance[i + 1] = ePointer->weight;
    }
  }

  return distance;
}

//
// For a graph
//
unsigned int GraphGetVertexDegree(Graph* g, unsigned int v) {
  assert(g->isDigraph == 0);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  struct _Vertex* p = ListGetCurrentItem(g->verticesList);

  return p->outDegree;
}

//
// For a digraph
//
unsigned int GraphGetVertexOutDegree(Graph* g, unsigned int v) {
  assert(g->isDigraph == 1);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  struct _Vertex* p = ListGetCurrentItem(g->verticesList);

  return p->outDegree;
}

//
// For a digraph
//
unsigned int GraphGetVertexInDegree(Graph* g, unsigned int v) {
  assert(g->isDigraph == 1);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  // Load vertex
  VERTMEM += 1;
  struct _Vertex* p = ListGetCurrentItem(g->verticesList);

  return p->inDegree;
}

// Edges

static int _addEdge(Graph* g, unsigned int v, unsigned int w, double weight) {
  struct _Edge* edge = (struct _Edge*)malloc(sizeof(struct _Edge));
  edge->adjVertex = w;
  edge->weight = weight;

  ListMove(g->verticesList, v);
  struct _Vertex* vertex = ListGetCurrentItem(g->verticesList);
  int result = ListInsert(vertex->edgesList, edge);

  if (result == -1) {
    free(edge);
    return 0;
  } else {
    g->numEdges++;
    vertex->outDegree++;

    ListMove(g->verticesList, w);
    struct _Vertex* destVertex = ListGetCurrentItem(g->verticesList);
    destVertex->inDegree++;
  }

  if (g->isDigraph == 0) {
    // Bidirectional edge
    struct _Edge* edge = (struct _Edge*)malloc(sizeof(struct _Edge));
    edge->adjVertex = v;
    edge->weight = weight;

    ListMove(g->verticesList, w);
    struct _Vertex* vertex = ListGetCurrentItem(g->verticesList);
    result = ListInsert(vertex->edgesList, edge);

    if (result == -1) {
      free(edge);
      return 0;
    } else {
      // g->numEdges++; // Do not count the same edge twice on a undirected
      // graph !!
      vertex->outDegree++;
    }
  }

  return 1;
}

int GraphAddEdge(Graph* g, unsigned int v, unsigned int w) {
  assert(g->isWeighted == 0);
  assert(v != w);
  assert(v < g->numVertices);
  assert(w < g->numVertices);

  return _addEdge(g, v, w, 1.0);
}

int GraphAddWeightedEdge(Graph* g, unsigned int v, unsigned int w,
                         double weight) {
  assert(g->isWeighted == 1);
  assert(v != w);
  assert(v < g->numVertices);
  assert(w < g->numVertices);

  return _addEdge(g, v, w, weight);
}

int GraphRemoveEdge(Graph* g, unsigned int v, unsigned int w) {
  assert(g != NULL);

  // Get source vertex
  ListMove(g->verticesList, v);
  // Load vertex
  VERTMEM += 1;
  struct _Vertex* vertex = ListGetCurrentItem(g->verticesList);
  List* edges = vertex->edgesList;

  // Vertex v has no edges
  if (ListGetSize(edges) == 0) {
    return 0;
  }

  int num_edges = ListGetSize(edges);

  // Find the edge in v edges list associated with the vertex w
  for (ListMoveToHead(edges); ListGetCurrentIndex(edges) < num_edges; ListMoveToNext(edges)) {
    // Load edge
    EDGEMEM += 1;
    struct _Edge* edge = ListGetCurrentItem(edges);

    // If edge found, finish the loop, current element of edges is the edge we want.
    if (edge->adjVertex == w) break;
  }

  if (ListGetCurrentIndex(edges) == num_edges) {
    // Edge not found!
    return 0;
  }

  // Remove the edge from v edges list and deallocates it
  free(ListRemoveCurrent(edges));
  // Decrease outDegree of vertex
  vertex->outDegree--;
  SUMS += 1;

  // Remove the edge from dest vertex (w)
  ListMove(g->verticesList, w);
  VERTMEM += 1;  // Load vertex
  vertex = ListGetCurrentItem(g->verticesList);
  edges = vertex->edgesList;

  vertex->inDegree--;
  SUMS += 1;

  // NOTE: We don't add operations count below, because we only want to analyze the time complexity of topological sorting
  //       functions, which only applies to directed graphs
  if (!g->isDigraph) {
    // Vertex w has no edges
    if (ListGetSize(edges) == 0) {
      return 0;
    }

    num_edges = ListGetSize(edges);
    // Find the edge in w edges list associated with the vertex v
    for (ListMoveToHead(edges); ListGetCurrentIndex(edges) < num_edges; ListMoveToNext(edges)) {
      struct _Edge* edge = ListGetCurrentItem(edges);

      // If edge found, finish the loop, current element of edges is the edge we want.
      if (edge->adjVertex == v) break;
    }

    if (ListGetCurrentIndex(edges) == num_edges) {
      // Edge not found!
      return 0;
    }

    // Remove the edge from w edges list and deallocates it
    free(ListRemoveCurrent(edges));
    // Decrease outDegree of vertex w
    vertex->outDegree--;
  }

  // Decrease the graph edge count
  g->numEdges--;
  // Note: the addEdge function does not check if the graph became complete after adding the edge
  //       but in this Graph implementation we assume that the graph cannot have self-loops nor parallel edges
  g->isComplete = 0;

  assert(GraphCheckInvariants(g));
  return 1;
}

// CHECKING

int GraphCheckInvariants(const Graph* g) {
  assert(g != NULL);

  // Check number of vertices
  if (g->numVertices != ListGetSize(g->verticesList)) {
    return 0;
  }

  unsigned int numEdges = 0, inDegrees = 0;

  ListMoveToHead(g->verticesList);
  for (unsigned int i = 0; i < g->numVertices; ListMoveToNext(g->verticesList), ++i) {
    struct _Vertex* v = ListGetCurrentItem(g->verticesList);

    // Check degree
    if (v->outDegree != ListGetSize(v->edgesList)) {
      return 0;
    }

    numEdges += v->outDegree;
    inDegrees += v->inDegree;

    // Check for self-loops
    List* edges = v->edgesList;
    ListMoveToHead(edges);

    for (unsigned int j = 0; j < ListGetSize(edges); ListMoveToNext(edges), ++j) {
      struct _Edge* edge = ListGetCurrentItem(edges);

      // Graph cannot have self-loops
      if (edge->adjVertex == v->id) return 0;
    }
  }

  if (g->isDigraph) {
    // Number of edges (sum of outDegrees) must be equal to inDegrees sum
    if (numEdges != inDegrees) return 0;
    // Sum of outDegrees must be equal to number of edges
    if (numEdges != g->numEdges) return 0;
  } else {
    // Divide by 2 because two adj vertices have the same edge in its edgeList
    numEdges /= 2;

    if (numEdges != g->numEdges) return 0;
  }

  // Check complete graph invariants
  if (g->isComplete) {
    if (g->isDigraph) {
      if (g->numEdges != g->numVertices * (g->numVertices - 1))
        return 0;
    } else {
      if (g->numEdges != g->numVertices * (g->numVertices - 1) / 2)
        return 0;
    }
  }

  return 1;
}

// DISPLAYING on the console

void GraphDisplay(const Graph* g) {
  printf("---\n");
  if (g->isWeighted) {
    printf("Weighted ");
  }
  if (g->isComplete) {
    printf("COMPLETE ");
  }
  if (g->isDigraph) {
    printf("Digraph\n");
    printf("Max Out-Degree = %d\n", GraphGetMaxOutDegree(g));
  } else {
    printf("Graph\n");
    printf("Max Degree = %d\n", GraphGetMaxDegree(g));
  }
  printf("Vertices = %2d | Edges = %2d\n", g->numVertices, g->numEdges);

  List* vertices = g->verticesList;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++) {
    printf("%2d ->", i);
    struct _Vertex* v = ListGetCurrentItem(vertices);
    if (ListIsEmpty(v->edgesList)) {
      printf("\n");
    } else {
      List* edges = v->edgesList;
      unsigned int i = 0;
      ListMoveToHead(edges);
      for (; i < ListGetSize(edges); ListMoveToNext(edges), i++) {
        struct _Edge* e = ListGetCurrentItem(edges);
        if (g->isWeighted) {
          printf("   %2d(%4.2f)", e->adjVertex, e->weight);
        } else {
          printf("   %2d", e->adjVertex);
        }
      }
      printf("\n");
    }
  }
  printf("---\n");
}

void GraphListAdjacents(const Graph* g, unsigned int v) {
  printf("---\n");

  unsigned int* array = GraphGetAdjacentsTo(g, v);

  printf("Vertex %d has %d adjacent vertices -> ", v, array[0]);

  for (unsigned int i = 1; i <= array[0]; i++) {
    printf("%d ", array[i]);
  }

  printf("\n");

  free(array);

  printf("---\n");
}
