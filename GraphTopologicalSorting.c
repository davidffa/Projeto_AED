//
// Algoritmos e Estruturas de Dados --- 2023/2024
//
// Topological Sorting
//

#include "GraphTopologicalSorting.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Graph.h"
#include "IntegersQueue.h"
#include "instrumentation.h"

struct _GraphTopoSort {
  int* marked;                     // Aux array
  unsigned int* numIncomingEdges;  // Aux array
  unsigned int* vertexSequence;    // The result
  int validResult;                 // 0 or 1
  unsigned int numVertices;        // From the graph
  Graph* graph;
};

// AUXILIARY FUNCTION
// Allocate memory for the struct
// And for its array fields
// Initialize all struct fields
//
static GraphTopoSort* _create(Graph* g) {
  assert(g != NULL);

  GraphTopoSort* p = (GraphTopoSort*)malloc(sizeof(GraphTopoSort));

  if (p == NULL) abort();

  p->numVertices = GraphGetNumVertices(g);

  p->marked = calloc(p->numVertices, sizeof(int));
  p->numIncomingEdges = calloc(p->numVertices, sizeof(unsigned int));
  p->vertexSequence = calloc(p->numVertices, sizeof(unsigned int));
  p->validResult = 0;
  p->graph = g;

  if (p->marked == NULL || p->numIncomingEdges == NULL || p->vertexSequence == NULL) abort();

  return p;
}

//
// Computing the topological sorting, if any, using the 1st algorithm:
// 1 - Create a copy of the graph
// 2 - Successively identify vertices without incoming edges and remove their
//     outgoing edges
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//
GraphTopoSort* GraphTopoSortComputeV1(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting

  Graph* copy = GraphCopy(g);
  // Computed sequence length
  unsigned int seq_len = 0;

  while (seq_len < topoSort->numVertices) {
    unsigned int v;
    // 1 if we can select a vertex with inDegree = 0, 0 otherwise
    int flag = 0;

    for (v = 0; v < GraphGetNumVertices(copy); ++v) {
      // Check if vertex belongs to the graph (marked == 0) and if its incoming edges are equal to 0
      if (!topoSort->marked[v] && GraphGetVertexInDegree(copy, v) == 0) {
        flag = 1;
        break;
      }
    }

    // Could not select a vertex, stop the loop
    if (!flag) break;

    // Add the vertex to the sorted sequence and increase its length
    topoSort->vertexSequence[seq_len++] = v;

    // Remove vertex from the graph by removing its edges and marking it
    topoSort->marked[v] = 1;
    unsigned int* adjs = GraphGetAdjacentsTo(copy, v);

    for (unsigned int i = 1; i <= adjs[0]; ++i) {
      GraphRemoveEdge(copy, v, adjs[i]);
    }

    free(adjs);
  }

  // Cleanup
  GraphDestroy(&copy);

  // Check if we could build a valid sequence
  if (seq_len == topoSort->numVertices) {
    topoSort->validResult = 1;
  }

  return topoSort;
}

//
// Computing the topological sorting, if any, using the 2nd algorithm
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//
GraphTopoSort* GraphTopoSortComputeV2(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting
  unsigned int v;

  // Registar num array auxiliar numEdgesPerVertex o InDegree de cada vértice
  for (v = 0; v < GraphGetNumVertices(g); v++) {
    topoSort->numIncomingEdges[v] = GraphGetVertexInDegree(topoSort->graph, v);
  }

  unsigned int seq_len = 0;

  // Enquanto for possível
  while (seq_len < topoSort->numVertices) {
    int flag = 0;
    for (v = 0; v < GraphGetNumVertices(g); v++) {
      if (topoSort->numIncomingEdges[v] == 0 && !topoSort->marked[v]) {
        flag = 1;
        break;
      }
    }

    // Could not select a vertex, stop the loop
    if (!flag) break;

    // Imprimir o seu ID
    topoSort->vertexSequence[seq_len++] = v;

    // Marcar o vertice como pertencente à ordenação
    topoSort->marked[v] = 1;

    // Para cada vértice w adjacente a v
    unsigned int* w = GraphGetAdjacentsTo(g, v);

    for (unsigned int i = 1; i <= w[0]; i++) {
      topoSort->numIncomingEdges[w[i]]--;
    }

    free(w);
  }

  // Check if we could build a valid sequence
  if (seq_len == topoSort->numVertices) {
    topoSort->validResult = 1;
  }

  return topoSort;
}

//
// Computing the topological sorting, if any, using the 3rd algorithm
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//
GraphTopoSort* GraphTopoSortComputeV3(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting
  unsigned int v;
  int numEdgesPerVertex[GraphGetNumVertices(g)];
  // Registar num array auxiliar numEdgesPerVertex o InDegree de cada vértice
  for (v = 0; v < GraphGetNumVertices(g); v++) {
    numEdgesPerVertex[v] = GraphGetVertexInDegree(topoSort->graph, v);
  }
  // Criar FILA vazia e inserir na FILA os vértices v com numEdgesPerVertex[v] == 0
  Queue* FILA = QueueCreate(1);
  QueueClear(FILA);
  for (v = 0; v < GraphGetNumVertices(g); v++) {
    if (!QueueIsFull(FILA)) {
      if (numEdgesPerVertex[v] == 0) {
        QueueEnqueue(FILA, v);
      }
    }
  }
  unsigned int seq_len = 0;
  // Enquanto a FILA não for vazia
  while (!QueueIsEmpty(FILA)) {
    // v = retirar próximo vértice da FILA
    v = QueueDequeue(FILA);
    // Imprimir o seu ID
    topoSort->vertexSequence[seq_len++] = v;
    // Para cada vértice w adjacente a v
    unsigned int* w = GraphGetAdjacentsTo(g, v);
    for (unsigned int i = 0; i < w[0]; i++) {
      numEdgesPerVertex[w[i + 1]]--;
      // Se numEdgesPerVertex[w] == 0 Então Inserir w na FILA
      if (!QueueIsFull(FILA)) {
        if (numEdgesPerVertex[w[i + 1]] == 0) {
          QueueEnqueue(FILA, w[i + 1]);
        }
      }
    }
  }

  // Check if we could build a valid sequence
  if (seq_len == topoSort->numVertices) {
    topoSort->validResult = 1;
  }
  QueueDestroy(&FILA);
  return topoSort;
}

void GraphTopoSortDestroy(GraphTopoSort** p) {
  assert(*p != NULL);

  GraphTopoSort* aux = *p;

  free(aux->marked);
  free(aux->numIncomingEdges);
  free(aux->vertexSequence);

  free(*p);
  *p = NULL;
}

//
// A valid sorting was computed?
//
int GraphTopoSortIsValid(const GraphTopoSort* p) { return p->validResult; }

//
// Getting an array containing the topological sequence of vertex indices
// Or NULL, if no sequence could be computed
// MEMORY IS ALLOCATED FOR THE RESULTING ARRAY
//
unsigned int* GraphTopoSortGetSequence(const GraphTopoSort* p) {
  assert(p != NULL);

  if (!p->validResult)
    return NULL;

  // Alloc memory for sequence array
  unsigned int* seq = malloc(p->numVertices * sizeof(unsigned int));
  if (seq == NULL) abort();

  // Copies the computed vertexSequence into the "seq" result array
  memcpy(seq, p->vertexSequence, p->numVertices * sizeof(unsigned int));

  return seq;
}

// DISPLAYING on the console

//
// The toplogical sequence of vertex indices, if it could be computed
//
void GraphTopoSortDisplaySequence(const GraphTopoSort* p) {
  assert(p != NULL);

  if (p->validResult == 0) {
    printf(" *** The topological sorting could not be computed!! *** \n");
    return;
  }

  printf("Topological Sorting - Vertex indices:\n");
  for (unsigned int i = 0; i < GraphGetNumVertices(p->graph); i++) {
    printf("%d ", p->vertexSequence[i]);
  }
  printf("\n");
}

//
// The toplogical sequence of vertex indices, if it could be computed
// Followed by the digraph displayed using the adjecency lists
// Adjacency lists are presented in topologic sorted order
//
void GraphTopoSortDisplay(const GraphTopoSort* p) {
  assert(p != NULL);

  // The topological order
  GraphTopoSortDisplaySequence(p);

  if (p->validResult == 0) {
    return;
  }

  // The Digraph
  for (unsigned int i = 0; i < GraphGetNumVertices(p->graph); i++) {
    GraphListAdjacents(p->graph, p->vertexSequence[i]);
  }
  printf("\n");
}
