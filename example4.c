//
// Algoritmos e Estruturas de Dados --- 2023/2024
//
// David Amorim - Dec 2023
//
// Graph EXAMPLE : Reading graph from file
//

#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "instrumentation.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "ERROR: No input graph file is provided!\n");
    exit(1);
  }

  char* filepath = argv[1];

  FILE* fp = fopen(filepath, "r");

  if (fp == NULL) {
    perror("ERROR: Could not open the file");
    exit(1);
  }

  Graph* g = GraphFromFile(fp);
  fclose(fp);

  GraphDisplay(g);

  InstrReset();
  Graph* cp = GraphCopy(g);
  InstrPrint();

  int valid = GraphCheckInvariants(g);
  printf("Is graph valid? %s\n", valid ? "Yes" : "No");

  GraphDestroy(&g);
  GraphDestroy(&cp);

  return 0;
}