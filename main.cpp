#include "Graph.h"

int main(int argc, char* argv[]){
    puts("\n-----------------------------------------------------------------------------------------");
	Graph *graph = new Graph(argv[1], atof(argv[2]));
    graph->read();
    graph->search();
    puts("-----------------------------------------------------------------------------------------\n");
    return 0;
}