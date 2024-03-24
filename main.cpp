#include "Graph.h"
int main(int argc, char *argv[]) {
	puts("\n-----------------------------------------------------------------------------------------");
	Graph *graph = new Graph(argv[1], atof(argv[2]));
// #ifdef _TEST_
	printf("#Filename=%s\n#K=%lf\n",argv[1], atof(argv[2]));
	// exit(0);
// #endif
	graph->read(); graph->setK();
	graph->print_info();
	// exit(0);
	graph->search(); 
	graph->write(); // delete graph;
	puts("-----------------------------------------------------------------------------------------\n");
}