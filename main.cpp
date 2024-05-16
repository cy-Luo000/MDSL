#include "Graph.h"

int main(int argc, char* argv[]){
    puts("\n-----------------------------------------------------------------------------------------");
	Graph *graph = new Graph(argv[1], atof(argv[2]));
#ifdef _TEST_
    printf("#Filename=%s\n#gamma=%lf\n",argv[1], atof(argv[2]));
#endif
    graph->read();
    graph->search();
    puts("-----------------------------------------------------------------------------------------\n");
    return 0;
}