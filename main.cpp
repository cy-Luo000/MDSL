#include "Graph.h"

int main(int argc, char* argv[]){
    puts("\n-----------------------------------------------------------------------------------------");
    switch (DSkind){
    case 1:
        Graph *graph = new Graph(argv[1], atof(argv[2]));
#ifdef _TEST_
    printf("#Filename=%s\n#gamma=%lf\n",argv[1], atof(argv[2]));
#endif
        graph->read();
        graph->search(DSkind);
        break;
    
    case 2:
          Graph *graph = new Graph(argv[1], atoi(argv[2]));
#ifdef _TEST_
    printf("#Filename=%s\n#K=%lf\n",argv[1], atoi(argv[2]));
#endif
        graph->read();
        graph->search(DSkind);
        break;
    }
	
    puts("-----------------------------------------------------------------------------------------\n");
    return 0;
}