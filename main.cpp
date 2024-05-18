#include "Graph.h"

int main(int argc, char* argv[]){
    puts("\n-----------------------------------------------------------------------------------------");
    ui DSkind = 2;//1: quasi-clique; 2: k-defective clique

    if (DSkind==1){
        Graph *graph = new Graph(argv[1], atof(argv[2]));
#ifdef _TEST_
    printf("#Filename=%s\n#gamma=%lf\n",argv[1], atof(argv[2]));
#endif
        graph->read();
        graph->search(DSkind);
    }
    else{
        Graph *graph = new Graph(argv[1], atoi(argv[2]));
#ifdef _TEST_
    printf("#Filename=%s\n#K=%u\n",argv[1], atoi(argv[2]));
#endif
        graph->read();
        graph->search(DSkind);
    }
    puts("-----------------------------------------------------------------------------------------\n");
    return 0;
}