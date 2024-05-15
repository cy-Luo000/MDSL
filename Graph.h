#ifndef _GRAPH_H_
#define _GRAPH_H_
#include "Timer.h"
#include "LinearHeap.h"

class Graph
{
private:
    std::string dir;
    ui n;
    ept m;
    double gamma;

    ept *pstart;
    ui *edges;
    ui* degree;
    ept *pstart2Hp;
    ui *neis2Hp;
    ui *deg2Hp;

    ui maxDeg; //The max degree of G
    ui max2HpDeg; // The max 2hop degree of G
    ui maxCore; // The max core(degeneracy) of G
    ui max2HpCore;  // The max 2hop core(weak degeneracy) of G
    ui maxSub;

    std::vector<ui> MDS;
public:
    Graph(const char *_dir, const double _GAMMA);
    ~Graph();
    void read();
    void search();
    void degen(ui n, ui *seq, ui *core, ept *pstart, ui *edges, ui* degree, bool output);
    void build2Hpdeg(ept* pstart, ui* edges, ui* deg2Hp, char* deleted);
    void get2HpNei(ept* pstart, ui* edges, ui u, std::vector<int>& u_2HpNei, char* vis,char* deleted);
    void weakdegen(ui n, ui* seq, ui* core2Hp, ept* pstart, ui* edges, ui *deg2Hp, bool output);
    void induceSubgraph(ui u, ui* seq, ept* pstart, ui* edges, char* deleted, char* exists,ui* ids, ui* rid, ui& sub_n, std::vector<std::pair<ui,ui>>& vp);
    //input the sequence, the origin graph and the deleted vertices; ouput the size of the subgraph and the edges
};






#endif