#ifndef _QUASI_CLQUE_BB_
#define _QUASI_CLQUE_BB_

#include "Utility.h"
#include "Timer.h"
using namespace std;

long long treeCnt=0;
class QuasiClique_BB
{
private:
    int n;
    int m;
    double gamma;
    int UB;
    int P_end;int C_end;
    int MEInP, MEInG;
    int *pstart;
    int *edges;
    // bool *adjmtx;
    bool* matrix;

    int *PC;
    int *PC_rid;
    int *neiInG;
    int *neiInP;
    

    vector<int> QC;// current best solution
    int LB;// current best size
public:
    QuasiClique_BB();
    ~QuasiClique_BB();
    void load_graph(int _n,int *_pstart, int *_pend, int *_edges);
    void load_subgraph(double _gamma, int _n, vector<pair<int,int>> &_vp, vector<int> &_QC, int _UB);
    void MQCSearch(double _gamma, int _UB, std::vector<int> &_QC);
    void MQCSearch2hop();
    void branch(int level);
    void CtoP(int u, int level);
    void PtoC(int u, int level);
    void PtoX(int u, int level);
    void XtoC(int u, int level);
    // bool verifyQC();
    void store(int newLB);
    void swapID(int i, int j);
    bool inP(int u);
    bool inC(int u);
    bool inX(int u);
    bool &isAdj(int u, int v);
};

QuasiClique_BB::QuasiClique_BB(){
    n=0; m=0; gamma=0.0;
    MEInP=MEInG=0;
    pstart=NULL;
    edges=NULL;
    matrix=NULL;
    PC=PC_rid=NULL;
    neiInG=neiInP=NULL;
    P_end=C_end=0;
    QC.clear();
    LB=0; UB=0;
    
}

QuasiClique_BB::~QuasiClique_BB(){
    if(pstart!=NULL){
        delete[] pstart;
        pstart=NULL;
    }
    if (edges!=NULL){
        delete[] edges;
        edges=NULL;
    }
    if(PC!=NULL){
        delete[] PC;
        PC=NULL;
    }
    if(PC_rid!=NULL){
        delete[] PC_rid;
        PC_rid=NULL;
    }
    if(neiInP!=NULL){
        delete[] neiInP;
        neiInP=NULL;
    }
    if(neiInG!=NULL){
        delete[] neiInG;
        neiInG=NULL;
    }
    if(!QC.empty()){
        QC.clear();
    }
}

void QuasiClique_BB::load_graph(int _n,int *_pstart, int *_pend, int *_edges){
    n=_n;
    C_end=n;
    //m is initially zero
    int MaxDeg=0;
    for (int i = 0; i < n; i++) m+=_pend[i]-_pstart[i];
    assert(pstart==NULL);
    pstart=new int[n+1]; edges=new int[m];
    neiInP=new int[n];neiInG=new int[n];
    PC=new int[n]; PC_rid=new int[n];
    m=0;
    memset(neiInP,0,n*sizeof(int));
    for (int i = 0; i < n; i++){
        PC[i]=PC_rid[i]=i;
        pstart[i]=m;
        neiInG[i]=_pend[i]-_pstart[i];
        MaxDeg=max(MaxDeg,neiInG[i]);
        for(int j = _pstart[i];j < _pend[i];j ++) edges[m ++] = _edges[j];
    }
    //renew the missing edges in G
    MEInG=n*(n-1)/2-m/2;
    printf("load graph of size n=%u, m=%u (undirected), density=%.5lf, max degree=%d\n", n, m/2, double(m)/n/(n-1), MaxDeg);
}

void QuasiClique_BB::load_subgraph(double _gamma, int _n, vector<pair<int,int>> &_vp, vector<int> &_QC, int _UB){
    gamma=_gamma;
    n=_n;
    UB=_UB;
    //initialize the current best solution
    QC.resize(_QC.size());
    for (int i = 0; i < QC.size(); i++) QC[i]=_QC[i];
    LB=QC.size();
    matrix=new bool[n*n];
    PC=new int[n];
    PC_rid=new int[n];
    pstart=new int[n+1];
    edges=new int[2*m];
    neiInP= new int[n];
    neiInG= new int[n];
    memset(matrix, false, (n*n)*sizeof(bool));
    memset(neiInG, 0, n*sizeof(int));

    for(auto pr:_vp){
        int u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
    }
    // m=0;
    int idx=0; 
    //construct the subgraph of pstart and edges
    for (int i = 0; i < n; i++){
        pstart[i]=idx;
        for (int j = 0; j < n; j++){
            if(isAdj(i,j)) edges[idx++]=j;
        }
        neiInG[i]=idx-pstart[i];
        PC[i]=i; PC_rid[i]=i;
    }
    pstart[n]=idx;
    m=idx/2;
}
void QuasiClique_BB::MQCSearch(double _gamma, int _UB, std::vector<int> &_QC){
    gamma=_gamma;
    UB=_UB;
    LB=_QC.size();
    //use branch and bound to search

    int u=PC[0];
    CtoP(u, 0);
    // printf("enter left branch\n");
    branch(1);
    // printf("exit left branch\n");
    PtoX(u,0);
    // printf("enter right branch\n");
    branch(1);
    // printf("exit right branch\n");
    XtoC(u, 0);
    if(LB>_QC.size()){
        //renew the best solution
        _QC.clear();
        for (int i = 0; i < LB; i++) _QC.push_back(QC[i]);
    }
}

void QuasiClique_BB::MQCSearch2hop(){

    int u=PC[0];
    CtoP(u, 0);
    printf("enter left branch\n");
    branch(1);
    PtoC(u,0);
}

void QuasiClique_BB::branch(int level){
    // if(treeCnt>80) exit(0);
    assert(level<=n+1);
    // printf("P: %d, C: %d, treeId: %lld, level: %d\n",P_end, C_end-P_end, treeCnt,level);
    int u=PC[P_end];
    // if(verifyQC()) 
    if(C_end <= LB) goto REC;
    if (P_end > LB && (double)( P_end*(P_end-1)-2*MEInP )/( P_end*(P_end-1) ) >= gamma ) store(P_end);
    if ((double)( C_end*(C_end-1)-2*MEInG )/( C_end*(C_end-1) ) >= gamma) {store(C_end); goto REC;}

    //if G[C] is a quasi clique
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return  
    // printf("enter break point\n");
    treeCnt++;
    CtoP(u,level);
    branch(level+1);// branch on adding the vertex u
    PtoX(u,level);
    branch(level+1);// branch on deleting the vertex u
    XtoC(u, level);

REC:
    return;
}

void QuasiClique_BB::PtoC(int u, int level){
    assert(inP(u));
    swapID(PC_rid[u],--P_end);
    int nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;//update the missing edges in P
    //update the neiInP of neighbors of u
    for (int j = pstart[u]; j < pstart[u+1]; j++){
        int v=edges[j];
        neiInP[v]--;
    }
    
}

void QuasiClique_BB::CtoP(int u, int level){
    assert(inC(u));
    swapID(PC_rid[u], P_end++);
    int nonNeiP=P_end-1-neiInP[u];
    MEInP+=nonNeiP;
    //update the neiInP of neighbors of u
    for (int j = pstart[u]; j < pstart[u+1]; j++){
        int v=edges[j];
        neiInP[v]++;
    }
}

void QuasiClique_BB::PtoX(int u, int level){
    assert(inP(u));
    //move u from P to C
    swapID(PC_rid[u],--P_end);
    int nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;
    //move u from C to X
    swapID(P_end, --C_end);
    int nonNeiG=C_end-neiInG[u];
    MEInG-=nonNeiG;
    //update the neiInP and neiInG of neighbors of u
    for (int j = pstart[u]; j < pstart[u+1]; j++){
        int v=edges[j];
        neiInP[v]--;
        neiInG[v]--;
    }
}

void QuasiClique_BB::XtoC(int u, int level){
    assert(inX(u));
    //move u from X to C
    swapID(PC_rid[u],C_end++);
    int nonNeiG=C_end-1-neiInG[u];
    MEInG+=nonNeiG;
     //update the neiInG of neighbors of u
    for (int j = pstart[u]; j < pstart[u+1]; j++){
        int v=edges[j];
        neiInG[v]++;
    }
}

void QuasiClique_BB::store(int newLB){
    QC.resize(LB=newLB);
    for (int i = 0; i < LB; i++) QC[i]=PC[i];
}

void QuasiClique_BB::swapID(int i,int j){
    swap(PC[i],PC[j]);
    PC_rid[PC[i]] = i;
	PC_rid[PC[j]] = j;
}

bool QuasiClique_BB::inP(int u){
    return (PC_rid[u]>=0 && PC_rid[u]<P_end);
}

bool QuasiClique_BB::inC(int u){
    return (PC_rid[u]>=P_end && PC_rid[u]<=C_end);
}

bool QuasiClique_BB::inX(int u){
    return (PC_rid[u]>=C_end && PC_rid[u]<n);
}

bool& QuasiClique_BB::isAdj(int u,int v){
    return matrix[u*n+v];
}

#endif