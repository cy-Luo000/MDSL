#ifndef _QUASI_CLQUE_BB_
#define _QUASI_CLQUE_BB_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
using namespace std;

long long treeCnt=0;
long long ub_prune=0;
long long prune1=0;
class QuasiClique_BB;
class QuasiClique_BB2hop;

class HeuriSearcher{
public:
	int n;
	int m;
	int* pstart;
	int* edges;
    // int k; //total missing edges
	double gamma;
	int maxDeg;
    int LB; // the lowerbound
    int UB; // the upperbound
	
    int n_s, m_s;
	std::vector<int> id_s;
	std::vector<int> rid_s;
	std::vector<int> deg_s;
    std::vector<int> core_s;
    std::vector<int> seq_s;
	bool* exist_s;
    char* vis_s;
	std::vector<int> heuQC;
	std::vector<std::vector<int>> adjList;
    // ListLinearHeap *s_heap;
	// vector<int> &KDC;

    // int* PC;
    // int* PC_rid;// to record the positon in PC of each vertex v of G
    // int P_end;// the P is from PC[0] to PC[P_end-1]
    // int C_end;// the C is from PC[P_end] to PC[C_end-1]

    // int* neiInP;// to record every vertex's degree in P
    // int* neiInG;// to record every vertex's degree in C
    // int MEInP;// the missing edges in P 
    // int MEInG;// the missing edges in C

	// Your code
	HeuriSearcher(int _n,int _m, int *_pstart, int *_edges, double _gamma, int _maxDeg, std::vector<int>& _QC){
		n=_n, m=_m, gamma=_gamma;
		this->pstart=_pstart;
		this->edges=_edges;
		maxDeg=_maxDeg;
		exist_s=new bool[n];
        vis_s=new char[1+maxDeg];
		memset(exist_s, false, n*sizeof(bool));
        LB=int(_QC.size());
        
	}
 	void induceNei(int u){
		n_s=m_s=0;
		//1. get the max id
		int maxId=u;
		for (int i = pstart[u]; i < pstart[u+1]; i++) maxId=max(maxId, edges[i]);
		rid_s.resize(1+maxId);
		//2. get the veritces
		id_s.push_back(u); rid_s[u]=id_s.size()-1; exist_s[u]=true;
		
		for (int i = pstart[u]; i < pstart[u+1]; i++){
			id_s.push_back(edges[i]);
            exist_s[edges[i]]=true;
			rid_s[edges[i]]=id_s.size()-1;
		}
        n_s=id_s.size(); 
        core_s.resize(n_s);//core_s initialization
        assert(seq_s.empty());
        for (int i = 0; i < n_s; i++) seq_s.push_back(i);  
		//3. get the edges
		for(auto v: id_s) exist_s[v]=true;
		adjList.resize(int(id_s.size()));
		for (int i = 0; i < id_s.size(); i++){
			int v=id_s[i];
			for (int j = pstart[v]; j < pstart[v+1]; j++){
				int w=edges[j];
				if(exist_s[w]) {
					adjList[rid_s[v]].push_back(rid_s[w]);
				}
			}
		}
		//4. construct the degrees
		for (int i = 0; i < n_s; i++) deg_s.push_back(int(adjList[i].size())), m_s+=deg_s[i];
        m_s/=2;
        //5. clear the exist
        for (int i = 0; i < n_s; i++) exist_s[id_s[i]]=false;   
	}
    int degenHeu(int u_0, int u_deg, ListLinearHeap *heap){
        int subSz=1+u_deg, s_maxcore=0;
        int edgeCnt_s=m_s, idx=subSz;
        // printf("u0: %d, subSz: %d\n",u_0, subSz);
        memset(vis_s, 0, subSz*sizeof(char));
        heap->init(subSz, subSz-1, seq_s.data(), deg_s.data());
        //delete the vertex by degeneracy order
        for (int i = 0; i < subSz; i++){
            if(idx == subSz && edgeCnt_s >= ceil(gamma*(1.0*(subSz - i)*(subSz - i - 1)/2.0)) ) idx = i;
            int u, key; heap->pop_min(u,key);
            if(key> s_maxcore) s_maxcore=key; 
            core_s[u]=s_maxcore; vis_s[u]=1; seq_s[i]=u;
            for(auto v: adjList[u]){
                if(vis_s[v]==0) heap->decrement(v,1);
            }
            edgeCnt_s-=key;
        }
        if(subSz-idx>heuQC.size()){
            printf("renew heuQC\n");
            LB=subSz-idx;
            heuQC.clear();
            for (int i = idx; i < subSz; i++) heuQC.push_back(id_s[seq_s[i]]);
            if(!checkExist(u_0, heuQC)) {printf("error, not include u: %d\n", u_0);}//this should be useless because u_0 is adjacent to all vertices 
        }

        seq_s.clear();
        return n;
    }
	void clear(){
		id_s.clear(), rid_s.clear();
		for(auto edgeVec: adjList) edgeVec.clear();
		adjList.clear(); deg_s.clear();
	}
    bool checkExist(int u, std::vector<int> &_QC){
        bool flag=false;
        for(auto v:_QC){
            if(v==u) {flag=true;break;}
        }
        return flag;
    }
	void search(std::vector<int> &_QC){
        int max_n=1+maxDeg;
        ListLinearHeap* s_heap=new ListLinearHeap(max_n, max_n-1);
        for (int u = 0; u < n; u++){
            //1. build the subgraph
            if(u==0){
                printf("check\n");
            }
            // printf("heuristic search subgraph: %d\n",u);
            induceNei(u);
            int u_deg=adjList[rid_s[u]].size();
            this->UB=degenHeu(u, u_deg, s_heap);
            this->clear();
        }
        if(LB>int(_QC.size())){
            _QC.resize(LB);
            for (int i = 0; i < LB; i++){
                _QC[i]=heuQC[i];
            }
            printf("renew the result in heuristic search\n");
        }
	}
};

class QuasiClique_BB
{
private:
    int n;
    int m;
    int maxDeg;
    int minDeg;
    double gamma;
    int UB;
    int P_end;int C_end;
    int MEInP, MEInG;
    int *pstart;
    int *edges;
    long long treeIdx;
    // bool *adjmtx;
    bool* matrix;

    int *PC;
    int *PC_rid;
    int *neiInG;
    int *neiInP;
    //graph coloring
    long long *colUseMtx;//to record if color is used
    int *colorSz;// the size of each color bucket
    int *colorVec;// the vertex in each color set
    std::vector<std::vector<int>> nonNeiB;
    std::vector<int> weiB;
    std::vector<int> weiPreSum;//the prefix sum of weight bucket
    std::vector<int> QC;// current best solution
    
public:
    int LB;// current best size
    QuasiClique_BB();
    ~QuasiClique_BB();
    void load_graph(int _n,int *_pstart, int *_pend, int *_edges);
    void load_subgraph(double _gamma, int _n, vector<pair<int,int>> &_vp, vector<int> &_QC, int _UB);
    void printInfo();
    void MQCSearch(double _gamma, int _UB, std::vector<int> &_QC);
    void MQCSearch2hop(vector<int> &_QC);
    int sortBound();
    bool verifyQC();
    bool verify2hop(int _end);
    void branch(int level);
    int vertexChoose();
    bool vertexCmp(int u, int v);
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
    maxDeg=0, minDeg=0;
    MEInP=MEInG=0;
    pstart=NULL;
    edges=NULL;
    matrix=NULL;
    PC=PC_rid=NULL;
    neiInG=neiInP=NULL;
    colUseMtx=NULL, colorSz=NULL, colorVec=NULL;
    P_end=C_end=0;
    QC.clear();
    nonNeiB.clear(), weiB.clear();
    weiPreSum.clear();
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
    if(!nonNeiB.empty()){
        nonNeiB.clear();
    }
    if(!weiB.empty()){
        weiB.clear();
    }
    if(!weiPreSum.empty()){
        weiPreSum.clear();
    }
    if(colorSz!=NULL){
        delete[] colorSz;
        colorSz=NULL;
    }
    if(colUseMtx!=NULL){
        delete[] colUseMtx;
        colUseMtx=NULL;
    }
    if(colorVec!=NULL){
        delete[] colorVec;
        colorVec=NULL;
    }
    if(matrix!=NULL){
        delete[] matrix;
        matrix=NULL;
    }
}
bool QuasiClique_BB::vertexCmp(int u, int v){
    if(neiInP[u]>neiInP[v]) return true;
	if(neiInP[u]==neiInP[v] && neiInG[u]>neiInG[v]) return true;
	return false;
}
int QuasiClique_BB::vertexChoose(){
    int u=PC[P_end];
    for(int id=P_end+1;id<C_end;id++){
        int v=PC[id];
        if(vertexCmp(v,u))u=v;
    }
    return u;
}
void QuasiClique_BB::load_graph(int _n,int *_pstart, int *_pend, int *_edges){
    n=_n;
    C_end=n;
    //m is initially zero
    minDeg=n;
    for (int i = 0; i < n; i++) m+=_pend[i]-_pstart[i];
    assert(pstart==NULL);
    pstart=new int[n+1]; edges=new int[m];
    neiInP=new int[n];neiInG=new int[n];
    PC=new int[n]; PC_rid=new int[n];
    colUseMtx=new long long[n*n]; colorSz=new int[n]; 
    colorVec=new int[n];
    m=0;
    memset(neiInP,0,n*sizeof(int));
    memset(colUseMtx, 0, n*n*sizeof(long long));
    for (int i = 0; i < n; i++){
        PC[i]=PC_rid[i]=i;
        pstart[i]=m;
        neiInG[i]=_pend[i]-_pstart[i];
        maxDeg=max(maxDeg,neiInG[i]); minDeg=min(minDeg, neiInG[i]);
        for(int j = _pstart[i];j < _pend[i];j ++) edges[m ++] = _edges[j];
    }
    //renew the missing edges in G
    MEInG=n*(n-1)/2-m/2;
    printf("load graph of size n=%u, m=%u (undirected), density=%.5lf, max degree=%d\n", n, m/2, double(m)/(n*(n-1)), maxDeg);
}

void QuasiClique_BB::load_subgraph(double _gamma, int _n, vector<pair<int,int>> &_vp, vector<int> &_QC, int _UB){
    gamma=_gamma;
    n=_n;
    UB=_UB;
    minDeg=n;
    C_end=n;
    treeIdx=0;
    m=_vp.size();
    //initialize the current best solution
    // QC.resize(_QC.size());
    // for (int i = 0; i < QC.size(); i++) QC[i]=_QC[i];
    // QC=_QC;
    LB=_QC.size();
    matrix=new bool[n*n];
    PC=new int[n];
    PC_rid=new int[n];
    pstart=new int[n+1];
    edges=new int[2*m];
    neiInP= new int[n];
    neiInG= new int[n];
    colUseMtx=new long long[n*n]; colorSz=new int[n]; 
    colorVec=new int[n];
    memset(matrix, false, (n*n)*sizeof(bool));
    memset(neiInP, 0, n*sizeof(int));
    memset(colUseMtx, 0, n*n*sizeof(long long));
    for(auto pr:_vp){
        int u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
    }
    int idx=0; 
    //construct the subgraph of pstart and edges
    for (int i = 0; i < n; i++){
        pstart[i]=idx;
        for (int j = 0; j < n; j++){
            if(isAdj(i,j)) edges[idx++]=j;
        }
        neiInG[i]=idx-pstart[i];
        maxDeg=max(maxDeg, neiInG[i]), minDeg=min(minDeg, neiInG[i]);
        PC[i]=i; PC_rid[i]=i;
    }
    pstart[n]=idx;
    m=idx/2;
    MEInG=n*(n-1)/2-m;
}
void QuasiClique_BB::printInfo(){
    printf("vertex num: %d, edge num: %d\n",n,m);
    printf("max degree: %d, min degree: %d\n", maxDeg,minDeg);
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

bool QuasiClique_BB::verifyQC(){
    bool flag=false;
    int m_qc=0, n_qc=this->QC.size();
    double density=0.0;
    if(n_qc==0){
        printf("trivial gamma quasi clique\n");
        return true;
    }
    for (int i = 0; i < QC.size(); i++){
        for (int j = i+1; j < QC.size(); j++){
            if(isAdj(QC[i], QC[j])) m_qc++;
        }
    }
    if(2.0*m_qc>= gamma*(double)(n_qc)*(n_qc-1)){
        return true;
    }
    density=2.0*m_qc/((double)n_qc*(n_qc-1));
    printf("QC vNum: %d, QC eNum: %d, QC density: %.2f\n",n_qc, m_qc, density);
    return false;
}
bool QuasiClique_BB::verify2hop(int _end){
    bool flag=true, curflag=false;
    int u_0=PC[0];
    vector<int> nonNei_u0;
    for (int i = 1; i < _end; i++){
        int v=PC[i];
        if(!isAdj(u_0, PC[i])) nonNei_u0.push_back(v);
    }
    for (int i = 0; i < _end; i++){
        int u=PC[i];
        for (int j = 0; j < nonNei_u0.size(); j++){
            int v=nonNei_u0[j];
            if(u==v || isAdj(u,v)) continue;
            curflag=false;
            for (int k = 0; k < _end; k++){
                int w=PC[k];
                if(w==u || w==v) continue;
                if(isAdj(u,w) && isAdj(v,w)) {
                    curflag=true;
                    break;
                }
            }
            if(!curflag) {
                nonNei_u0.clear();
                return false;
            }
        }   
    }
    nonNei_u0.clear();
    return flag;
}
void QuasiClique_BB::MQCSearch2hop(vector<int> &_QC){
    Timer t;
    int u=PC[0];
    CtoP(u, 0);
    printf("enter left branch\n");
    branch(1);
    PtoC(u,0);
    if(LB>_QC.size()){
        printf("renew result\n");
        _QC.resize(LB);
        for (int i = 0; i < QC.size(); i++) _QC[i]=QC[i];
    }
    printf("subgraph search complete, search time: %.2f, treeCnt: %d\n",double(t.elapsed())/1000000, treeCnt);

}
int QuasiClique_BB::sortBound(){
    int UB=0; int maxWei=0;
    nonNeiB.resize(P_end+1);
    weiB.resize(C_end);
    weiPreSum.resize(C_end-P_end);
    int colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (int i = P_end; i < C_end; i++){
        int u= PC[i];
		// if(P_end-neiInP[i]>k) continue;
		nonNeiB[int(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    for (int i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            int col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (int j = pstart[u]; j < pstart[u+1]; j++){
                int nei=edges[j];
                // if(!isAdj(u,nei)) continue;
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorSz[col]++;
            int wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
        
    }
    //3. calculation of the prefix sum of the weight bucket 
    int vIdx=0;
    for (int w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
        // if(weiB[w]==0) continue;
    }
    //4. calculating the upper bound
    for (int i = C_end-1; i >= P_end; i--){
        if(i*(i+1)/2-MEInP- weiPreSum[i-P_end]>= gamma*i*(i+1)/2.0){
            UB=i+1;
            break;
        }
    }
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear();
    return UB;
}
void QuasiClique_BB::branch(int level){
    // if(treeCnt>80) exit(0);
    treeIdx++;
    assert(level<=n+1);
    // printf("P: %d, C: %d, treeId: %lld, level: %d\n",P_end, C_end-P_end, treeCnt,level);
    int u=n;bool must_include=false;
    // u=vertexChoose();
    // if(verifyQC()) 
    if(C_end <= LB) goto REC;
    if (P_end > LB && ( (double)P_end*(P_end-1)-2*MEInP ) >= gamma* (double)P_end*(P_end-1) ) {
        if(verify2hop(P_end)){
            store(P_end);
            printf("P_end: %d, MEInP: %d\n", P_end, MEInP);
            if(verifyQC()==false){
                printf("tree cnt: %d\n", treeCnt);
            }
            assert(verifyQC());
        }
    }
    if (( (double)(C_end)*(C_end-1)-2*MEInG ) >= gamma*( (double)(C_end)*(C_end-1) )) {
        prune1++;
        // printf("G: %d, MEInG: %d, %.3f\n",C_end, MEInG, gamma);
        if(verify2hop(C_end)){
            store(C_end); 
            assert(verifyQC());
            goto REC;
        }

    }
    if(sortBound()<=LB){
        ub_prune++;
        goto REC;
    }
    //if G[C] is a quasi clique
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return  
    // printf("enter break point\n");
    // u=vertexChoose();
    treeCnt++;
    
    for (int i = P_end; i < C_end; i++){
        int v=PC[i];
        if(neiInG[v]>=C_end-2){
            double edge_sum=(double)P_end*(P_end-1)-2*MEInP+2*neiInP[v];
            if(edge_sum >= gamma * (double)P_end*(P_end+1)){
                must_include=true;
                u=v;
                break;
            }
        }
    }
    if(u==n)u=PC[P_end];
    u=PC[P_end];
    CtoP(u,level);
    branch(level+1);// branch on adding the vertex u
    if(must_include&&false){
        PtoC(u,level);
        goto REC;
    }
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