#ifndef _MDSL_BB_
#define _MDSL_BB_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
#include "ColorPacker.h"
using namespace std;

long long treeCnt=0;
long long confTreeCnt=0;
long long ub_prune=0;
long long prune1=0;
bool usePrune1=true;
bool useUB = true;
bool useHeu = true;
bool useRed = true;

long long gapUBSum = 0;
long long totConflictNum = 0;
ui maxSubSz=0;
class QuasiClique_BB;
class kDefectiveClique_BB;
// class QuasiClique_BB2hop;
class QuasiClique_BB
{
private:
    ui n;
    ept m;
    ui maxDeg;
    ui minDeg;
    double gamma;
    bool subsearch;
    ui UB;
    ui P_end;ui C_end;
    ept MEInP, MEInG;
    ept *pstart;
    ui *edges;
    long long treeIdx;
    // bool *adjmtx;
    bool* matrix;

    ui *PC;
    ui *PC_rid;
    ui *neiInG;
    ui *neiInP;
    //graph coloring
    long long *colUseMtx;//to record if color is used
    ui *colorSz;// the size of each color bucket
    ui *colorVec;// the vertex in each color set
    int *colorLabel;
    bool *colvis;
    std::vector<std::vector<ui>> nonNeiB;
    std::vector<ui> weiB;
    std::vector<ui> weiPreSum;//the prefix sum of weight bucket
    std::vector<ui> QC;// current best solution

    ui* weight;
    ept *colExNum;
    vector<pair<ui,ui>> conflictPairs;
    vector<vector<bool>> MuEx;
    colorPacker *colPacker;
    // ept treeId;

    vector<vector<ui>> colorVecL; //the vertex in each color sets
public:
    ui LB;// current best size
    QuasiClique_BB();
    ~QuasiClique_BB();
    void load_graph(ui _n,ept *_pstart, ept *_pend, ui *_edges);
    void load_subgraph(double _gamma, ui _n, vector<pair<ui,ui>> &_vp, vector<ui> &_QC, ui _UB);
    void printInfo();
    void MQCSearch(double _gamma, ui _UB, std::vector<ui> &_QC);
    void MQCSearch2hop(vector<ui> &_QC);
    ui sortBound();
    ui sortConfBound();
    ui sortConfBoundL();
    ui complexBound();
    ui complexBoundL();
    ui simpleBound();
    ui sortBoundL();
    bool verifyQC();
    bool verify2hop(ui _end);
    void branchSubG(ui level);
    void branch(ui level);
    void branch2(ui level);
    void CtoP(ui u, ui level);
    void PtoC(ui u, ui level);
    void PtoX(ui u, ui level);
    void XtoC(ui u, ui level);
    // bool verifyQC();
    void store(ui newLB);
    bool is2hop(ui u, ui v);
    bool is2D(ui newLB);
    // bool is2hopAdj(ui u, ui v);
    void swapID(ui i, ui j);
    bool inP(ui u);
    bool inC(ui u);
    bool inX(ui u);
    bool &isAdj(ui u, ui v);
};

QuasiClique_BB::QuasiClique_BB(){
    n=0; m=0; gamma=0.0;
    maxDeg=0, minDeg=0;
    MEInP=MEInG=0;
    pstart=NULL;
    edges=NULL;
    weight = NULL;
    colExNum = NULL;
    matrix=NULL;
    PC=PC_rid=NULL;
    neiInG=neiInP=NULL;
    colUseMtx=NULL, colorSz=NULL, colorVec=NULL;
    colvis=NULL;
    colorLabel=NULL;
    P_end=C_end=0;
    QC.clear();
    nonNeiB.clear(), weiB.clear();
    weiPreSum.clear();
    LB=0; UB=0;
    treeIdx = 0;
    
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
    if(colvis!=NULL){
        delete[] colvis;
        colvis=NULL;
    }
    if(colorLabel!=NULL){
        delete[] colorLabel;
        colorLabel=NULL;
    }
    if(matrix!=NULL){
        delete[] matrix;
        matrix=NULL;
    }
    if(weight!=NULL){
        delete[] weight;
        weight = NULL;
    }
    if(colExNum!=NULL){
        delete[] colExNum;
        colExNum = NULL;
    }
    if(!conflictPairs.empty()) conflictPairs.clear();
    if(!MuEx.empty()) MuEx.clear();
    if(!colorVecL.empty()) colorVecL.clear();
}

void QuasiClique_BB::load_graph(ui _n,ept *_pstart, ept *_pend, ui *_edges){
    n=_n;
    C_end=n;
    //m is initially zero
    minDeg=n;
    for (ui i = 0; i < n; i++) m+=_pend[i]-_pstart[i];
    assert(pstart==NULL);
    pstart=new ept[n+1]; edges=new ui[m];
    neiInP=new ui[n];neiInG=new ui[n];
    PC=new ui[n]; PC_rid=new ui[n];
    colorSz=new ui[n]; 
    colorVec=new ui[n];
    colvis=new bool[n];
    colorLabel=new int[n];
    colExNum = new ept[n];
    weight =  new ui[n];
    m=0;
    memset(colvis, false, n*sizeof(bool));
    memset(colExNum, 0, sizeof(ept)*n);
    fill(colorLabel, colorLabel+n, -1);
    memset(neiInP,0,n*sizeof(ui));
    memset(weight, 0, n*sizeof(ui));
    for (ui i = 0; i < n; i++){
        PC[i]=PC_rid[i]=i;
        pstart[i]=m;
        neiInG[i]=_pend[i]-_pstart[i];
        maxDeg=max(maxDeg,neiInG[i]); minDeg=min(minDeg, neiInG[i]);
        for(ept j = _pstart[i];j < _pend[i];j ++) edges[m ++] = _edges[j];
    }
    pstart[n]=m;
    //renew the missing edges in G
    long long meing=(long long)n*(n-1)/2-m/2;//the number of missing edges in G
    printf("meing:%lld\n",meing);
    MEInG=meing;
   
    printf("load graph of size n=%u, m=%u (undirected), density=%.5lf, max degree=%d\n", n, m/2, double(m)/n/(n-1), maxDeg);
    m/=2;
}

void QuasiClique_BB::load_subgraph(double _gamma, ui _n, vector<pair<ui,ui>> &_vp, vector<ui> &_QC, ui _UB){
    bool onlyUB=true;
    gamma=_gamma;
    n=_n;
    maxSubSz=max(maxSubSz, n);
    UB=_UB;
    minDeg=n;
    C_end=n;
    treeIdx=0;
    m=_vp.size();
    //initialize the current best solution
    weight = new ui[n];
    colExNum = new ept[n];
    colPacker=new colorPacker(n);
    memset(weight, 0, sizeof(ui)*n);
    memset(colExNum, 0, sizeof(ept)*n);

    LB=_QC.size();
    matrix=new bool[n*n];
    PC=new ui[n];
    PC_rid=new ui[n];
    pstart=new ept[n+1];
    edges=new ui[2*m];
    neiInP= new ui[n];
    neiInG= new ui[n];
    colUseMtx=new long long[n*n]; colorSz=new ui[n]; 
    colorVec=new ui[n*n];
    memset(matrix, false, (n*n)*sizeof(bool));
    memset(neiInP, 0, n*sizeof(ui));
    if(onlyUB) fill(colUseMtx, colUseMtx+n*n, -1);
    else memset(colUseMtx, 0, n*n*sizeof(long long));
    MuEx.resize(n);
    for (int i = 0; i < n; i++) MuEx[i].resize(n, false);
    for(auto pr:_vp){
        ui u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
    }
    ept idx=0; 
    //construct the subgraph of pstart and edges
    for (ui i = 0; i < n; i++){
        pstart[i]=idx;
        for (ui j = 0; j < n; j++){
            if(isAdj(i,j)) edges[idx++]=j;
        }
        neiInG[i]=idx-pstart[i];
        maxDeg=max(maxDeg, neiInG[i]), minDeg=min(minDeg, neiInG[i]);
        PC[i]=i; PC_rid[i]=i;
    }
    pstart[n]=idx;
    m=idx/2;
    MEInG=n*(n-1)/2-m;
    // treeId=0;
}
void QuasiClique_BB::printInfo(){
    printf("vertex num: %d, edge num: %d\n",n,m);
    printf("max degree: %d, min degree: %d\n", maxDeg,minDeg);
}
void QuasiClique_BB::MQCSearch(double _gamma, ui _UB, std::vector<ui> &_QC){
    gamma=_gamma;
    UB=_UB;
    LB=_QC.size();
    ui u=PC[0];
    CtoP(u, 0);
    branch2(1);
    PtoX(u,0);
    branch2(1);
    XtoC(u, 0);
    if(LB>_QC.size()){
        //renew the best solution
        _QC.clear();
        for (ui i = 0; i < LB; i++) _QC.push_back(QC[i]);
    }
}

bool QuasiClique_BB::verifyQC(){
    bool flag=false;
    ui m_qc=0, n_qc=this->QC.size();
    double density=0.0;
    if(n_qc==0){
        printf("trivial gamma quasi clique\n");
        return true;
    }
    for (ui i = 0; i < QC.size(); i++){
        for (ui j = i+1; j < QC.size(); j++){
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
bool QuasiClique_BB::verify2hop(ui _end){
    bool flag=true, curflag=false;
    ui u_0=PC[0];
    vector<ui> nonNei_u0;
    for (ui i = 1; i < _end; i++){
        ui v=PC[i];
        if(!isAdj(u_0, PC[i])) nonNei_u0.push_back(v);
    }
    for (ui i = 0; i < _end; i++){
        ui u=PC[i];
        for (ui j = 0; j < nonNei_u0.size(); j++){
            ui v=nonNei_u0[j];
            if(u==v || isAdj(u,v)) continue;
            curflag=false;
            for (ui k = 0; k < _end; k++){
                ui w=PC[k];
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
void QuasiClique_BB::MQCSearch2hop(vector<ui> &_QC){
    Timer t;
    ui u=PC[0];
    CtoP(u, 0);
    branchSubG(1);
    PtoC(u,0);
    if(LB>_QC.size()){
        _QC.resize(LB);
        for (ui i = 0; i < QC.size(); i++) _QC[i]=QC[i];
    }
}
ui QuasiClique_BB::simpleBound(){
     // printf("enter MQC sort bound\n");
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1);
    weiB.resize(C_end);
    weiPreSum.resize(C_end-P_end);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		// if(P_end-neiInP[i]>k) continue;
		nonNeiB[ui(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    ui col=0;
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            col++;
        }
    }
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }
    //4. calculating the upper bound
    ui i = C_end;
    while (i>P_end && i*(i-1)/2 < gamma*i*(i-1)/2.0+MEInP+weiPreSum[i-P_end-1]){
        i--;
    }
    UB=i;
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear();
    return UB;
}
ui QuasiClique_BB::sortBound(){
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1);
    weiB.resize(C_end);
    weiPreSum.resize(C_end-P_end);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		nonNeiB[ui(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }
    //4. calculating the upper bound
    ui i = C_end;
    while (i>P_end && i*(i-1)/2 < gamma*i*(i-1)/2.0+MEInP+weiPreSum[i-P_end-1]){
        i--;
    }
    UB=i;
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear();
    return UB;
}
ui QuasiClique_BB::sortConfBoundL(){
    UB=n; ui maxWei=0; ui CMax = 0;//the max num put into the optimal solution
    ui maxColSz = 0;
    ui S_UB = 0, UB_gap = 0;
    ui colNum=0; int maxCol=0; // record the max color label and the color num
    ui newCsz = 0; // the new size of the C set
    nonNeiB.resize(P_end+1); for(auto& b:nonNeiB) b.clear(); // the max num of non nei in P is P_end+1
    weiB.resize(C_end); weiB.assign(C_end,0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end-P_end,0);
    
    colorSz[0]=0; // the size of color 0 is initialized as 0
    vector<ui> colvec0;
    colorVecL.push_back(colvec0);//push back the color 0 vector

    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            int col=0;
            //to record which color is used in the neighbors
            for (ept j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        int lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            
            while(colvis[col]) col++; //to find the minimum color not use
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
                vector<ui> colvec;
                colorVecL.push_back(colvec);//push back the color 0 vector
			}
            colorVecL[col].push_back(u);
            colorLabel[u]=col; colorSz[col]++; 
            maxColSz = max(maxColSz, colorSz[col]);
            ui wei=i+colorSz[col]-1;
            // if(u<0 || u>=n)
            weight[u]=wei;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
    }
    
    // 3. get the sort bound
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            vIdx++; weiB[w]--;
        }
    }

    //4. calculating the upper bound
    int ub = C_end;
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]) ub--;
    S_UB=ub; UB=S_UB;
    if(S_UB <= LB || maxColSz<= 0.5*C_end) goto RET4;
    
    weiB.assign(C_end, 0); //reset the weight bucket
    weiPreSum.assign(C_end-P_end,0);
    confTreeCnt++;
    for (int col = 0; col <= maxCol; col++){
        ui rem = colorSz[col];
        vector<ui> R;
        vector<vector<ui>> RVec;
        list<ui> colVecCpy(colorVecL[col].begin(), colorVecL[col].end());
        while(rem > 0){
            auto it = colVecCpy.begin();
            while (it != colVecCpy.end()){
                bool allconflict = true;  //to check if v in Pi_i is conflict to all other vertices in 
                ui v = *it;
                for(auto u: R){
                    if(is2hop(u, v)){
                        allconflict = false;
                        break;
                    }
                }
                if (allconflict){
                    R.push_back(*it);
                    it  = colVecCpy.erase(it);
                    rem--;
                }else ++it;
            }
            RVec.emplace_back(R);
            ui wei = P_end - neiInP[R[0]] + RVec.size() - 1;
            weiB[wei]++; newCsz++;
            R.clear();
        }
        RVec.clear(); colVecCpy.clear();
    }

    vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            vIdx++; weiB[w]--;
        }
    }
    //4. calculating the upper bound
    ub = P_end+newCsz;
    if(newCsz==0){
        printf("error in sort conflict bound\n");
        exit(0);
    }
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]) ub--;
    UB=ub;
    RET4:
    //recover the color label
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    }
    colorVecL.clear(); nonNeiB.clear(); weiB.clear();
    return UB;
}
ui QuasiClique_BB::sortConfBound(){
    ui UB=0; ui maxWei=0;
    // ui MeUB = floor((1.0-gamma)*C_end*(C_end-1)/2)-MEInP;
    nonNeiB.resize(P_end+1); for(auto &b: nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end-P_end, 0);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    ept totConlictPairs=0;
    ui newCsz = 0; // the new size of the C set
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		// if(P_end-neiInP[i]>k) continue;
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }

    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ept j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorVec[n*col+colorSz[col]] = u;//error
            //weight
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weight[u] = wei; weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }
    // ui kCur = K - MEInP, CMax = 0;
    ui S_UB = 0, UB_gap = 0;
    
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }

    //4. calculating the upper bound
    ui ub = C_end;
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]){
        ub--;
    }
    S_UB=ub;
    UB=S_UB;
    if(S_UB <= LB) goto RET2;

    // 3. get the conflict relationships in each color set
    if(true){
        for (int col = 0; col <= maxCol; col++){
            colExNum[col] = 0;
            for (ui i = 0; i < colorSz[col]; i++){
                int u = colorVec[n*col+i];
                for (ui j = i+1; j < colorSz[col]; j++){
                    int v = colorVec[n*col+j];
                    if(!is2hop(u, v)){
                        MuEx[u][v]=true;
                        MuEx[v][u]=true;
                        conflictPairs.emplace_back(make_pair(u, v));
                        totConlictPairs++; // record the num of conflict pairs
                        if(col<0 || col >= n) {
                            printf("error\n");
                            exit(0);
                        }
                        colExNum[col]++;
                    }
                }
            }
        }
    }
    
    //4. get the vertices in each colorset
    if(totConlictPairs>0){
        maxWei=0;
        weiB.assign(C_end,0); weiPreSum.assign(C_end-P_end, 0);
        for (int col = 0; col <= maxCol; col++){
            if(colExNum[col]>0){
                ui* colvec=colorVec+n*col;
				colPacker->set(MuEx, colvec, colorSz[col],C_end-1);
				colPacker->coloring(colvec,neiInP, P_end);
                newCsz+=colPacker->vChose.size();
				colPacker->getNewWeights(colvec, colorSz[col], neiInP, P_end, n);
                maxWei = max(maxWei, colPacker->maxWeight);
                for(ui wei = 0; wei <= colPacker->maxWeight; wei++){
                    weiB[wei]+=colPacker->newWeiB[wei];
                }
				colPacker->reset();
            }else{
                for (int i = 0; i < colorSz[col]; i++){
                    ui u = colorVec[n*col+i];
                    if(weight[u]<=C_end-1){
                        ui wei = weight[u];
                        weiB[wei]++;
                        maxWei = max(maxWei, wei);
                    }
                }
                newCsz+=colorSz[col];
            }
        }
        //5. get the upper bound
        vIdx=0;
        for (ui w = 0; w <=  maxWei; w++){
            while (weiB[w]!=0){
                if(vIdx==0) weiPreSum[vIdx]=w;
                else{
                    weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
                }
                vIdx++;
                weiB[w]--;
            }
        }
        if(vIdx!=newCsz) {
            printf("bound calculation err!\n");
            exit(0);
        }
        //4. calculating the upper bound
        ub = P_end+newCsz;
        if(newCsz==0){
            printf("error in sort conflict bound\n");
            exit(0);
        }
        while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]){
            ub--;
        }
        UB=ub;
    }
    if (true){
        if(UB>S_UB){
            printf("gap alg error!---------sort bound: %u, complex bound: %u\n", S_UB, UB);
            printf("colNum=%u, Psz=%u, Csz=%u\n",maxCol+1, P_end, C_end-P_end);
            for (ui c = 0; c <= maxCol; c++){
                printf("color: %u\n", c);
                for (ui i = 0; i < colorSz[c]; i++){
                    ui u = colorVec[n*c+i];
                    printf("%u: w=%u, w_s=%u\n", colorVec[n*c+i], P_end-neiInP[u], weight[u]);
                }
                
            }
            for (ui w = 0; w <= C_end-1; w++) printf("%u, ", weiB[w]);
            printf("\n");
            exit(0);
        }
        UB_gap  = S_UB - UB;
        if(UB_gap>n){
            printf("ub gap err, s_ub=%u, ub=%d, ub_gap=%u\n", S_UB, UB, UB_gap);
            exit(0);
        }
        gapUBSum += UB_gap;
    }
RET2: 
    nonNeiB.clear(); weiB.clear();weiPreSum.clear();
    if(!conflictPairs.empty()){
        for(auto p:conflictPairs){
            int u = p.first, v = p.second;
            MuEx[u][v]=false; MuEx[v][u]=false;
        }
        conflictPairs.clear();
    }
    return UB;
}
ui QuasiClique_BB::complexBoundL(){
    ui UB=0; ui maxWei=0; ui CMax = 0;//the max num put into the optimal solution
    // ui kCur=K-MEInP;
    nonNeiB.resize(P_end+1); for(auto& b:nonNeiB) b.clear(); // the max num of non nei in P is P_end+1
    weiB.resize(C_end); weiB.assign(C_end,0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end-P_end,0);
    ui colNum=0, maxCol=0; // record the max color label and the color num
    colorSz[0]=0; // the size of color 0 is initialized as 0
    vector<ui> colvec0;
    colorVecL.push_back(colvec0);//push back the color 0 vector
    ept totConlictPairs=0; // the total num of conflict pairs
    ui newCsz = 0; // the new size of the C set
    vector<pair<ui, ui>> conflictPairLoc; //to record the location of the start index and the end index in the confilctPairs vec
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		// if(P_end-neiInP[i]>k) continue;
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            //to record which color is used in the neighbors
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        ui lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            
            while(colvis[col]) col++; //to find the minimum color not use
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
                vector<ui> colvec;
                colorVecL.push_back(colvec);//push back the color 0 vector
			}
            colorVecL[col].push_back(u);
            colorLabel[u]=col; colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weight[u]=wei;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
        
    }
    ui S_UB = 0, UB_gap = 0;
    //3. get the sort bound
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }

    //4. calculating the upper bound
    ui ub = C_end;
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]){
        ub--;
    }
    S_UB=ub;
    UB=S_UB;
    //4. get the conflict relationships in each color set
    if(true){
        for (int col = 0; col <= maxCol; col++){
            colExNum[col] = 0;
            ui st = (ui)conflictPairs.size(), end = 0;
            for (int i = 0; i < colorSz[col]; i++){
                ui u = colorVecL[col][i];;
                for (int j = i+1; j < colorSz[col]; j++){
                    ui v = colorVecL[col][j];;
                    if(!is2hop(u, v)){
                        conflictPairs.emplace_back(make_pair(u, v));
                        totConlictPairs++; // record the num of conflict pairs
                        if(col<0 || col >= n) {
                            printf("error\n");
                            exit(0);
                        }
                        colExNum[col]++;
                    }
                }
            }
            end = (ui)conflictPairs.size();
            conflictPairLoc.emplace_back(make_pair(st, end));
        }
    }
    //5. get the vertices in each colorset
    if(totConlictPairs>0){
        weiB.assign(C_end,0); weiPreSum.assign(C_end-P_end, 0);
        maxWei = 0;
        for (int col = 0; col <= maxCol; col++){
            if(colExNum[col]>0){
                // ui* colvec=colorVec+n*col;
                pair<ui, ui> loc = conflictPairLoc[col];
				colPacker->set(conflictPairs, loc.first, loc.second, colorVecL[col],colorSz[col],C_end-1);
				colPacker->coloring(colorVecL[col],neiInP, P_end);
                newCsz+=colPacker->vChose.size();
				colPacker->getNewWeights(colorVecL[col], colorSz[col], neiInP, P_end, n);
                maxWei = max(maxWei, colPacker->maxWeight);
                for(ui wei = 0; wei <= colPacker->maxWeight; wei++){
                    weiB[wei]+=colPacker->newWeiB[wei];
                }
				colPacker->resetL();
            }else{
                for (int i = 0; i < colorSz[col]; i++){
                    ui u = colorVecL[col][i];
                    if(weight[u]<=C_end-1){
                        ui wei = weight[u];
                        weiB[wei]++;
                        maxWei = max(maxWei, wei);
                    }
                }
                newCsz+=colorSz[col];
            }
        }
        //5.2. get the upper bound
        //5. get the upper bound
        vIdx=0;
        for (ui w = 0; w <=  maxWei; w++){
            while (weiB[w]!=0){
                if(vIdx==0) weiPreSum[vIdx]=w;
                else{
                    weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
                }
                vIdx++;
                weiB[w]--;
            }
        }
        //4. calculating the upper bound
        ub = P_end+newCsz;
        if(newCsz==0){
            printf("error in sort conflict bound\n");
            exit(0);
        }
        while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]){
            ub--;
        }
        UB=ub;
    }

    RET4:
    //recover the color label
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    }
    colorVecL.clear(); nonNeiB.clear(); weiB.clear();
    if(!conflictPairs.empty()) conflictPairs.clear();
    if(!conflictPairLoc.empty()) conflictPairLoc.clear();
    return UB;
}
ui QuasiClique_BB::complexBound(){
    // printf("enter the quasi complex bound\n");
    ui UB=0; ui maxWei=0;
    ui colNum=0, maxCol=0;
    ui S_UB = 0, UB_gap = 0;
    ui maxColSz = 0;
    nonNeiB.resize(P_end+1); for(auto &b: nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end-P_end, 0);
    colorSz[0]=0;
    ept totConlictPairs=0;
    ui newCsz = 0; // the new size of the C se
    
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }

    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ept j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorVec[n*col+colorSz[col]] = u; colorSz[col]++;
            maxColSz = max(maxColSz, colorSz[col]);
            ui wei=i+colorSz[col]-1; weight[u] = wei; weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            vIdx++; weiB[w]--;
        }
    }

    //4. calculating the upper bound
    ui ub = C_end;
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]) ub--;
    S_UB=ub; UB=S_UB;
    if(S_UB <= LB || maxColSz <= 0.5*C_end) goto RET2;
    

    weiB.assign(C_end, 0); weiPreSum.assign(C_end-P_end,0);
    confTreeCnt++;
    for (ui col = 0; col <= maxCol; col++){   
        ui rem = colorSz[col]; // the remain num of \Pi_i
        vector<ui> R; vector<vector<ui>> RVec;
        ui* colvec = colorVec+col*n;
        list<ui> colVecCpy(colvec, colvec+colorSz[col]);
        while(rem > 0){
            auto it = colVecCpy.begin();
            while(it!=colVecCpy.end()){
                bool allconflict = true;  //to check if v in Pi_i is conflict to all other vertices in 
                ui v = *it;
                for(auto u: R){
                    if(is2hop(u, v)){
                        allconflict = false;
                        break;
                    }
                }
                if(allconflict){
                    R.push_back(*it);
                    it = colVecCpy.erase(it);
                    rem--;// delete the vertices in Pi_i
                }else ++it;
            }
            RVec.emplace_back(R);
            ui wei = P_end - neiInP[R[0]] + RVec.size() - 1;
            weiB[wei]++; newCsz++;
            R.clear();
        }
        RVec.clear(); colVecCpy.clear();
    }

    vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            vIdx++; weiB[w]--;
        }
    }
    //4. calculating the upper bound
    ub = P_end+newCsz;
    while (ub>P_end && ub*(ub-1)/2 < gamma*ub*(ub-1)/2.0+MEInP+weiPreSum[ub-P_end-1]) ub--;
    UB=ub;
    if (true){
        if(UB>S_UB){
            printf("gap alg error!---------sort bound: %u, complex bound: %u\n", S_UB, UB);
            printf("colNum=%u, Psz=%u, Csz=%u\n",maxCol+1, P_end, C_end-P_end);
            for (ui c = 0; c <= maxCol; c++){
                printf("color: %u\n", c);
                for (ui i = 0; i < colorSz[c]; i++){
                    ui u = colorVec[n*c+i];
                    printf("%u: w=%u, w_s=%u\n", colorVec[n*c+i], P_end-neiInP[u], weight[u]);
                }
                
            }
            for (ui w = 0; w <= C_end-1; w++) printf("%u, ", weiB[w]);
            printf("\n");
            exit(0);
        }
        UB_gap  = S_UB - UB;
        if(UB_gap>n){
            printf("ub gap err, s_ub=%u, ub=%d, ub_gap=%u\n", S_UB, UB, UB_gap);
            exit(0);
        }
        gapUBSum += UB_gap;
    }
RET2: 
    nonNeiB.clear(); weiB.clear();weiPreSum.clear();
    return UB;
}
ui QuasiClique_BB::sortBoundL(){
    UB=n; ui maxWei=0;
    nonNeiB.resize(P_end+1); for(auto &b: nonNeiB) b.clear();
    weiB.resize(C_end);  weiB.assign(C_end, 0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end-P_end, 0);
    ui colNum=0;
    int maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = ui(P_end-neiInP[u]);
		nonNeiB[nonNei].push_back(u);
    }
    // 2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            int col=0;
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        int lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            
            while(colvis[col]) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            colorLabel[u]=col;
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
        
    }
    // 3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }
    int i = C_end;
    while (i>P_end && i*(i-1)/2 < gamma*i*(i-1)/2.0+MEInP+weiPreSum[i-P_end-1]){
        i--;
    }
    UB=i;
    printf("complete ub cal: treeIdx=%lld, ub=%d\n", treeIdx, UB);
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear(); 
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    } 
    return UB;
}
void QuasiClique_BB::branch(ui level){
    assert(level<=n+1);
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if (P_end > LB && ( (double)P_end*(P_end-1)-2*MEInP ) >= gamma* (double)P_end*(P_end-1) ) store(P_end);
    if( (double(m)/C_end/(C_end-1))*2.0>= gamma){
        if(C_end<=1) printf("Error!\n"),exit(0);
        prune1++; store(C_end); goto REC;
    }
    if(sortBoundL()<=LB&&useUB){ub_prune++; goto REC;}
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return 
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                double edge_sum=(double)P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= gamma * (double)P_end*(P_end+1)){
                    must_include=true;
                    u=v; break;
                }
            }
        }
    }
    if(u==n) u=PC[P_end];

    CtoP(u,level);
    branch(level+1);// branch on adding the vertex u
    if(must_include){
        prune1++;
        PtoC(u,level);
        goto REC;
    }
    PtoX(u,level);
    branch(level+1);// branch on deleting the vertex u
    XtoC(u, level);
REC:
    return;
}

void QuasiClique_BB::branch2(ui level){
    
    assert(level<=n+1);
    treeIdx++;
    printf("enter branch2, treeIdx=%ld\n", treeIdx);
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if (P_end > LB && ( (double)P_end*(P_end-1)-2*MEInP ) >= gamma* (double)P_end*(P_end-1) && is2D(P_end)) store(P_end);
    if( (double(m)/C_end/(C_end-1))*2.0>= gamma && is2D(C_end)){
        if(C_end<=1) printf("Error!\n"),exit(0);
        printf("G: %d, MEInG: %d, %.3f\n",C_end, MEInG, gamma);
        prune1++; store(C_end); goto REC;
    }
    printf("enter the sortBound, treeId = %lld\n", treeIdx);
    if(sortBoundL()<=LB&&useUB){ub_prune++; 
        goto REC;}
    printf("out the sortBound, treeId = %lld\n", treeIdx);
    
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return 
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                double edge_sum=(double)P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= gamma * (double)P_end*(P_end+1)){
                    must_include=true;
                    u=v; break;
                }
            }
        }
    }
    if(u==n) u=PC[P_end];
    
    CtoP(u,level);
    branch2(level+1);// branch on adding the vertex u
    if(must_include){
        if(verify2hop(P_end)){
            prune1++;
            PtoC(u,level);
            goto REC;
        }
    }
    PtoX(u,level);
    branch2(level+1);// branch on deleting the vertex u
    XtoC(u, level);
REC:
    printf("recover the branch: %ld\n", treeIdx);
    return;
}
void QuasiClique_BB::branchSubG(ui level){
    treeIdx++;
    assert(level<=n+1);
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if (P_end > LB && ( (double)P_end*(P_end-1)-2*MEInP ) >= gamma* (double)P_end*(P_end-1) ) {
         if(verify2hop(P_end)){
            store(P_end);
            assert(verifyQC());
        }
    }
    if (( (double)(C_end)*(C_end-1)-2*MEInG ) >= gamma*( (double)(C_end)*(C_end-1) )) {
        printf("G: %d, MEInG: %d, %.3f\n",C_end, MEInG, gamma);
        if(verify2hop(C_end)){
            prune1++; store(C_end); 
            assert(verifyQC());
            goto REC;
        }
    }
    
    if(complexBound()<=LB&&useUB){ub_prune++;goto REC;}
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return  
    // printf("enter break point\n");
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                double edge_sum=(double)P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= gamma * (double)P_end*(P_end+1)){
                    must_include=true;
                    u=v;
                    break;
                }
            }
        }
    }
    
    if(u==n) u=PC[P_end];

    CtoP(u,level);
    branchSubG(level+1);// branch on adding the vertex u
    if(must_include){
        if(verify2hop(P_end)){
            prune1++;
            PtoC(u,level);
            goto REC;
        }
    }
    PtoX(u,level);
    branchSubG(level+1);// branch on deleting the vertex u
    XtoC(u, level);

REC:

    return;
}
void QuasiClique_BB::PtoC(ui u, ui level){
    assert(inP(u));
    swapID(PC_rid[u],--P_end);
    ui nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;//update the missing edges in P
    //update the neiInP of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]--;
    }
    
}
void QuasiClique_BB::CtoP(ui u, ui level){
    assert(inC(u));
    swapID(PC_rid[u], P_end++);
    ui nonNeiP=P_end-1-neiInP[u];
    MEInP+=nonNeiP;
    //update the neiInP of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]++;
    }
}
void QuasiClique_BB::PtoX(ui u, ui level){
    assert(inP(u));
    //move u from P to C
    swapID(PC_rid[u],--P_end);
    ui nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;
    //move u from C to X
    swapID(P_end, --C_end);
    ui nonNeiG=C_end-neiInG[u];
    m-=neiInG[u];
    MEInG-=nonNeiG;
    //update the neiInP and neiInG of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]--;
        neiInG[v]--;
    }
}
void QuasiClique_BB::XtoC(ui u, ui level){
    assert(inX(u));
    //move u from X to C
    swapID(PC_rid[u],C_end++);
    ui nonNeiG=C_end-1-neiInG[u];
    m+=neiInG[u];
    MEInG+=nonNeiG;
     //update the neiInG of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInG[v]++;
    }
}
void QuasiClique_BB::store(ui newLB){
    QC.resize(LB=newLB);
    for (ui i = 0; i < LB; i++) QC[i]=PC[i];
}

bool QuasiClique_BB::is2hop(ui u, ui v){
    // bool is2hop = false;
    if (neiInG[u]+neiInG[v]>C_end-2) return true;// if the sum of degree of u and v larger than n-2, then u adj v/ dist(u,v)==2
    unordered_set<ui> u_neis;
    for (ept i = pstart[u]; i < pstart[u+1]; i++){//find the neighbors of the vertex u
        ui u_nei = edges[i];
        if(v==u_nei) return true;
        if(PC_rid[u_nei] >= 0 && PC_rid[u_nei] < C_end) u_neis.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pstart[v+1]; i++){
        ui v_nei = edges[i];
        if(u_neis.find(v_nei) != u_neis.end()) return true;
    }
    return false;
}

bool QuasiClique_BB::is2D(ui newLB){
    bool is2hop=true;
    for (int i = 0; i < newLB; i++){
        int u = PC[i];
        unordered_set<int> neis;
        //get the neighbors of u in PC[0]~PC[newLB]
        for (int j = pstart[u]; j < pstart[u+1]; j++){
            int v = edges[j];
            if(PC_rid[v] >= 0 && PC_rid[v] < newLB) neis.insert(v);
        }
        for (int j = i+1; j < newLB; j++){
            int v = PC[j];
            is2hop = false;
            //check if v is in the neis of u
            if(neis.find(v)!=neis.end()){
                is2hop=true;
            }else{
                for (int k = pstart[v]; k < pstart[v+1]; k++){
                    int w = edges[k];
                    if(neis.find(w)!=neis.end()) {
                        is2hop = true;
                        break;
                    }
                }   
            }
            if(!is2hop) return false;
        }
        neis.clear();
    }
    return true;
}

void QuasiClique_BB::swapID(ui i,ui j){
    swap(PC[i],PC[j]);
    PC_rid[PC[i]] = i;
	PC_rid[PC[j]] = j;
}
bool QuasiClique_BB::inP(ui u){
    return (PC_rid[u]>=0 && PC_rid[u]<P_end);
}
bool QuasiClique_BB::inC(ui u){
    return (PC_rid[u]>=P_end && PC_rid[u] < C_end);
}
bool QuasiClique_BB::inX(ui u){
    return (PC_rid[u]>=C_end && PC_rid[u]<n);
}
bool& QuasiClique_BB::isAdj(ui u,ui v){
    return matrix[u*n+v];
}

class kDefectiveClique_BB
{
private:
    ui n;
    ept m;
    ui maxDeg;
    ui minDeg;
    ui K;
    ui UB;
    ui P_end;ui C_end;
    ept MEInP, MEInG;
    ept *pstart;
    ui *edges;
    ui* weight;
    long long treeIdx;
    // bool *adjmtx;
    bool* matrix;

    ui *PC;
    ui *PC_rid;
    ui *neiInG;
    ui *neiInP;
    //graph coloring
    long long *colUseMtx;//to record if color is used
    ui *colorSz;// the size of each color bucket
    ui *colorVec;// the vertex in each color set
    int *colorLabel;
    bool *colvis;
    ept *colExNum;
    std::vector<std::vector<ui>> nonNeiB;
    std::vector<ui> weiB;
    std::vector<ui> weiPreSum;//the prefix sum of weight bucket
    std::vector<ui> KDC;// current best solution
    vector<vector<bool>> MuEx;
    vector<pair<ui,ui>> conflictPairs;
    colorPacker *colPacker;

    //data structure used in complexL
    vector<vector<ui>> colorVecL; //the vertex in each color sets
    

    
public:
    ui LB;// current best size
    kDefectiveClique_BB();
    ~kDefectiveClique_BB();
    void load_graph(ui _n,ept *_pstart, ept *_pend, ui *_edges);
    void load_subgraph(ui _K, ui _n, vector<pair<ui,ui>> &_vp, vector<ui> &_KDC, ui _UB);
    void printInfo();
    void MKDCSearch(ui _K, ui _UB, std::vector<ui> &_KDC);
    void MKDCSearch2hop(vector<ui> &_KDC);
    ui sortBound();
    ui sortConfBoundL();
    ui sortConfBound();
    ui complexBound();
    ui complexBoundL();
    ui simpleBound();
    ui sortBoundL();
    bool verifyKDC();
    bool verify2hop(ui _end);
    ui vtxSelect();
    ui reduction(ui level);
    ui recover(ui pushCnt, ui level);
    void branchSubG(ui level);
    void branch(ui level);
    void branch2(ui level);
    void CtoP(ui u, ui level);
    void PtoC(ui u, ui level);
    void PtoX(ui u, ui level);
    void XtoC(ui u, ui level);
    void CtoX(ui u, ui level);
    // bool verifyKDC();
    void store(ui newLB);
    bool is2D(ui newLB);
    bool is2hop(ui u, ui v);
    bool is2hop2(ui u, ui v);
    void swapID(ui i, ui j);
    bool inP(ui u);
    bool inC(ui u);
    bool inX(ui u);
    bool &isAdj(ui u, ui v);
};

kDefectiveClique_BB::kDefectiveClique_BB(){
    n=0; m=0; K=0;
    maxDeg=0, minDeg=0;
    treeIdx = 0; //the tree Id in branch is initializad as 0
    MEInP=MEInG=0;
    pstart=NULL;
    edges=NULL;
    weight = NULL;
    colExNum = NULL;
    matrix=NULL;
    PC=PC_rid=NULL;
    neiInG=neiInP=NULL;
    colUseMtx=NULL, colorSz=NULL, colorVec=NULL;
    colvis=NULL;
    colorLabel=NULL;
    P_end=C_end=0;
    KDC.clear();
    nonNeiB.clear(), weiB.clear();
    weiPreSum.clear();
    MuEx.clear();
    LB=0; UB=0;
    
}

kDefectiveClique_BB::~kDefectiveClique_BB(){
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
    if(weight!=NULL){
        delete[] weight;
        weight = NULL;
    }
    if(neiInG!=NULL){
        delete[] neiInG;
        neiInG=NULL;
    }
    if(!KDC.empty()){
        KDC.clear();
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
    if(!MuEx.empty()) MuEx.clear();
    if(!conflictPairs.empty()) conflictPairs.clear();
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
    if(colvis!=NULL){
        delete[] colvis;
        colvis=NULL;
    }
    if(colorLabel!=NULL){
        delete[] colorLabel;
        colorLabel=NULL;
    }
    if(colExNum!=NULL){
        delete[] colExNum;
        colExNum = NULL;
    }
    if(matrix!=NULL){
        delete[] matrix;
        matrix=NULL;
    }
    if(!colorVecL.empty()) colorVecL.clear();
}

void kDefectiveClique_BB::load_graph(ui _n,ept *_pstart, ept *_pend, ui *_edges){
    n=_n;
    C_end=n;
    //m is initially zero
    minDeg=n;
    for (ui i = 0; i < n; i++) m+=_pend[i]-_pstart[i];
    assert(pstart==NULL);
    pstart=new ept[n+1]; edges=new ui[m];
    neiInP=new ui[n];neiInG=new ui[n];
    PC=new ui[n]; PC_rid=new ui[n];
    colorSz=new ui[n]; 
    weight = new ui[n];
    colvis=new bool[n];
    colorLabel=new int[n];
    colExNum = new ept[n];
    m=0;
    fill(colvis, colvis+n, false);
    memset(weight, 0, sizeof(ui)*n);
    memset(colExNum, 0, sizeof(ept)*n);
    fill(colorLabel, colorLabel+n, -1);
    memset(neiInP,0,n*sizeof(ui));
    for (ui i = 0; i < n; i++){
        PC[i]=PC_rid[i]=i;
        pstart[i]=m;
        neiInG[i]=_pend[i]-_pstart[i];
        maxDeg=max(maxDeg,neiInG[i]); minDeg=min(minDeg, neiInG[i]);
        for(ept j = _pstart[i];j < _pend[i];j ++) edges[m ++] = _edges[j];
    }
    pstart[n]=m;
    //renew the missing edges in G
    long long meing=(long long)n*(n-1)/2-m/2;//the number of missing edges in G
    printf("meing:%lld\n",meing);
    MEInG=meing;
   
    printf("load graph of size n=%u, m=%u (undirected), density=%.5lf, max degree=%d\n", n, m/2, double(m)/n/(n-1), maxDeg);
    m/=2;
}

void kDefectiveClique_BB::load_subgraph(ui _K, ui _n, vector<pair<ui,ui>> &_vp, vector<ui> &_KDC, ui _UB){
    bool onlyUB=true;
    K=_K;
    n=_n;
    maxSubSz=max(maxSubSz, n);
    UB=_UB;
    minDeg=n;
    C_end=n;
    treeIdx=0;
    m=_vp.size();
    //initialize the current best solution
    LB=_KDC.size();
    matrix=new bool[n*n];
    PC=new ui[n];
    PC_rid=new ui[n];
    pstart=new ept[n+1];
    edges=new ui[2*m];
    neiInP= new ui[n];
    neiInG= new ui[n];
    weight = new ui[n];
    colExNum = new ept[n];
    colUseMtx=new long long[n*n]; colorSz=new ui[n]; 
    colorVec=new ui[n*n];
    colPacker=new colorPacker(n);
    memset(matrix, false, (n*n)*sizeof(bool));
    memset(neiInP, 0, n*sizeof(ui));
    memset(weight, 0, sizeof(ui)*n);
    memset(colExNum, 0, sizeof(ept)*n);
    memset(colorVec, 0, sizeof(ui)*n*n);
    if (onlyUB) fill(colUseMtx, colUseMtx+n*n, -1);
    else memset(colUseMtx, 0, n*n*sizeof(long long));
    MuEx.resize(n);
    for (int i = 0; i < n; i++) MuEx[i].resize(n, false);
    for(auto pr:_vp){
        ui u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
    }
    ept idx=0; 
    //construct the subgraph of pstart and edges
    for (ui i = 0; i < n; i++){
        pstart[i]=idx;
        for (ui j = 0; j < n; j++){
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
void kDefectiveClique_BB::printInfo(){
    printf("vertex num: %d, edge num: %d\n",n,m);
    printf("max degree: %d, min degree: %d\n", maxDeg,minDeg);
}
void kDefectiveClique_BB::MKDCSearch(ui _K, ui _UB, std::vector<ui> &_KDC){
    K=_K;
    UB=_UB;
    LB=_KDC.size();
    // subsearch=false;
    //use branch and bound to search
    printf("enter MKDC search 0114\n");
    ui u=PC[0];
    CtoP(u, 0);
    branch2(1);
    PtoX(u,0);
    branch2(1);
    XtoC(u, 0);
    if(LB>_KDC.size()){
        //renew the best solution
        _KDC.clear();
        for (ui i = 0; i < LB; i++) _KDC.push_back(KDC[i]);
    }
}

bool kDefectiveClique_BB::verifyKDC(){
    //need to rephrase
    bool flag=false;
    ui m_qc=0, n_qc=this->KDC.size();
    ui MeNum=0;
    if(n_qc==0){
        printf("trivial K quasi clique\n");
        return true;
    }
    for (ui i = 0; i < KDC.size(); i++){
        for (ui j = i+1; j < KDC.size(); j++){
            if(isAdj(KDC[i], KDC[j])) m_qc++;
        }
    }
    if(m_qc >= n_qc * (n_qc-1)/2 - K) return true;
    MeNum=n_qc * (n_qc - 1)/2 - m_qc;
    printf("KDC vNum: %d, KDC eNum: %d, KDC MeNum: %.2f\n",n_qc, m_qc, MeNum);
    return false;
}
bool kDefectiveClique_BB::verify2hop(ui _end){
    if(_end > K+1) return true;
    bool flag=true, curflag=false;
    ui u_0=PC[0];
    vector<ui> nonNei_u0;
    for (ui i = 1; i < _end; i++){
        ui v=PC[i];
        if(!isAdj(u_0, PC[i])) nonNei_u0.push_back(v);
    }
    for (ui i = 0; i < _end; i++){
        ui u=PC[i];
        for (ui j = 0; j < nonNei_u0.size(); j++){
            ui v=nonNei_u0[j];
            if(u==v || isAdj(u,v)) continue;
            curflag=false;
            for (ui k = 0; k < _end; k++){
                ui w=PC[k];
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
ui kDefectiveClique_BB::vtxSelect(){
    ui u = PC[P_end];
    for (ui i = P_end+1; i < C_end; i++){
        ui v = PC[i];
        if(neiInP[u]<neiInP[v]) u = v;
    }
    return u;
}
ui kDefectiveClique_BB::reduction(ui level){
    ui pushCnt=0;
    for (ui i = P_end; i < C_end;){
        ui v = PC[i];
        if(P_end - neiInP[v] >K - MEInP) CtoX(v, level), pushCnt++;
        else i++;
    }
    return pushCnt;
}
ui kDefectiveClique_BB::recover(ui pushCnt, ui level){
    while (pushCnt--){
        ui v = PC[C_end];
        XtoC(v, level);
    }
    return 0;
}
void kDefectiveClique_BB::MKDCSearch2hop(vector<ui> &_KDC){
    Timer t;
    ui u=PC[0];
    // subsearch=true;
    CtoP(u, 0);
    branchSubG(1);
    PtoC(u,0);
    if(LB>_KDC.size()){
        _KDC.resize(LB);
        for (ui i = 0; i < KDC.size(); i++) _KDC[i]=KDC[i];
    }

}
ui kDefectiveClique_BB::simpleBound(){
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1);
    weiB.resize(C_end);
    weiPreSum.resize(C_end-P_end);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		nonNeiB[ui(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    ui col=0;
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            col++;
        }
        
    }
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        while (weiB[w]!=0){
            if(vIdx==0) weiPreSum[vIdx]=w;
            else{
                weiPreSum[vIdx]+=weiPreSum[vIdx-1]+w;
            }
            vIdx++;
            weiB[w]--;
        }
    }
    //4. calculating the upper bound
    ui i = C_end;
    while (i>P_end && K < MEInP+weiPreSum[i-P_end-1]){
        i--;
    }
    UB=i;
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear();
    return UB;
}

ui kDefectiveClique_BB::sortBound(){
    // printf("enter MKDC sort bound\n");
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1);
    weiB.resize(C_end);
    weiPreSum.resize(C_end-P_end);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                // if(!isAdj(u,nei)) continue;
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }
    ui kCur = K - MEInP; ui CMax = 0;
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            // printf("remain=%u\n", kCur-w*weiB[w]);
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
    }
    UB=P_end+CMax;
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear();
    return UB;
}
ui kDefectiveClique_BB::complexBoundL(){
    printf("enter the complex bound, treeId=%ld, kCur=%u, Psz=%u, Csz=%u\n", treeIdx,K-MEInP, P_end, C_end-P_end);
    ui UB=0; ui maxWei=0; ui CMax = 0;//the max num put into the optimal solution
    ui kCur=K-MEInP; // the current k value
    nonNeiB.resize(P_end+1); // the max num of non nei in P is P_end+1
    for(auto& b: nonNeiB) b.clear();
    weiB.resize(C_end);//the max weight <=k
    weiB.assign(C_end, 0);
    ui colNum=0, maxCol=0; // record the max color label and the color num
    colorSz[0]=0; // the size of color 0 is initialized as 0
    vector<ui> colvec0;
    colorVecL.push_back(colvec0);//push back the color 0 vector
    ept totConlictPairs=0; // the total num of conflict pairs
    ui newCsz = 0; // the new size of the C set
    vector<pair<ui, ui>> conflictPairLoc; //to record the location of the start index and the end index in the confilctPairs vec
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            //to record which color is used in the neighbors
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        ui lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            
            while(colvis[col]) col++; //to find the minimum color not use
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
                vector<ui> colvec;
                colorVecL.push_back(colvec);//push back the color 0 vector
			}
            colorVecL[col].push_back(u);
            colorLabel[u]=col; colorSz[col]++;
            ui wei=i+colorSz[col]-1; weight[u]=wei;
            maxWei=max(maxWei, wei);
            weiB[wei]++;
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
        
    }
    //3. get the sort bound
    ui S_UB = 0, UB_gap = 0;
    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
    }
    S_UB = CMax+P_end;
    UB=S_UB;
    if(S_UB <= LB) goto RET3;
    kCur = K - MEInP, CMax = 0;

    //4. get the conflict relationships in each color set
    if(true){
        for (int col = 0; col <= maxCol; col++){
            colExNum[col] = 0;
            ui st = (ui)conflictPairs.size(), end = 0;
            for (int i = 0; i < colorSz[col]; i++){
                ui u = colorVecL[col][i];
                for (int j = i+1; j < colorSz[col]; j++){
                    ui v = colorVecL[col][j];
                    if(!is2hop(u, v)){
                        conflictPairs.emplace_back(make_pair(u, v));
                        totConlictPairs++; // record the num of conflict pairs
                        if(col<0 || col >= n) {
                            printf("error\n");
                            exit(0);
                        }
                        colExNum[col]++;
                    }
                }
            }
            end = (ui)conflictPairs.size();
            conflictPairLoc.emplace_back(make_pair(st, end));
        }
    }
    //5. get the vertices in each colorset
    if(totConlictPairs>0){
        weiB.resize(kCur+1,0), maxWei = 0; weiB.assign(kCur+1, 0);
        for (int col = 0; col <= maxCol; col++){
            if(colExNum[col]>0){
                pair<ui, ui> loc = conflictPairLoc[col];
				colPacker->set(conflictPairs, loc.first, loc.second, colorVecL[col],colorSz[col],kCur);
				colPacker->coloring(colorVecL[col],neiInP, P_end);
                newCsz+=colPacker->vChose.size();
				colPacker->getNewWeights(colorVecL[col], colorSz[col], neiInP, P_end, kCur);
                maxWei = max(maxWei, colPacker->maxWeight);
                for(ui wei = 0; wei <= colPacker->maxWeight; wei++){
                    weiB[wei]+=colPacker->newWeiB[wei];
                }
				colPacker->resetL();
            }else{
                for (int i = 0; i < colorSz[col]; i++){
                    ui u = colorVecL[col][i];
                    if(weight[u]<=kCur){
                        ui wei = weight[u];
                        weiB[wei]++;
                        maxWei = max(maxWei, wei);
                        newCsz ++;
                    }
                }
            }
        }
        //5.2. get the upper bound
        printf("begin the sortConfBound calculation\n");
        for (ui w = 0; w <=  maxWei; w++){
            if(kCur>=w*weiB[w]) {
                CMax+=weiB[w], kCur=kCur-w*weiB[w];
            }
            else {
                CMax+=kCur/w;
                break;
            }
        }
        UB=P_end+CMax;
        printf("end the sortConfBound calculation\n");
    }
    if (true){
        if(S_UB< P_end+CMax){
            printf("gap alg error!---------sort bound: %u, complex bound: %u\n", S_UB, P_end+CMax);
            printf("colNum=%u, Psz=%u, Csz=%u\n",maxCol+1, P_end, C_end-P_end);
            for (ui c = 0; c <= maxCol; c++){
                printf("color: %u\n", c);
                for (ui i = 0; i < colorSz[c]; i++){
                    ui u = colorVecL[c][i];
                    printf("%u: w=%u, w_s=%u\n", u, P_end-neiInP[u], weight[u]);
                }
            }
            for (ui w = 0; w <= K-MEInP; w++){
                printf("%u, ", weiB[w]);
            }
            printf("\n");
            exit(0);
        }
        UB_gap  = S_UB - (P_end+CMax);
        gapUBSum += UB_gap;
    }
    RET3:
    // return P_end+CMax;
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    }
    colorVecL.clear(); nonNeiB.clear(); weiB.clear();
    if(!conflictPairs.empty()) conflictPairs.clear();
    if(!conflictPairLoc.empty()) conflictPairLoc.clear();
    return UB;
}
ui kDefectiveClique_BB::sortConfBound(){
    //1. initialize
    bool debug = false;//if enter the debug mode
    UB=0; 
    ui maxWei=0, kCur = K - MEInP;
    ui colNum=0, maxCol=0; 
    ui S_UB = 0, CMax = 0, UB_gap = 0;
    nonNeiB.resize(kCur+1); //the non nei no more than kCur
    for(auto& b: nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    colorSz[0]=0;

    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
        if(nonNei <= kCur) nonNeiB[nonNei].push_back(u);
    }

    //2. coloring the vertices in C
    for (ui i = 0; i <= kCur; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorVec[n*col+colorSz[col]] = u;
            //weight
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weight[u] = wei;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }

    //3. calculate the sort bound
    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
    }
    S_UB = CMax+P_end;
    UB = S_UB;
    if(S_UB <= LB) goto RETSCBOUND;

    //4. calculate the sortconfbound
    weiB.assign(C_end, 0); //reset the weight bucket
    confTreeCnt++;
    for (ui col = 0; col <= maxCol; col++){
        for (ui i = 0; i < colorSz[col]; i++){
            ui u = colorVec[n*col+i];
            for (ui j = i+1; j < colorSz[col]; j++){
                ui v = colorVec[n*col+j];
                if(is2hop(u, v)!= is2hop2(u, v)){
                    printf("conflict error, (u, v)=(%u, %u)\n", u, v);
                    if(is2hop(u, v)) printf("uv not conflict\n");
                    else printf("uv conflict\n");
                    if(is2hop2(u, v)) printf("uv not conflict2\n");
                    else printf("uv conflict2\n");
                    printf("u nei: ");
                    for(ept k = pstart[u]; k < pstart[u+1]; k++){
                        ui nei = edges[k];
                        if(PC_rid[nei]>=0 && PC_rid[nei] < C_end) printf("%u, ", nei);
                    }
                    printf("\n");
                    printf("v nei: ");
                    for(ept k = pstart[v]; k < pstart[v+1]; k++){
                        ui nei = edges[k];
                        if(PC_rid[nei]>=0 && PC_rid[nei] < C_end) printf("%u, ", nei);
                    }
                    printf("\n");
                    printf("uv cmNei: ");
                    for(ui i = 0; i < C_end; i++){
                        ui uu = PC[i];
                        if(isAdj(uu, v) && isAdj(uu, u)) printf("%u, ", uu);
                    }
                    printf("\n");
                    exit(0);
                }
                if(!is2hop(u, v)) totConflictNum++;
            }
        }
        
        ui rem = colorSz[col]; // the remain num of \Pi_i
        vector<ui> R;
        vector<vector<ui>> RVec;
        ui* colvec = colorVec+col*n;
        list<ui> colVecCpy(colvec, colvec+colorSz[col]);
        while(rem > 0){
            auto it = colVecCpy.begin();
            while(it!=colVecCpy.end()){
                bool allconflict = true;  //to check if v in Pi_i is conflict to all other vertices in 
                ui v = *it;
                for(auto u: R){
                    if(is2hop(u, v)){
                        allconflict = false;
                        break;
                    }
                }
                if(allconflict){
                    R.push_back(*it);
                    it = colVecCpy.erase(it);
                    rem--;// delete the vertices in Pi_i
                }else ++it;
            }
            RVec.emplace_back(R);
            ui wei = P_end - neiInP[R[0]] + RVec.size() - 1;
            weiB[wei]++;
            R.clear();
        }
        RVec.clear(); colVecCpy.clear();
    }
    kCur = K - MEInP, CMax = 0;
    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
    }
    UB = P_end + CMax;

    RETSCBOUND:
    UB_gap  = S_UB - UB; gapUBSum += UB_gap;
    nonNeiB.clear(); weiB.clear();
    return UB;
}
ui kDefectiveClique_BB::complexBound(){
    // printf("enter MKDC sort bound\n");
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1); 
    for(auto &b: nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    ept totConlictPairs=0; 
    ui newCsz = 0; // the new size of the C set
    //1. use bucket sorting to sort
    if(false)//verify the correctness
    {
        if(!conflictPairs.empty()) {
            printf("error: non empty conflict pairs\n");
            exit(0);
        }
        for(ui u=0; u<n; u++){
            for(ui v=0; v<n; v++){
                if(MuEx[u][v]){
                    printf("error : init muex fail\n");
                }
            }
        }
        for (ui i = 0; i < P_end+1; i++){
            if(!nonNeiB[i].empty()){
                printf("error nonNeiB\n");
                exit(0);
            }
        }
        for(ui i=0; i<C_end; i++){
            if(weiB[i]!=0){
                printf("error weiB\n");
                exit(0);
            }
        }
    }
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
        ui nonNei = P_end-neiInP[u];
		nonNeiB[nonNei].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            while(colUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                colUseMtx[nei*n+col]=treeIdx;
            }
            colorVec[n*col+colorSz[col]] = u;
            //weight
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weight[u] = wei;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
        }
    }
    ui kCur = K - MEInP, CMax = 0;
    ui S_UB = 0, UB_gap = 0;

    if(true){
        for (ui w = 0; w <=  maxWei; w++){
            if(kCur>=w*weiB[w]) {
                CMax+=weiB[w], kCur=kCur-w*weiB[w];
            }
            else {
                CMax+=kCur/w;
                break;
            }
        }
        S_UB = CMax+P_end;
        UB = S_UB;
        if(S_UB <= LB && false) goto RETCOMPLEX;
        kCur = K - MEInP, CMax = 0;
    }

    //3. get the conflict relationships in each color set
    if(true){
        for (int col = 0; col <= maxCol; col++){
            colExNum[col] = 0;
            for (int i = 0; i < colorSz[col]; i++){
                ui u = colorVec[n*col+i];
                for (int j = i+1; j < colorSz[col]; j++){
                    ui v = colorVec[n*col+j];
                    if(!is2hop(u, v)){
                        MuEx[u][v]=MuEx[v][u]=true;
                        conflictPairs.emplace_back(make_pair(u, v));
                        totConlictPairs++; // record the num of conflict pairs
                        if(col<0 || col >= n) {
                            printf("error\n");
                            exit(0);
                        }
                        colExNum[col]++;
                    }
                }
            }
        }
    }
    
    
    //4. get the vertices in each colorset
    if(totConlictPairs>0){
        weiB.resize(kCur+1), maxWei = 0; weiB.assign(kCur+1, 0);
        for (int col = 0; col <= maxCol; col++){
            if(colExNum[col]>0){
                ui* colvec=colorVec+n*col;
				colPacker->set(MuEx, colvec, colorSz[col],kCur);
				colPacker->coloring(colvec,neiInP, P_end);
                newCsz+=colPacker->vChose.size();
				colPacker->getNewWeights(colvec, colorSz[col], neiInP, P_end, kCur);
                maxWei = max(maxWei, colPacker->maxWeight);
                for(ui wei = 0; wei <= colPacker->maxWeight; wei++){
                    weiB[wei]+=colPacker->newWeiB[wei];
                }
				colPacker->reset();
            }else{
                for (int i = 0; i < colorSz[col]; i++){
                    ui u = colorVec[n*col+i];
                    if(weight[u]<=kCur){
                        ui wei = weight[u];
                        weiB[wei]++;
                        maxWei = max(maxWei, wei);
                    }
                }
                newCsz+=colorSz[col];
            }
        }
        //5. get the upper bound
        for (ui w = 0; w <=  maxWei; w++){
            if(kCur>=w*weiB[w]) {
                CMax+=weiB[w], kCur=kCur-w*weiB[w];
            }
            else {
                CMax+=kCur/w;
                break;
            }
        }
        UB=P_end+CMax;
    }
    
    
    if (true){
        if(S_UB< UB){
            printf("gap alg error!---------sort bound: %u, complex bound: %u\n", S_UB, P_end+CMax);
            printf("colNum=%u, Psz=%u, Csz=%u\n",maxCol+1, P_end, C_end-P_end);
            for (ui c = 0; c <= maxCol; c++){
                printf("color: %u\n", c);
                for (ui i = 0; i < colorSz[c]; i++){
                    ui u = colorVec[n*c+i];
                    printf("%u: w=%u, w_s=%u\n", colorVec[n*c+i], P_end-neiInP[u], weight[u]);
                }
            }
            for (ui w = 0; w <= K-MEInP; w++) printf("%u, ", weiB[w]);
            printf("\n");
            exit(0);
        }
        UB_gap  = S_UB - UB;
        gapUBSum += UB_gap;
    }
RETCOMPLEX: 
    if(!conflictPairs.empty()){
        for(auto p:conflictPairs){
            int u = p.first, v = p.second;
            MuEx[u][v]=false; MuEx[v][u]=false;
        }
        conflictPairs.clear();
    }
    nonNeiB.clear(); weiB.clear();weiPreSum.clear();
    return UB;
}
ui kDefectiveClique_BB::sortConfBoundL(){
    UB=0; ui maxWei=0, kCur = K-MEInP;
    ui maxColSz = 0;
    ui CMax = 0, S_UB = 0, UB_gap = 0;
    ui colNum=0; int maxCol=0;
    nonNeiB.resize(kCur+1); for(auto&b : nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    vector<ui> colvec0; colorVecL.push_back(colvec0);//push back the color 0 vector
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		if(ui(P_end-neiInP[u]) <= kCur)nonNeiB[ui(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= kCur; i++){
        for (auto u:nonNeiB[i]){
            int col=0;
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        ui lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            while(colvis[col]) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
                vector<ui> colvec; colorVecL.push_back(colvec);//push back the color 0 vector
			}
            colorLabel[u]=col; colorSz[col]++; colorVecL[col].push_back(u); maxColSz = max(maxColSz, colorSz[col]);
            ui wei=i+colorSz[col]-1; weiB[wei]++; maxWei=max(maxWei, wei);
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
    }
    //3. calculation of the prefix sum of the weight bucket 

    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
        // printf("w=%d, kCur=%d, weiB=%d\n", w, kCur, weiB[w]);
    }
    UB = P_end + CMax;
    S_UB = UB;
    if(S_UB <= LB || maxColSz <= 0.5*C_end ) goto RETSCBOUNDL;

    //4. calculate the sortconfbound
    weiB.assign(C_end, 0);
    confTreeCnt++;
    for (ui col = 0; col <= maxCol; col++){
        ui rem = colorSz[col];
        vector<ui> R;
        vector<vector<ui>> RVec;
        list<ui> colVecCpy(colorVecL[col].begin(), colorVecL[col].end());
        while(rem > 0){
            auto it = colVecCpy.begin();
            while (it != colVecCpy.end()){
                bool allconflict = true;  //to check if v in Pi_i is conflict to all other vertices in 
                ui v = *it;
                for(auto u: R){
                    if(is2hop(u, v) || true){
                        allconflict = false;
                        break;
                    }
                }
                if (allconflict){
                    R.push_back(*it);
                    it  = colVecCpy.erase(it);
                    rem--;
                }else ++it;
            }
            RVec.emplace_back(R);
            ui wei = P_end - neiInP[R[0]] + RVec.size() - 1;
            weiB[wei]++;
            R.clear();
        }
        RVec.clear(); colVecCpy.clear();
    }
    kCur = K - MEInP, CMax = 0;
    for (ui w = 0; w <=  maxWei; w++){
        if(kCur>=w*weiB[w]) {
            CMax+=weiB[w], kCur=kCur-w*weiB[w];
        }
        else {
            CMax+=kCur/w;
            break;
        }
    }
    UB = P_end + CMax;

    RETSCBOUNDL:
    UB_gap  = S_UB - UB; gapUBSum += UB_gap;
    nonNeiB.clear(); weiB.clear(); colorVecL.clear();
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    }
    return UB;
}
ui kDefectiveClique_BB::sortBoundL(){
    ui UB=0; ui maxWei=0;
    nonNeiB.resize(P_end+1); for(auto&b : nonNeiB) b.clear();
    weiB.resize(C_end); weiB.assign(C_end, 0);
    weiPreSum.resize(C_end-P_end); weiPreSum.assign(C_end- P_end, 0);
    ui colNum=0, maxCol=0; 
    colorSz[0]=0;
    //1. use bucket sorting to sort
    for (ui i = P_end; i < C_end; i++){
        ui u= PC[i];
		nonNeiB[ui(P_end-neiInP[u])].push_back(u);
    }
    //2. coloring the vertices in C
    for (ui i = 0; i <= P_end; i++){
        for (auto u:nonNeiB[i]){
            ui col=0;
            for (ui j = pstart[u]; j < pstart[u+1]; j++){
                ui nei=edges[j];
                if(inC(nei)){
                    if(colorLabel[nei]>=0){
                        ui lab=colorLabel[nei];
                        if(lab>=n){
                            printf("error label: %u\n",lab);
                        }
                        colvis[lab]=true;
                    }
                }
            }
            
            while(colvis[col]) col++;
			if(col>maxCol){
				maxCol=max(maxCol,col);
				colorSz[maxCol]=0;
			}
            colorLabel[u]=col;
            colorSz[col]++;
            ui wei=i+colorSz[col]-1;
            weiB[wei]++;
            maxWei=max(maxWei, wei);
            for (ui c = 0; c <= maxCol; c++) colvis[c]=false;
        }
        
    }
    //3. calculation of the prefix sum of the weight bucket 
    ui vIdx=0;
    for (ui w = 0; w <=  maxWei; w++){
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
    ui i = C_end;
    while (i>P_end && K < MEInP+weiPreSum[i-P_end-1]){
        i--;
    }
    UB=i;
    nonNeiB.clear(); weiB.clear(); weiPreSum.clear(); 
    for (ui i = P_end; i < C_end; i++){
        ui u=PC[i];
        colorLabel[u]=-1;
    }
    return UB;
}
void kDefectiveClique_BB::branch(ui level){
    assert(level<=n+1);
    ui pushCnt = (useRed==true? reduction(level):0) ;
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if(MEInP > K) goto REC;
    if (P_end > LB) store(P_end);
    if( MEInG <= K){
        //because LB>=1, we do not need to consider C_end<=1
        if(C_end<=1) printf("Error!\n"),exit(0);
        printf("G: %d, MEInG: %d, %.3f\n",C_end, MEInG, K);
        prune1++; store(C_end); goto REC;
    }
    if(sortBoundL()<=LB&&useUB){ub_prune++; goto REC;}
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return 
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                ept edge_sum=P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= P_end * (P_end + 1) - 2 * K){
                    must_include=true;
                    u=v; break;
                }
            }
        }
    }
    if(u==n) u=vtxSelect();

    CtoP(u,level);
    branch(level+1);// branch on adding the vertex u
    if(must_include){
        prune1++;
        PtoC(u,level);
        goto REC;
    }
    PtoX(u,level);
    branch(level+1);// branch on deleting the vertex u
    XtoC(u, level);
REC:
    recover(pushCnt, level);
    return;
}
void kDefectiveClique_BB::branch2(ui level){
    treeIdx++;
    // printf("enter the branch %ld\n", treeIdx);
    assert(level<=n+1);
    ui pushCnt = (useRed==true? reduction(level):0) ;
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if(MEInP > K) goto REC;
    if (P_end > LB && is2D(P_end)) store(P_end);
    if( MEInG <= K && is2D(C_end)){
        //because LB>=1, we do not need to consider C_end<=1
        if(C_end<=1) printf("Error!\n"),exit(0);
        prune1++; store(C_end); goto REC;
    }
    if(sortConfBoundL()<=LB&&useUB){ub_prune++; goto REC;}
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return 
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                ept edge_sum=P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= P_end * (P_end + 1) - 2 * K){
                    must_include=true;
                    u=v; break;
                }
            }
        }
    }
    if(u==n) u=vtxSelect();

    CtoP(u,level);
    branch2(level+1);// branch on adding the vertex u
    if(must_include){
        if(verify2hop(P_end)){
            prune1++;
            PtoC(u,level);
            goto REC;
        }
    }
    PtoX(u,level);
    branch2(level+1);// branch on deleting the vertex u
    XtoC(u, level);
REC:
    recover(pushCnt, level);
    return;
}
void kDefectiveClique_BB::branchSubG(ui level){
    treeIdx++;
    assert(level<=n+1);
    ui pushCnt = (useRed==true? reduction(level):0) ;
    ui u=n; bool must_include=false;
    if(C_end <= LB) goto REC;
    if(MEInP>K) goto REC;
    if (P_end > LB && ( P_end*(P_end-1)/2-MEInP ) >= P_end*(P_end-1)/2 - K ) {
         if(verify2hop(P_end)){
            store(P_end);
            assert(verifyKDC());
        }
    }
      //if G[P+C] is a k-defective clique
    if (( C_end * (C_end-1)/2 - MEInG ) >= ( C_end * (C_end - 1)/2 - K )) {
        printf("G: %d, MEInG: %d\n",C_end, MEInG);
        if(verify2hop(C_end)){
            prune1++; store(C_end); 
            assert(verifyKDC());
            goto REC;
        }
    }
    if(sortConfBound()<=LB&&useUB){ub_prune++;goto REC;}
    if(P_end>=C_end) goto REC; //if candidate set is empty, then return  
    treeCnt++;
    if (usePrune1){
        for (ui i = P_end; i < C_end; i++){
            ui v=PC[i];
            if(neiInG[v]>=C_end-2){
                double edge_sum=(double)P_end*(P_end-1)-2*MEInP+2*neiInP[v];
                if(edge_sum >= K * (double)P_end*(P_end+1)){
                    must_include=true;
                    u=v;
                    break;
                }
            }
        }
    }
    
    if(u==n) u=vtxSelect();

    CtoP(u,level);
    branchSubG(level+1);// branch on adding the vertex u
    if(must_include){
        if(verify2hop(P_end)){
            prune1++;
            PtoC(u,level);
            goto REC;
        }
    }
    PtoX(u,level);
    branchSubG(level+1);// branch on deleting the vertex u
    XtoC(u, level);

REC:
    recover(pushCnt, level);
    return;
}

void kDefectiveClique_BB::PtoC(ui u, ui level){
    assert(inP(u));
    swapID(PC_rid[u],--P_end);
    ui nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;//update the missing edges in P
    //update the neiInP of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]--;
    }
    
}

void kDefectiveClique_BB::CtoP(ui u, ui level){
    assert(inC(u));
    swapID(PC_rid[u], P_end++);
    ui nonNeiP=P_end-1-neiInP[u];
    MEInP+=nonNeiP;
    //update the neiInP of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]++;
    }
}

void kDefectiveClique_BB::PtoX(ui u, ui level){
    assert(inP(u));
    //move u from P to C
    swapID(PC_rid[u],--P_end);
    ui nonNeiP=P_end-neiInP[u];
    MEInP-=nonNeiP;
    //move u from C to X
    swapID(P_end, --C_end);
    ui nonNeiG=C_end-neiInG[u];
    m-=neiInG[u];
    MEInG-=nonNeiG;
    //update the neiInP and neiInG of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInP[v]--;
        neiInG[v]--;
    }
}
void kDefectiveClique_BB::CtoX(ui u, ui level){
    assert(inC(u));
    swapID(PC_rid[u],--C_end);
    ui nonNeiG = C_end - neiInG[u];
    m-=neiInG[u];
    MEInG-=nonNeiG;
    //update the neiInG of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v = edges[j];
        neiInG[v]--;
    }
    
}
void kDefectiveClique_BB::XtoC(ui u, ui level){
    assert(inX(u));
    //move u from X to C
    swapID(PC_rid[u],C_end++);
    ui nonNeiG=C_end-1-neiInG[u];
    m+=neiInG[u];
    MEInG+=nonNeiG;
     //update the neiInG of neighbors of u
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v=edges[j];
        neiInG[v]++;
    }
}

void kDefectiveClique_BB::store(ui newLB){
    KDC.resize(LB=newLB);
    for (ui i = 0; i < LB; i++) KDC[i]=PC[i];
}

bool kDefectiveClique_BB::is2D(ui newLB){
    if(newLB>=K+2) return true;
    bool is2hop=true;
    for (int i = 0; i < newLB; i++){
        int u = PC[i];
        unordered_set<int> neis;
        //get the neighbors of u in PC[0]~PC[newLB]
        for (int j = pstart[u]; j < pstart[u+1]; j++){
            int v = edges[j];
            if(PC_rid[v] >= 0 && PC_rid[v] < newLB) neis.insert(v);
        }
        for (int j = i+1; j < newLB; j++){
            int v = PC[j];
            is2hop = false;
            //check if v is in the neis of u
            if(neis.find(v)!=neis.end()){
                is2hop=true;
            }else{
                for (int k = pstart[v]; k < pstart[v+1]; k++){
                    int w = edges[k];
                    if(neis.find(w)!=neis.end()) {
                        is2hop = true;
                        break;
                    }
                }   
            }
            if(!is2hop) return false;
        }
        neis.clear();
    }
    return true;
}

bool kDefectiveClique_BB::is2hop(ui u, ui v){
    // if (neiInG[u]+neiInG[v]>C_end-2) return true;
    unordered_set<ui> u_neis;
    for (ept i = pstart[u]; i < pstart[u+1]; i++){//find the neighbors of the vertex u
        ui u_nei = edges[i];
        if(v==u_nei) return true;
        if(PC_rid[u_nei] >= 0 && PC_rid[u_nei] < C_end) u_neis.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pstart[v+1]; i++){
        ui v_nei = edges[i];
        if(u_neis.find(v_nei) != u_neis.end()) return true;
    }
    return false;
}
bool kDefectiveClique_BB::is2hop2(ui u, ui v){
    if(isAdj(u, v)) return true;
    for(int i = 0; i < C_end; i++){
        ui w = PC[i];
        if(isAdj(u, w) && isAdj(v, w)) return true;
    }
    return false;
}
void kDefectiveClique_BB::swapID(ui i,ui j){
    swap(PC[i],PC[j]);
    PC_rid[PC[i]] = i;
	PC_rid[PC[j]] = j;
}

bool kDefectiveClique_BB::inP(ui u){
    return (PC_rid[u]>=0 && PC_rid[u]<P_end);
}

bool kDefectiveClique_BB::inC(ui u){
    return (PC_rid[u]>=P_end && PC_rid[u]<C_end);
}

bool kDefectiveClique_BB::inX(ui u){
    return (PC_rid[u]>=C_end && PC_rid[u]<n);
}

bool& kDefectiveClique_BB::isAdj(ui u,ui v){
    return matrix[u*n+v];
}





#endif