#include "Graph.h"
#include "Heuristc.h"
using namespace std;

Graph::Graph(const char *_dir, const double _GAMMA){
    dir=string(_dir);
    gamma=_GAMMA;
    n=m=0;
    maxDeg=maxCore=max2HpDeg=max2HpCore=0;
    maxSub = 0;
    pstart=nullptr;
    edges=nullptr;
    degree=nullptr;
    pstart2Hp=nullptr;
    neis2Hp=nullptr;
    deg2Hp=nullptr;
    MDS.clear();

}
void Graph::read(){
    FILE *f=fopen(dir.c_str(), "rb");
    fread(&n, sizeof(int), 1, f); fread(&n, sizeof(int), 1, f); fread(&m, sizeof(int), 1,f);
    printf("\t(orgin): n = %u; m = %llu (undirected)\n",n , m/2);
    degree=new ui[n];
    // printf("read the degree begin\n");
    fread(degree, sizeof(int), n,f);
    // printf("read the degree end\n");
    if(pstart == nullptr) pstart = new ept[n+1];
    if(edges == nullptr) edges= new ui[m];

    pstart[0]=0;
    for(ui u=0; u<n; u++){
        if(degree[u]>0){
            fread(edges+pstart[u], sizeof(int), degree[u],f);
            ui* edge_u=edges+pstart[u];
            sort(edge_u, edge_u+degree[u]);
            ui idx=0;
            for(ui vLoc=0; vLoc<degree[u];vLoc++){
                if(edge_u[vLoc]>= n) printf("vertex id %u out of range\n", edge_u[vLoc]);
                if(edge_u[vLoc] == u || (vLoc>0 && edge_u[vLoc]==edge_u[vLoc-1])) continue;//delete the self-loops and the loops
                edge_u[idx++]=edge_u[vLoc];
            }
            degree[u]=idx;
        }
        pstart[u+1]=pstart[u]+degree[u];
    }
    m=pstart[n];
    printf("\t(processed): n = %u; m = %llu (undirected)\n",n , m/2);
    fclose(f);

}
void Graph::search(){
    if(n<=0){
        printf("Max Dense Subgraph is %u\n",0);
    }
    bool weakDegen=true;
    bool turingRed=true;
    ui *core=new ui[n];
    ui *seq = new ui[n];
    char* deleted = new char[n]; memset(deleted, 0, n*sizeof(char));
    for (ui u = 0; u < n; u++) seq[u]=u; 
    degen(n, seq, core, pstart, edges, degree, true);
    HeuriSearcher *heuri_solver=new HeuriSearcher(n,m, pstart,edges,gamma, maxDeg, MDS);
    heuri_solver->degenSearch(MDS, seq);
    if(weakDegen){
        deg2Hp = new ui[n];
        ui* core2Hp = new ui[n];
        
        build2Hpdeg(pstart, edges, deg2Hp, deleted);
        weakdegen(n, seq, core2Hp, pstart, edges, deg2Hp, true);
        // reorder the seq with weak degeneracy
        delete[] deg2Hp;
        delete[] core2Hp;
    }
    delete[] core;
    if(true){
        //enter the search
        
        if(turingRed){
            ui* ids = new ui[n]; ui sub_n=0;
            ui* rid = new ui[n];
            vector<pair<ui, ui>> vp;
            
            char* exists = new char[n]; memset(exists, 0, n*sizeof(char));
            for (ui i = 0; i < n; i++){
                ui u = seq[i];
                ui pre_size = ui(MDS.size());
                induceSubgraph(u, seq, pstart, edges, deleted, exists, ids, rid, sub_n, vp);
                maxSub = max(maxSub, sub_n);
                // printf("vtx: %u, edg: %llu\n", sub_n, ept(vp.size()));
                // exit(0);
                // if(i == 5) exit(0);

                deleted[u] = 1;
            }
            printf("max size of subgraph: %u\n", maxSub);

            delete[] ids;
            delete[] rid;
            
            delete[] exists;
        }
        delete[] deleted;

        
    }
}
void Graph::induceSubgraph(ui u, ui* seq, ept* pstart, ui* edges, char* deleted, char* exists,ui* ids, ui* rid, ui& sub_n, std::vector<std::pair<ui,ui>>& vp){
    vp.clear();

    //find the vertices of the subgraph
    sub_n=0; ids[sub_n++]=u; exists[u]=1;rid[u]=sub_n-1;
    for (ept i = pstart[u]; i < pstart[u+1]; i++){
        ui v = edges[i];
        if(!deleted[v]){
            ids[sub_n++]=v; exists[v]=2; rid[v]=sub_n-1;
        }
    }
    ui nb_size=sub_n-1;
    for (ui i = 1; i < 1+nb_size; i++){
        ui v = ids[i];
        for (ept j = pstart[v]; j < pstart[v+1]; j++){
            ui w = edges[j];
            if (!exists[w] && !deleted[w]){
                ids[sub_n++]=w;
                exists[w]= 3;
                rid[w]=sub_n-1;
            }
        }
    }
    
    //find the edges of the subgraph
    for (ui i = 0; i < sub_n; i++){
        ui v = ids[i];
        for (ept j = pstart[v]; j < pstart[v+1]; j++){
            ui w = edges[j];
            if(exists[w] && w>v) vp.pb(make_pair(rid[v], rid[w]));
        }
    }
    for (int i = 0; i < sub_n; i++) exists[ids[i]] = 0;
}
void Graph::degen(ui n, ui *seq, ui *core, ept *pstart, ui *edges, ui* degree, bool output){
    Timer t;
    if(n > 0){
        char* vis=new char[n];
        memset(vis, 0, n*sizeof(char));
        ListLinearHeap *heap=new ListLinearHeap(n, n-1);
        heap->init(n, n-1, seq, degree);
        for (ui i = 0; i < n; i++){
            ui u, key; heap->pop_min(u, key);
            if(key> maxCore) maxCore = key;
            core[u]=maxCore;
            seq[i]=u, vis[u]=1;
            for (ept j = pstart[u]; j < pstart[u+1]; j++){
                ui v=edges[j];
                if(vis[v]) continue;
                heap->decrement(v, 1);
            }
        }
    }
    if(output) printf("*** MaxCore: %d, degen Time: %.2f ***\n",  maxCore, double(t.elapsed())/1000000);
    return;
}
void Graph::get2HpNei(ept* pstart, ui* edges, ui u, std::vector<int>& u_2HpNei, char* vis, char* deleted){
    vis[u]=1;
    for (ept j = pstart[u]; j < pstart[u+1]; j++){
        ui v = edges[j];
        if(!deleted[v]) u_2HpNei.push_back(v), vis[v]=1;
    }
    ui deg_u=u_2HpNei.size();
    for (int i = 0; i < deg_u; i++){
        ui v = u_2HpNei[i];
        // find the neighbors of v
        for (ept k = pstart[v]; k < pstart[v+1]; k++){
            ui w = edges[k];
            if(vis[w] || deleted[w]) continue;
            u_2HpNei.push_back(w), vis[w]=1;
        }
    }
}
void Graph::build2Hpdeg(ept* pstart, ui* edges, ui* deg2Hp, char* deleted){
    Timer t;
    char* vis=new char[n];
    memset(vis, 0, n*sizeof(char));
    std::vector<int> nei2Hp;
    for (ui u = 0; u < n; u++){
        get2HpNei(pstart, edges, u, nei2Hp, vis, deleted);
        deg2Hp[u] = ui(nei2Hp.size());
        max2HpDeg=max(max2HpDeg, deg2Hp[u]);
        //recover the vis array
        for(auto v: nei2Hp) vis[v]=0;
        vis[u]=0;
        nei2Hp.clear();
    }
    printf("The max 2hop is: %d, the time of 2hop construction is: %.2f\n", max2HpDeg,  double(t.elapsed())/1000000);
}
void Graph:: weakdegen(ui n, ui* seq, ui* core2Hp, ept* pstart, ui* edges, ui *deg2Hp, bool output){
    Timer t;
    char* deleted=new char[n]; memset(deleted, 0, n*sizeof(char));
    char* vis = new char[n]; memset(vis, 0, n*sizeof(char));
    std::vector<int> nei2Hp;
    ListLinearHeap *heap=new ListLinearHeap(n, n-1);
    for (ui i = 0; i < n; i++) seq[i]=i;
    heap->init(n, n-1, seq, deg2Hp);
    for (ui i = 0; i < n; i++){
        ui u, key; heap->pop_min(u, key);
        if(key > max2HpCore) max2HpCore = key;
        core2Hp[u]= max2HpCore;
        seq[i] = u; deleted[u]=1;
        //renew the deg2Hp of the each 2hop neighbors of u
        get2HpNei(pstart, edges, u, nei2Hp, vis, deleted);
        for(auto v: nei2Hp){
           if(!deleted[v]) heap->decrement(v, 1);
           vis[v]=0;
        }
        vis[u]=0;
        nei2Hp.clear();
        
    }
    ui UB=max2HpCore+1;
    if(output) printf("MaxCore2hop: %d, UB: %d, Order Time: %.2f\n", max2HpCore, UB, double(t.elapsed())/1000000);
}
Graph::~Graph()
{
    if(pstart!=nullptr){
        delete[] pstart;
        pstart=nullptr; 
    }
    if(edges!= nullptr){
        delete[] edges;
        edges=nullptr;
    }
    if(degree != nullptr){
        delete[] degree;
        degree = nullptr;
    }
    if(pstart2Hp != nullptr){
        delete[] pstart2Hp;
        pstart2Hp = nullptr;
    }
    if(neis2Hp != nullptr){
        delete[] neis2Hp;
        neis2Hp = nullptr;
    }
    if(deg2Hp != nullptr){
        delete[] deg2Hp;
        deg2Hp = nullptr;
    }
    if(!MDS.empty()) MDS.clear();
}