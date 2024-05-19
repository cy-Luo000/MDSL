#include "Graph.h"
#include "Heuristc.h"
#include "QuasiClique_BB.h"
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
Graph::Graph(const char *_dir, const int _K){
    dir=string(_dir);
    K=_K;
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
    // printf("\t(orgin): n = %u; m = %llu (undirected)\n",n , m/2);
#ifdef _TEST_
	printf("#n=%u\n#m=%llu\n", n, m/2);
#else 
	printf("\t(orgin): n = %u; m = %llu (undirected)\n",n , m/2);
#endif
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
            degree[u]=idx; maxDeg=max(maxDeg, degree[u]);

        }
        pstart[u+1]=pstart[u]+degree[u]; 
    }
    m=pstart[n];
    printf("\t(processed): n = %u; m = %llu (undirected)\n",n , m/2);
    fclose(f);

}
void Graph::search(ui _DSkind){
    Timer t;
    if(n<=0) printf("Max Dense Subgraph is %u\n",0);

    bool weakDegen=true;
    bool turingRed=true;
    bool useHeu=true;
    ui UB = n;
    ui *core=new ui[n];
    ui *seq = new ui[n];
    char* deleted = new char[n]; memset(deleted, 0, n*sizeof(char));// the array to record the deleted vertices
    char* vis = new char[n]; memset(vis, 0, n*sizeof(char));// the array to record the visited vertices
    for (ui u = 0; u < n; u++) seq[u]=u; 
    degen(n, seq, core, pstart, edges, degree, vis,true);
    if(useHeu){
        if(_DSkind==1){
            HeuriSearcher *heuri_solver=new HeuriSearcher(n,m, pstart,edges,gamma, maxDeg, MDS);
            heuri_solver->degenSearch(MDS, seq, _DSkind);
            delete heuri_solver;
        }else{
            HeuriSearcher *heuri_solver=new HeuriSearcher(n,m, pstart,edges, K, maxDeg, MDS);
            heuri_solver->degenSearch(MDS, seq, _DSkind);
            delete heuri_solver;
        }
    }
    // if(_DSkind==2) exit(0);
    if(weakDegen){
        deg2Hp = new ui[n];
        ui* core2Hp = new ui[n];
        
        build2Hpdeg(pstart, edges, deg2Hp, vis);
        UB = weakdegen(n, seq, core2Hp, pstart, edges, deg2Hp, deleted, vis, true);
        // reorder the seq with weak degeneracy
        delete[] deg2Hp;
        delete[] core2Hp;
    }
    delete[] core;
    Timer tt;
    if(UB>(ui)MDS.size()){
        //enter the search
        if (_DSkind==1){
            if(turingRed){
                ui* ids = new ui[n]; ui sub_n=0;
                ui* rid = new ui[n];
                vector<pair<ui, ui>> vp;
                // char* exists = new char[n]; memset(exists, 0, n*sizeof(char));
                for (ui i = 0; i < n; i++){
                    ui u = seq[i];
                    ui pre_size = ui(MDS.size());
                    induceSubgraph(u, seq, pstart, edges, deleted, vis, ids, rid, sub_n, vp);
                    maxSub = max(maxSub, sub_n);
                    //begin the search of the subgraphs
                    if(sub_n>pre_size){
                        QuasiClique_BB *MQCSolver=new QuasiClique_BB();
                        MQCSolver->load_subgraph(gamma, sub_n, vp, MDS,UB);
                        MQCSolver->MQCSearch2hop(MDS);
                        //update the best solution
                        if(MDS.size() > pre_size) for(ui j = 0; j < (ui)MDS.size(); j++) MDS[j] = ids[MDS[j]];
                        delete MQCSolver;
                        
                        
                    }
                    deleted[u] = 1;
                }
                // printf("max size of subgraph: %u\n", maxSub);
                delete[] ids;
                delete[] rid;
                
                // delete[] exists;
            }else{
                QuasiClique_BB *MQCSolver=new QuasiClique_BB();
                MQCSolver->load_graph(n,pstart,pstart+1,edges);
    #ifndef _TEST_
                printf("graph load success!\n");
    #endif
                MQCSolver->MQCSearch(gamma, UB, MDS);
    #ifndef _TEST_
                printf("BBSearch complete, tree count: %lld\n", treeCnt);
    #endif
                    // exit(0);
                delete MQCSolver;
            }


#ifdef _TEST_
		// printf("#MaxQCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", MDS.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("#maxP=%d\n#minPUB=%d\n#maxME=%d\n", max_P_end, P_UBMin,maxME);
		// printf("#NodeCount=%lld\n",treeCnt);
		printf("#MaxSG=%d\n",maxSub);
		// printf("#FeasibleSubgraph=%d\n",feasible);
		// printf("#prune1=%lld\n#ubprune=%lld\n", prune1,ub_prune);
		// printSubInfo();
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("Max Mutual Exclusive: %d, Mutual exclusive sum: %lld, Mutual exclusive avg: %lf, Mutual exclusive dense avg: %lf\n",MaxMuExNum, MuExSum, (double)MuExSum/(tree_cnt+0.1), denSum/(denNum+0.001));
		// printf("Sum of Mutual Exclusive in color set: %lld, max mutual exclusive in color set: %d\n", colMuExSum, maxColMuNum);
		printf("Search Tree Size: %lld\n",treeCnt);
		printf("Feasible Subgraph: %d\n",feasible);
#endif
        }
        else{
            if(turingRed){
                ui* ids = new ui[n]; ui sub_n=0;
                ui* rid = new ui[n];
                vector<pair<ui, ui>> vp;
                // char* exists = new char[n]; memset(exists, 0, n*sizeof(char));
                for (ui i = 0; i < n; i++){
                    ui u = seq[i];
                    ui pre_size = ui(MDS.size());
                    induceSubgraph(u, seq, pstart, edges, deleted, vis, ids, rid, sub_n, vp);
                    //begin the search of the subgraphs
                    if(sub_n>pre_size){
                        kDefectiveClique_BB *MKDCSolver=new kDefectiveClique_BB();
                        MKDCSolver->load_subgraph(K, sub_n, vp, MDS,UB);
                        MKDCSolver->MKDCSearch2hop(MDS);
                        //update the best solution
                        if(MDS.size() > pre_size) for(ui j = 0; j < (ui)MDS.size(); j++) MDS[j] = ids[MDS[j]];
                        delete MKDCSolver;
                        maxSub = max(maxSub, sub_n);
                        
                    }
                    deleted[u] = 1;
                }
                // printf("max size of subgraph: %u\n", maxSub);
                delete[] ids;
                delete[] rid;
                
                // delete[] exists;
            }else{
                kDefectiveClique_BB *MKDCSolver=new kDefectiveClique_BB();
                MKDCSolver->load_graph(n,pstart,pstart+1,edges);
    #ifndef _TEST_
                printf("graph load success!\n");
    #endif
                MKDCSolver->MKDCSearch(K, UB, MDS);
    #ifndef _TEST_
                printf("BBSearch complete, tree count: %lld\n", treeCnt);
    #endif
                    // exit(0);
                delete MKDCSolver;
            }


#ifdef _TEST_
		// printf("#MaxQCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", MDS.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("#maxP=%d\n#minPUB=%d\n#maxME=%d\n", max_P_end, P_UBMin,maxME);
		// printf("#NodeCount=%lld\n",treeCnt);
		printf("#MaxSG=%d\n",maxSub);
		// printf("#FeasibleSubgraph=%d\n",feasible);
		// printf("#prune1=%lld\n#ubprune=%lld\n", prune1,ub_prune);
		// printSubInfo();
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("Max Mutual Exclusive: %d, Mutual exclusive sum: %lld, Mutual exclusive avg: %lf, Mutual exclusive dense avg: %lf\n",MaxMuExNum, MuExSum, (double)MuExSum/(tree_cnt+0.1), denSum/(denNum+0.001));
		// printf("Sum of Mutual Exclusive in color set: %lld, max mutual exclusive in color set: %d\n", colMuExSum, maxColMuNum);
		printf("Search Tree Size: %lld\n",treeCnt);
		printf("Feasible Subgraph: %d\n",feasible);
#endif
        }
        delete[] deleted;
        delete[] vis;
        
    }
    #ifdef _TEST_
        if(_DSkind==1) printf("#MaxQCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", MDS.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
        else printf("#MaxKDCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", MDS.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("#maxP=%d\n#minPUB=%d\n#maxME=%d\n", max_P_end, P_UBMin,maxME);
		printf("#NodeCount=%lld\n",treeCnt);
		// printf("#MaxSG=%d\n",maxSubSz);
		// printf("#FeasibleSubgraph=%d\n",feasible);
		// printf("#prune1=%lld\n#ubprune=%lld\n", prune1,ub_prune);
		// printSubInfo();
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("Max Mutual Exclusive: %d, Mutual exclusive sum: %lld, Mutual exclusive avg: %lf, Mutual exclusive dense avg: %lf\n",MaxMuExNum, MuExSum, (double)MuExSum/(tree_cnt+0.1), denSum/(denNum+0.001));
		// printf("Sum of Mutual Exclusive in color set: %lld, max mutual exclusive in color set: %d\n", colMuExSum, maxColMuNum);
		printf("Search Tree Size: %lld\n",treeCnt);
		printf("Feasible Subgraph: %d\n",feasible);
#endif
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
void Graph::degen(ui n, ui *seq, ui *core, ept *pstart, ui *edges, ui* degree, char* vis,bool output){
    Timer t;
    if(n > 0){
        // char* vis=new char[n];
        // memset(vis, 0, n*sizeof(char));
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
    // if(output) printf("*** MaxCore: %d, degen Time: %.2f ***\n",  maxCore, double(t.elapsed())/1000000);
    memset(vis, 0, n*sizeof(char));// recover the visited array
    #ifdef _TEST_
		if(output) printf("#MaxCore=%d\n#DegenTime=%.2f\n", maxCore, double(t.elapsed())/1000000);
		printf("#MaxDeg=%d\n",maxDeg);
    #else
        if(output) printf("*** MaxCore: %d, degen Time: %.2f ***\n",  maxCore, double(t.elapsed())/1000000);
    #endif
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
void Graph::build2Hpdeg(ept* pstart, ui* edges, ui* deg2Hp, char* vis){
    Timer t;
    std::vector<int> nei2Hp;
    ui min2HpDeg=n+1;
    for (ui u = 0; u < n; u++){
        vis[u]=1;
        for (ept j = pstart[u]; j < pstart[u+1]; j++){
            ui v = edges[j];
            // if(vis[v]) continue;
            nei2Hp.push_back(v), vis[v]=1;
        }
        ui deg_u=nei2Hp.size();
        for (int i = 0; i < deg_u; i++){
            ui v = nei2Hp[i];
            // find the neighbors of v
            for (ept k = pstart[v]; k < pstart[v+1]; k++){
                ui w = edges[k];
                if(vis[w]) continue;
                nei2Hp.push_back(w), vis[w]=1;
            }
        }

        deg2Hp[u] = ui(nei2Hp.size());
        max2HpDeg=max(max2HpDeg, deg2Hp[u]);
        min2HpDeg=min(min2HpDeg, deg2Hp[u]);
        //recover the vis array
        for(auto v: nei2Hp) vis[v]=0;
        vis[u]=0;
        nei2Hp.clear();
    }
    // printf("The max 2hop is: %d, the time of 2hop construction is: %.2f\n", max2HpDeg,  double(t.elapsed())/1000000);
#ifdef _TEST_
	printf("#max2hop=%d\n#min2hop=%d\n#2hopbuildtime=%.2f\n", max2HpDeg, min2HpDeg,double(t.elapsed())/1000000);
#elif
	printf("The max 2hop is: %d, the min 2hop is: %d, the time of 2hop construction is: %.2f\n", max2hop, min2hop, double(t.elapsed())/1000000);
#endif
}
ui Graph:: weakdegen(ui n, ui* seq, ui* core2Hp, ept* pstart, ui* edges, ui *deg2Hp, char* deleted, char* vis,bool output){
    Timer t;
    // char* deleted=new char[n]; memset(deleted, 0, n*sizeof(char));
    // char* vis = new char[n]; memset(vis, 0, n*sizeof(char));
    // char* 
    std::vector<ui> nei2Hp;//2 hop neighbors
    std::vector<ui> neis;//neighbors(1 hop neighbors)
    std::vector<ui> nei2Hp_v;// the 2hop neighbors of v(v is the neighbor of u)
    ListLinearHeap *heap=new ListLinearHeap(n, n-1);
    for (ui i = 0; i < n; i++) seq[i]=i;
    heap->init(n, n-1, seq, deg2Hp);
    for (ui i = 0; i < n; i++){
        ui u, key; heap->pop_min(u, key);
        if(key > max2HpCore) max2HpCore = key;
        core2Hp[u]= max2HpCore;
        seq[i] = u; deleted[u]=1;
        //renew the deg2Hp of the each 2hop neighbors of u
        // get2HpNei(pstart, edges, u, nei2Hp, vis, deleted);
        assert(neis.empty());
        assert(nei2Hp.empty());
        assert(nei2Hp_v.empty());
        vis[u]=1;
        for (ept j = pstart[u]; j < pstart[u+1]; j++){
            ui v = edges[j];
            if(!deleted[v]) neis.push_back(v), vis[v]=1;
        }
        // ui deg_u=(ui)nei2Hp.size();
        for (ui j = 0; j < neis.size(); j++){
            ui v = neis[j];
            // find the neighbors of v
            for (ept k = pstart[v]; k < pstart[v+1]; k++){
                ui w = edges[k];
                if(vis[w] || deleted[w]) continue;
                nei2Hp.push_back(w), vis[w]=1;
            }
        }
        vis[u]=0;
        for(auto v:nei2Hp) vis[v]=0;
        for(auto v:neis) vis[v]=0;
        for(auto v: nei2Hp){
           if(!deleted[v]) heap->decrement(v, 1);
        //    vis[v]=0;
        }
        // calculate the 2hop neighbors of the neighbor of u from the beginning
        for(auto v:neis){
            vis[v] =1;
            for (ept j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(!deleted[w]){
                    nei2Hp_v.push_back(w), vis[w]=1;
                }
            }
            ui deg_v = (ui)nei2Hp_v.size();
            for (ui j = 0; j < deg_v; j++){
                ui w = nei2Hp_v[j];
                for (ept k = pstart[w]; k < pstart[w+1]; k++){
                    ui x = edges[k];
                    if(!deleted[x] && !vis[x]){
                        nei2Hp_v.push_back(x), vis[x]=1;
                    }
                }
            }
            ui curDeg2Hp=heap->get_key(v);
            ui dec = curDeg2Hp - (ui)nei2Hp_v.size();
            heap->decrement(v, dec);
            for(auto w: nei2Hp_v) vis[w]=0;
            vis[v]=0;
            nei2Hp_v.clear();
        }
        
        nei2Hp.clear();
        neis.clear();
        
    }
    ui UB=max2HpCore+1;
    memset(deleted, 0, n*sizeof(char));//recover the deleted array
    // if(output) printf("MaxCore2hop: %d, UB: %d, Order Time: %.2f\n", max2HpCore, UB, double(t.elapsed())/1000000);
#ifdef _TEST_
	if(output) printf("#MaxCore2hop=%d\n#UB=%d\n#Order2hopTime=%.2f\n", max2HpCore, UB, double(t.elapsed())/1000000);
#elif
	if(output) printf("MaxCore2hop: %d, UB: %d, Order Time: %.2f\n", max2HpCore, UB, double(t.elapsed())/1000000);
#endif
    return UB;
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