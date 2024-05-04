#include "Graph.h"
#include "QuasiClique_BB.h"
#include "MIP.h"

#define _TEST_
using namespace std;
Graph::Graph(const char *_dir, const long double _GAMMA) {
	dir = string(_dir);
	K = INT32_MAX;
	gamma=_GAMMA;
	n = m = 0;
	max2hopCore=0;
	maxCore=0;
	maxDeg=0;
	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	pstart2hop=nullptr; adj2hop.clear();
	KDC.clear();
	vector<int> subSz;
	s_degree = s_edges = NULL;
	s_pstart = s_pend = NULL;
	s_peel_sequence = s_core = NULL;
	s_vis = NULL;
	s_heap = NULL;

	s_edgelist_pointer = NULL;
	s_tri_cnt = s_edge_list = NULL;
	s_active_edgelist = NULL;
	s_deleted = NULL;
	
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pend_buf != NULL) {
		delete[] pend_buf;
		pend_buf = NULL;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(pstart2hop!= nullptr){
		delete[] pstart2hop;
		pstart2hop=nullptr;
	}
	if(!adj2hop.empty()){
		adj2hop.clear();
	}
	if(!subSz.empty()){
		subSz.clear();
	}
	if(s_degree != NULL) {
		delete[] s_degree;
		s_degree = NULL;
	}
	if(s_pstart != NULL) {
		delete[] s_pstart;
		s_pstart = NULL;
	}
	if(s_pend != NULL) {
		delete[] s_pend;
		s_pend = NULL;
	}
	if(s_edges != NULL) {
		delete[] s_edges;
		s_edges = NULL;
	}
	if(s_peel_sequence != NULL) {
		delete[] s_peel_sequence;
		s_peel_sequence = NULL;
	}
	if(s_core != NULL) {
		delete[] s_core;
		s_core = NULL;
	}
	if(s_vis != NULL) {
		delete[] s_vis;
		s_vis = NULL;
	}
	if(s_heap != NULL) {
		delete s_heap;
		s_heap = NULL;
	}
	if(s_edgelist_pointer != NULL) {
		delete[] s_edgelist_pointer;
		s_edgelist_pointer = NULL;
	}
	if(s_active_edgelist != NULL) {
		delete[] s_active_edgelist;
		s_active_edgelist = NULL;
	}
	if(s_deleted != NULL) {
		delete[] s_deleted;
		s_deleted = NULL;
	}
}

void Graph::read() {
	FILE *f = fopen(dir.c_str(), "rb");
	fread(&n, sizeof(int), 1, f); fread(&n, sizeof(int), 1, f); fread(&m, sizeof(int), 1, f);
#ifdef _TEST_
	printf("#n=%d\n#m=%d\n", n, m/2);
#else 
	printf("\tn = %d; m = %d (undirected)\n", n, m/2);
#endif
	

	int *degree = new int[n];
	fread(degree, sizeof(int), n, f);
	if(pstart == nullptr) pstart = new int[n+1];
	if(edges == nullptr) edges = new int[m];

	pstart[0] = 0;
	for(int i = 0;i < n;i ++) {
		if(degree[i] > 0) {
			fread(edges+pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			int *buff = edges+pstart[i];
			sort(buff, buff+degree[i]);
			int idx = 0;
			for(int j = 0;j < degree[i];j ++) {
				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;//buff[j]==i-->exist selfloop, buff[j]==buff[j-1]-->exist parallel edges
				buff[idx ++] = buff[j];
			}
			degree[i] = idx; maxDeg=max(maxDeg, degree[i]);
		}

		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);
	delete[] degree;
}
int Graph::twoHopDegConstruct(int *pstart, int *edges, int *deg2hop){
	Timer t;

	// std::vector<int> degvec;

	pstart2hop=new int[n+1];
	bool *vis=new bool[n];
	memset(vis, false, n*sizeof(bool));
	
	int max2hop=0, min2hop=n;
	vector<int> nei2hop;
	int idx=0;
	for (int u = 0; u < n; u++){
		pstart2hop[u]=idx;
		//move vertex i's neighbors into nei2hop
		vis[u]=true;
		for (int j = pstart[u]; j < pstart[u+1]; j++){
			int v=edges[j];
			if(vis[v]) continue;
			// assert(vis[v]==false);
			nei2hop.push_back(v);
			vis[v]=true;
		}
		//move vertex i's 2hop neighbors into nei2hop
		int curDeg=pstart[u+1]-pstart[u];
		for (int j = 0; j < curDeg; j++){
			int v=nei2hop[j];
			//find the neighbors of u
			for (int k = pstart[v]; k < pstart[v+1]; k++){
				int w=edges[k];
				if(vis[w]) continue;
				nei2hop.push_back(w);
				vis[w]=true;
			}
		}
		for (int i = 0; i < nei2hop.size(); i++){
			int v=nei2hop[i];
			adj2hop.push_back(v);
			idx++;
		}
		deg2hop[u]=idx-pstart2hop[u];
		max2hop=max(max2hop, deg2hop[u]), min2hop=min(min2hop, deg2hop[u]);
		// printf("%d: %d\n", u, deg2hop[u]);
		for(auto v:nei2hop) vis[v]=false;
		vis[u]=false;
		nei2hop.clear();
	}
	pstart2hop[n]=idx;
	// for (int i = 0; i < n; i++) printf("deg2hop: %d\n",deg2hop[i]);
	assert(idx==adj2hop.size());
#ifdef _TEST_
	printf("#max2hop=%d\n#min2hop=%d\nadj2hop=%ld\n#2hopbuildtime=%.2f\n", max2hop, min2hop, adj2hop.size(),double(t.elapsed())/1000000);
#elif
	printf("The max 2hop is: %d, the min 2hop is: %d, the size of adj2hop is: %ld, the time of 2hop construction is: %.2f\n", max2hop, min2hop, adj2hop.size(),double(t.elapsed())/1000000);
#endif
	return 0;
}
void Graph::write() {
	FILE *fout = fopen("KDC.txt", "w");
	fprintf(fout, "%d\n", int(KDC.size()));
	sort(KDC.begin(), KDC.end());
	for(int i = 0;i < KDC.size();i ++) fprintf(fout, "%d ", KDC[i]);
	fclose(fout);
}

void Graph::search() {
	// printf("enter search\n");

	Timer t;

	KDC.resize(0); //screen out trivial cases
	int *seq = new int[n];
	int *core = new int[n];
	int *deg = new int[n];
	int *deg2hop=new int[n];
	// int *deg2hop2=new int[n];
	int *core2hop=new int[n];
	// int *core2hop2=new int[n];
	char *vis = new char[n];

	ListLinearHeap *heap = new ListLinearHeap(n, n-1);
	ListLinearHeap *heap2 = new ListLinearHeap(n, n-1);
	// int UB = degen(n, seq, core, pstart, edges, deg, vis, heap, true);
	int UB=n; bool use2hop=true; bool useDegen=false;
	bool useMIP=false;
#ifdef _TEST_
	if (useMIP)
	{
		if(use2hop){
			printf("#Solver=GF3-2\n");
		}else{
			printf("#Solver=GF3\n");
		}
	}else{
		if (use2hop)
		{
			if(!useHeu){
				printf("#Solver=GA2NH\n");
			}else{
				if(useDegen){
					printf("#Solver=GA-2D\n");
				}
				else{
					if(useUB && !usePrune1){
						printf("#Solver=GA-2\n");
					}else if (!useUB && usePrune1)
					{
						printf("#Solver=GA2NBD\n");
					}else if(useUB && usePrune1){
						printf("#Solver=GA2+\n");
					}
				}
			}
		}else{
			printf("#Solver=GA\n");
		}
	}
	
#endif
	if(use2hop) {
		// printf("2 hop degen order\n");
		if (useDegen)
		{
			// UB=degen(n, seq, core, pstart, edges, deg, vis, heap, true);
			degen1hop(n, seq, core, pstart, edges, deg, vis, heap2, true);
			if(use2hop && useHeu){
				HeuriSearcher* heuri_solver=new HeuriSearcher(n,m,pstart, edges, gamma, maxDeg, KDC);
				heuri_solver->degenSearch(KDC); // Your code Here
			}
		}else{
			if(n<=5000){
				// printf("enter 2hop construct\n");
				degen1hop(n, seq, core, pstart, edges, deg, vis, heap2, true);
				if(use2hop && useHeu){
					HeuriSearcher* heuri_solver=new HeuriSearcher(n,m,pstart, edges, gamma, maxDeg, KDC);
					heuri_solver->degenSearch(KDC); // Your code Here
				}
				twoHopDegConstruct(pstart, edges, deg2hop);
				// printf("enter 2hop construct\n");
				// exit(0);
				UB=degenBy2hopDeg(n, seq, core2hop, pstart2hop, adj2hop.data(), deg2hop,vis, heap, true);//return ordered seq
			}
			else{
				// twoHopDegConstruct(pstart, edges, deg2hop2);
				// degenBy2hopDeg(n, seq, core2hop2, pstart2hop, adj2hop.data(), deg2hop2,vis, heap2, true);
				degen1hop(n, seq, core, pstart, edges, deg, vis, heap2, true);
				if(use2hop && useHeu){
					HeuriSearcher* heuri_solver=new HeuriSearcher(n,m,pstart, edges, gamma, maxDeg, KDC);
					heuri_solver->degenSearch(KDC); // Your code Here
				}
				UB=degen2hop4Large(n,seq, core2hop,pstart, edges, deg2hop, vis, heap, true);
			}
			// exit(0);
		}
		
	}else{
		// printf("degen order\n");
		UB=degen(n, seq, core, pstart, edges, deg, vis, heap, true);
		if(UB<0){
					printf("degen error\n");
		}
	}
	// exit(0);
	delete heap;
	delete[] vis;
	delete[] deg;
	adj2hop.clear();
	// delete[] deg2hop;
	// printf("enter heuristc search\n");
	// printf("the bestSz is: %d\n", int(KDC.size()));
	if(false && use2hop && useHeu){
		HeuriSearcher* heuri_solver=new HeuriSearcher(n,m,pstart, edges, gamma, maxDeg, KDC);
		heuri_solver->search(KDC); // Your code Here
	}
	// printf("bestSz after the heu search: %d\n",int(KDC.size()));
	// exit(0);
	if(KDC.size() < UB) {		
		int old_size = KDC.size();
		int *out_mapping = new int[n];
		int *rid = new int[n];
		// bool twoHopSearch=false;

		shrink_graph(n, m, seq, core, out_mapping, NULL, rid, pstart, edges);
		// exit(0);
		int *deg = new int[n]; for(int i = 0;i < n;i ++) deg[i] = pstart[i+1] - pstart[i];

		ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
		linear_heap->init(n, n-1, seq, deg);

		pend = new int[n];
		pend_buf = new int[n];
		int *edge_list = new int[m];
		int *edgelist_pointer = new int[m];
		int *tri_cnt = new int[m/2];
		oriented_triangle_counting(n, m, seq, pstart, pend, edges, edgelist_pointer, rid); // edgelist_pointer currently stores triangle_counts
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer, rid);
		int *active_edgelist = new int[m>>1]; int active_edgelist_n = m>>1;
		for(int i = 0;i < (m>>1);i ++) active_edgelist[i] = i;
		for(int i = 0;i < n;i ++) pend[i] = pstart[i+1];
		int *Qe = new int[m>>1];
		char *deleted = new char[m>>1];
		memset(deleted, 0, sizeof(char)*(m>>1));
		char *exists = new char[n];
		memset(exists, 0, sizeof(char)*n);
		int *Qv = new int[n]; int Qv_n = 0;
		
		m -= 2*peeling(n, linear_heap, Qv, Qv_n, -1, Qe, -1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
		printf("*** Core-Truss Shrink: n = %d, m = %d, Density = %.2f%%\n", n-Qv_n, m/2, double(m)/(n-Qv_n)/(n-Qv_n-1)*100);
		// exit(0);
		Timer tt;

		// exit(0);
		// extract_graph(n,m,deg,id)
		
		int *t_degree = new int[n];
		if(use2hop){
			printf("enter two hop search\n");
			int max_n = n - Qv_n;
			//This is the structure for subgraph
			s_degree = new int[max_n];
			s_pstart = new int[max_n+1];
			s_pend = new int[max_n];
			s_edges = new int[m];
			s_peel_sequence = new int[max_n];
			s_core = new int[max_n];
			s_vis = new char[max_n];
			s_heap = new ListLinearHeap(max_n,max_n-1);
			s_edgelist_pointer = new int[m];
			s_tri_cnt = new int[m/2];
			s_edge_list = new int[m];
			s_active_edgelist = new int[m/2];
			s_deleted = new char[m/2];

			vector<pair<int,int> > vp; vp.reserve(m/2);
			// printf("current QC size: %d\n", int(KDC.size()));
			for(int i = 0;i < n&&m&&KDC.size() < UB;i ++) {
				// printf("enter subgraph search: %d\n", i);
				int u, key; 
				if (!useDegen)
				{
					u=seq[i],key=n;
				}else{
					linear_heap->pop_min(u, key);
				}
				
				
				if(false && key < int(KDC.size()-K)) {
					if(deg[u] != 0) {
						Qv[0] = u; Qv_n = 1;
						m -= 2*peeling(n, linear_heap, Qv, Qv_n, KDC.size()-K, Qe, KDC.size()-K-1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
					}
					continue;
				}
				// printf("m: %d\n",m);
				if(!m) break;
				// printf("enter subgraph search-lcy\n");
				int *ids = Qv; int ids_n = 0; int subUB=UB;
				int pre_size0;
				// do{
				// 	pre_size0=KDC.size();
				// 	induceSubgraph(u, ids, ids_n, rid, vp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
				// 	if(ids_n > KDC.size()&& vp.size()*2 < m) subUB=min(UB,subgraph_heuri(ids, ids_n, vp, rid, Qv, Qe, exists));
				// }
				// while(KDC.size()!=pre_size0);
				induceSubgraph(u, ids, ids_n, rid, vp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
				// printf("vp: %d\n", vp.size());
				// exit(0);
				// if(ids_n > KDC.size()&& vp.size()*2 < m) subUB=min(UB,subgraph_heuri(ids, ids_n, vp, rid, Qv, Qe, exists));
				int pre_size = KDC.size();
				// printf("pre_size: %d, ids_n: %d, subUB: %d\n", pre_size, ids_n, subUB);
				if(ids_n > pre_size && subUB> pre_size) { 
#ifndef _TEST_
					printf("enter subgraph %d search\n", i);
#endif
					if(useMIP){
#ifndef _TEST_
						printf("using MIP 2hop\n");
#endif
						MIP *MIPSolver=new MIP();
						MIPSolver->load_subgraph(ids_n, vp, int(KDC.size()),gamma);
						int maxSz=MIPSolver->MIPSolve2hop();
#ifndef _TEST_
						printf("MIP finish computing subgraph %d\n",i);
#endif
						if(maxSz>KDC.size()) KDC.resize(maxSz);
						delete MIPSolver;
					}else{
						// printf("use 2hop bb_search\n");
						QuasiClique_BB *MQCSolver=new QuasiClique_BB();
						subSz.push_back(ids_n);
						// if(core2hop[u]!=ids_n-1 && core2hop[u]==max2hopCore){
						// 	printf("error u: %d, core2hop: %d, subsz: %d\n",u, core2hop[u], ids_n-1);
						// 	exit(0);
						// }
						MQCSolver->load_subgraph(gamma, ids_n, vp, KDC,UB);

						// MQCSolver->printInfo();
						// exit(0);
						// printf("enter bb search\n");
						MQCSolver->MQCSearch2hop(KDC);
						
						// ExactSearcher exact_solver(K, ids_n, vp, KDC, UB); exact_solver.search();
						if(KDC.size() != pre_size) for(int j = 0;j < KDC.size();j ++) KDC[j] = ids[KDC[j]];
						// printf("current QC size: %d\n", int(KDC.size()));
						// exit(0);
						delete MQCSolver;
					}
					// exit(0);
				}
				Qv[0] = u; Qv_n = 1;
				m -= 2*peeling(n, linear_heap, Qv, Qv_n, -1, Qe, -1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
			}
		
		}else{
			if(useMIP){
#ifndef _TEST_
				printf("enter MIP");
#endif
				MIP *MIPSolver= new MIP();
				MIPSolver->load_graph(n, pstart, pstart+1, edges, int(KDC.size()),gamma);
				int maxSz=MIPSolver->MIPSolve();
				if(maxSz>KDC.size()){
					KDC.resize(maxSz);
				}
				delete MIPSolver;
			}else{
				if(UB<0){
					printf("error\n");
				}
				QuasiClique_BB *MQCSolver=new QuasiClique_BB();
				MQCSolver->load_graph(n,pstart,pstart+1,edges);
#ifndef _TEST_
				printf("graph load success!\n");
#endif
				// exit(0);
				// KDC.clear(); old_size=0;
				// if(UB<0){
				// 	printf("error\n");
				// }
				MQCSolver->MQCSearch(gamma, UB, KDC);
				// if(KDC.size()>old_size){
				// 	for (int i = 0; i < KDC.size(); i++) KDC[i]=out_mapping[KDC[i]];
				// }
#ifndef _TEST_
				printf("BBSearch complete, tree count: %lld\n", treeCnt);
#endif
				// exit(0);
				delete MQCSolver;
			}
		}
		if(KDC.size() > old_size) for(int i = 0;i < KDC.size();i ++) KDC[i] = out_mapping[KDC[i]];
		
#ifdef _TEST_
		printf("#MaxQCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("#maxP=%d\n#minPUB=%d\n#maxME=%d\n", max_P_end, P_UBMin,maxME);
		printf("#NodeCount=%lld\n",treeCnt);
		printf("#MaxSG=%d\n",maxSubSz);
		printf("#FeasibleSubgraph=%d\n",feasible);
		printf("#prune1=%lld\n#ubprune=%lld\n", prune1,ub_prune);
		printSubInfo();
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		// printf("Max Mutual Exclusive: %d, Mutual exclusive sum: %lld, Mutual exclusive avg: %lf, Mutual exclusive dense avg: %lf\n",MaxMuExNum, MuExSum, (double)MuExSum/(tree_cnt+0.1), denSum/(denNum+0.001));
		// printf("Sum of Mutual Exclusive in color set: %lld, max mutual exclusive in color set: %d\n", colMuExSum, maxColMuNum);
		printf("Search Tree Size: %lld\n",treeCnt);
		printf("Feasible Subgraph: %d\n",feasible);
#endif
		
		
// 		if(treeCnt>0){
// #ifdef _TEST_
// 			printf("#Bound-Prune=%lld %.2f%%\n#Color-Bound-Prune=%lld %.2f%%\n", boundPrune,double(boundPrune)/(boundPrune+colorBndPrune+tree_cnt)*100, colorBndPrune, double(colorBndPrune)/(boundPrune+colorBndPrune+tree_cnt)*100);
// #else
// 			// printf("Bound Prune: %lld %.2f%%, Color Bound Prune: %lld %.2f%%\n", boundPrune,double(boundPrune)/(boundPrune+colorBndPrune+tree_cnt)*100, colorBndPrune, double(colorBndPrune)/(boundPrune+colorBndPrune+tree_cnt)*100);
// #endif
// 		}
		delete linear_heap;
		delete[] t_degree;
		delete[] deg2hop;
		delete[] core2hop;
		delete[] exists;
		delete[] out_mapping;
		delete[] rid;
		delete[] deg;
		delete[] edgelist_pointer;
		delete[] tri_cnt;
		delete[] active_edgelist;
		delete[] Qe;
		delete[] Qv;
		delete[] deleted;
	}else if(KDC.size()>=UB){
#ifdef _TEST_
		printf("#MaxQCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", KDC.size(), 0.0, double(t.elapsed())/1000000);
		printf("#NodeCount=%lld\n",treeCnt);
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), 0.0, double(t.elapsed())/1000000);
		printf("Search Tree Size: %lld, the sum of induced graphs: %lld\n",0,0);
#endif
	}
	delete[] core;
	delete[] seq;

}

void Graph::load_graph_from_edgelist(int _n, const vector<pair<int,int> > &edge_list, int &n, int &m, int *degree, int *pstart, int *edges) {
	n = _n;
	m = (int)edge_list.size()*2;
	for(int i = 0; i < n; i++) degree[i] = 0;
	for(int i = 0;i < m/2;i ++) {
		assert(edge_list[i].first >= 0&&edge_list[i].first < n&&edge_list[i].second >= 0&&edge_list[i].second < n);
		degree[edge_list[i].first] ++;
		degree[edge_list[i].second] ++;
	}

	pstart[0] = 0;
	for(int i = 0;i < n;i ++) pstart[i+1] = pstart[i]+degree[i];
	for(int i = 0;i < m/2;i ++) {
		int a = edge_list[i].first, b = edge_list[i].second;
		edges[pstart[a]++] = b;
		edges[pstart[b]++] = a;
	}
	for(int i = 0;i < n;i ++) pstart[i] -= degree[i];
}

void Graph::extract_graph(int n, int m, int *degree, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) {
	ids_n = 0; vp.clear();
	for(int i=0; i<n; ++i){
		if(degree[i]){
			ids[ids_n] = i; rid[i] = ids_n++;
		}
	}
	for(int i = 0; i<ids_n ; i++) {
		int u = ids[i];
		for(int j = pstart[u];j < pend[u]; j++) if(!deleted[edgelist_pointer[j]] && u < edges[j]) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
}

void Graph::extract_subgraph(int u, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) {
	ids_n = 0; vp.clear();
	ids[ids_n++] = u; exists[u] = 1; rid[u] = 0;
	int u_n = pstart[u];
	for(int i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		int v = edges[i];
		rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
	}
	pend[u] = u_n;
	int old_size = ids_n;
	for(int i = 1;i < old_size;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			int v = edges[j];
			if(exists[v]) continue;
			rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < old_size;i ++) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) if(edges[j] > u) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(int i = old_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
	
}

	
void Graph::induceSubgraph(int u, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, int *Q, int* degree, char *exists, int *pend, char *deleted, int *edgelist_pointer) {
	vp.clear();
	ids_n = 0; ids[ids_n++] = u; exists[u] = 1;
	int u_n = pstart[u];
	for(int i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		int v = edges[i];
		ids[ids_n++] = v; exists[v] = 2;
	}
	pend[u] = u_n;
	
	int Q_n = 0;
	for(int i = 1;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(exists[edges[j]] == 2) ++ degree[u];
		}
		pend[u] = u_n;
		if(degree[u]+2+K <= KDC.size()) Q[Q_n++] = u;
	}
	for(int i = 0;i < Q_n;i ++) {
		u = Q[i];
		exists[u] = 10;
		for(int j = pstart[u];j < pend[u];j ++) if(exists[edges[j]] == 2) {
			if( (--degree[edges[j]])+2+K == KDC.size()) {
				assert(Q_n < m/2);
				Q[Q_n++] = edges[j];
			}
		}
	}
	assert(Q_n <= ids_n);
	if(ids_n-Q_n+K<=KDC.size()) {
		for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
		ids_n = 0;
		return ;
	}
	
	int nr_size = ids_n;
	for(int i = 1;i < nr_size;i ++) if(exists[ids[i]] == 2) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) {
			if(!exists[edges[j]]) {
				ids[ids_n++] = edges[j];
				exists[edges[j]] = 3;
				degree[edges[j]] = 1;
			}
			else if(exists[edges[j]] == 3) ++ degree[edges[j]];
		}
	}

#ifndef NDEBUG
	//printf("Entire list: ");
	//for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
	//printf("\n");
#endif

	int new_size = 1;
	for(int i = 1;i < nr_size;i ++) {
		if(exists[ids[i]] == 10) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if(new_size + Q_n != nr_size) {
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for(int i = 0;i < new_size;i ++) printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for(int i = 0;i < Q_n;i ++) printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	int old_nr_size = nr_size;
	nr_size = new_size;
	for(int i = old_nr_size;i < ids_n;i ++) {
		if(degree[ids[i]]+2+K-1<= KDC.size()) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for(int i = 1;i < nr_size;i ++) assert(exists[ids[i]] == 2);
	for(int i = nr_size;i < ids_n;i ++) assert(exists[ids[i]] == 3);
#endif

	//for(ui i = 0;i < ids_n;i ++) printf(" %u", ids[i]);
	//printf("\n");

	for(int i = 0;i < ids_n;i ++) {
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for(int i = 0;i < nr_size;i ++) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) if(exists[edges[j]]&&edges[j] > u) {
			assert(!deleted[edgelist_pointer[j]]);
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(int i = nr_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
#ifndef NDEBUG
	for(int i = 0;i < n;i ++) assert(exists[i] == 0);
#endif
}
void Graph::printSubInfo(){
	vector<int> subSzNum;
	int subtot=(int)subSz.size();
	double subAvg=0.0;
	double subMid=0.0;
	int subMost=0, subMostNum=0;
	int subMax=maxSubSz;
	double subMostRatio=0.0;
	
	subSzNum.resize(subMax+1);
	if(subtot<=0){
		printf("#subtot=%d\n",subtot);
		printf("#subMax=NA\n");
		printf("#subAvg=NA\n");
		printf("#subMid=NA\n");
		printf("#subMost=NA\n");
		printf("#subMostRatio=NA\n");
	}else{
		sort(subSz.begin(),subSz.end());
		for (int i = 0; i < subtot; i++){
			// subMax=max(subMax, subSz[i]);//update the subMax
			int subsz=subSz[i];
			// printf("ori----i: %d, subAvg: %.2f\n",i,subAvg);
			double a=double(i)/(double(i+1)), b=double(subsz)/(double(i+1));
			// subAvg= subAvg*(double(i)/(i+1))+double(subsz)/(i+1);//build the average
			// printf("i: %d, a: %.lf, b: %.lf\n", i,a,b);
			subAvg=subAvg*a+b;
			// printf("i: %d, subsz: %d, subAvg: %.2f\n",i, subsz,subAvg);
			// int cur_subSz=0, cur_subNum=0;
			
			subSzNum[subsz]++;
		}
		//compute the medium number of the size
		if((subtot-1)%2==0){
			int midId=subtot/2;
			subMid=subSz[midId];
		}else{
			int midId1=int(floor(double(subtot-1)/2.0));
			int midId2=int(ceil(double(subtot-1)/2.0));
			subMid=double(subSz[midId1])/2.0+double(subSz[midId2])/2.0;//update the mid number
		}
		//compute the most number of the size
		for (int sz = 0; sz <= subMax; sz++){
			int cursz=0, curnum=0;
			if(subSzNum[sz] < subMostNum) continue; 
			subMost=sz, subMostNum=subSzNum[sz];
		}
		subMostRatio=(double)subMostNum/subtot;
		printf("#subtot=%d\n",subtot);
		printf("#subMax=%d\n",subMax);
		printf("#subAvg=%.2f\n", subAvg);
		printf("#subMid=%.2f\n", subMid);
		printf("#subMost=%d\n", subMost);
		printf("#subMostRatio=%.2f\n", subMostRatio);
	}
}
// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
int Graph::degen(int n, int *seq, int *core, int *pstart, int *edges, int *degree, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;
	int threshold = KDC.size()-K;
	int edgeCnt=0;
	for(int i = 0;i < n;i++) degree[i] = pstart[i+1] - pstart[i], edgeCnt+=degree[i]; edgeCnt=edgeCnt/2;
	int queue_n = 0, new_size = 0;
	for(int i = 0;i < n;i++) if(degree[i] < threshold) seq[queue_n ++] = i;//first order reduction
	for(int i = 0;i < queue_n;i ++) {
		int u = seq[i]; degree[u] = 0;
		for(int j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) seq[queue_n ++] = edges[j];
			edgeCnt--;
		}
	}
	int UB = n; if(queue_n == n) UB = KDC.size();
	memset(vis, 0, sizeof(char)*n);
	for(int i = 0;i < n;i ++) {
		if(degree[i] >= threshold) seq[queue_n + (new_size ++)] = i;
		else vis[i] = 1, core[i] = 0;
	}
	assert(queue_n + new_size == n);

	if(new_size != 0) {
		heap->init(new_size, new_size-1, seq+queue_n, degree);
		int maxCore = 0; int idx = n;
		int t_UB=0;
		for(int i = 0;i < new_size;i ++) {
			if(idx == n && edgeCnt >= ceil(gamma*(1.0*(new_size - i)*(new_size - i - 1)/2.0)) ) idx = i;
			int u, key; heap->pop_min(u, key);
			if(key > maxCore) maxCore = key; core[u] = maxCore;
			seq[queue_n + i] = u; vis[u] = 1;
			t_UB=max(t_UB,min(core[u]+K+1,new_size-i));
			for(int j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) {
				heap->decrement(edges[j], 1);
			}
			edgeCnt-=key;
		}
		// UB=min(UB,min(t_UB,maxCore+(int)floor((1+sqrt(1+8*K))/2)));
		UB=n;
		// cout << "new size: "<<new_size<<endl<<"idx: "<<idx<<endl;
		if(new_size - idx > KDC.size()) { 
			// cout<<"enter\n";
			KDC.clear(); for(int i = idx;i < new_size;i ++) KDC.pb(seq[queue_n + i]);
			if(!output) printf("Find a QC of size: %u\n", new_size - idx);
		}
#ifdef _TEST_
		if(output) printf("#HeuSize=%u\n#MaxCore=%d\n#UB=%d\n#HeuTime=%.2f\n", KDC.size(), maxCore, UB, double(t.elapsed())/1000000);
		printf("#MaxDeg=%d\n",maxDeg);
#else
		if(output) printf("*** HeuriQDC size: %u, MaxCore: %d, UB: %d, Heuri Time: %.2f\n", (unsigned int)KDC.size(), maxCore, UB, double(t.elapsed())/1000000);
#endif
	}
	return UB;
}
int Graph::degen1hop(int n, int *seq, int *core, int *pstart, int *edges, int *degree, char *vis, ListLinearHeap *heap, bool output){
	Timer t;
	int UB=n;
	for (int i = 0; i < n; i++) seq[i]=i;
	for (int i = 0; i < n; i++) degree[i]=pstart[i+1]-pstart[i];
	
	memset(vis,0,n*sizeof(char));
	if(n>0){
		heap->init(n,n-1,seq,degree);
		maxCore=0;
		for (int i = 0; i < n; i++){
			int u,key; heap->pop_min(u,key);
			if(key>maxCore) maxCore=key; core[u]=maxCore;
			// printf("maxCore i:%d, %d\n",i, maxCore); 
			seq[i]=u, vis[u]=1;
			for (int j = pstart[u]; j < pstart[u+1]; j++){
				int v=edges[j];
				if (vis[v]==0) heap->decrement(v,1);
			}
		}
	}
	// for (int i = 0; i < n; i++)
	// {
	// 	printf("core[%d]: %d\n",i,core[i]);
	// }
	
#ifdef _TEST_
		if(output) printf("#MaxCore=%d\n#DegenTime=%.2f\n", maxCore, double(t.elapsed())/1000000);
		printf("#MaxDeg=%d\n",maxDeg);
#else
		if(output) printf("*** MaxCore: %d, UB: %d, Degen Time: %.2f\n", maxCore, UB, double(t.elapsed())/1000000);
#endif
	return n;
}
int Graph::degenBy2hopDeg(int n, int *seq, int *core2hop, int *pstart2hop, int *adj2hop, int *deg2hop, char *vis, ListLinearHeap *heap, bool output){
	// printf("enter dengeracy ordering\n");
	Timer t;
	int UB=n;
	int del_n=0;
	for (int i = 0; i < n; i++) seq[i]=i;
	memset(vis,0,n*sizeof(char));
	if(n>0){
		heap->init(n,n-1,seq,deg2hop);
		max2hopCore=0;
		for (int i = 0; i < n; i++){
			int u, key; heap->pop_min(u,key);
			if(key > max2hopCore) max2hopCore=key; core2hop[u]=max2hopCore;
			seq[i]=u; vis[u]=1;//vis should be false at start
			for (int j = pstart2hop[u]; j < pstart2hop[u+1]; j++){
				if(vis[adj2hop[j]]==0) heap->decrement(adj2hop[j],1);
			}
		}
	}
	UB=max2hopCore+1;
#ifdef _TEST_
	if(output) printf("#MaxCore2hop=%d\n#UB=%d\n#Order2hopTime=%.2f\n", max2hopCore, UB, double(t.elapsed())/1000000);
#elif
	if(output) printf("MaxCore2hop: %d, UB: %d, Order Time: %.2f\n", max2hopCore, UB, double(t.elapsed())/1000000);
#endif
	return UB;
	// printf("complete dengeracy ordering\n");
}
int Graph::degen2hop4Large(int n, int *seq, int *core2hop, int *pstart, int* edges, int *deg2hop, char *vis, ListLinearHeap *heap, bool output){
	Timer t;
	int UB=n;
	int del_n=0;
	int max2hop=0, min2hop=n;
	bool *visNei=new bool[n];
	memset(visNei, false, n*sizeof(bool));
	memset(vis, 0, n*sizeof(char));
	for (int i = 0; i < n; i++) seq[i]=i;//init sequence
	
	vector<int> nei2hop;
	//1. construct the deg2hop
	// printf("enter construction\n");
	for (int u = 0; u < n; u++){
		visNei[u]=true;
		for (int j = pstart[u]; j < pstart[u+1]; j++){
			int v=edges[j];
			if(visNei[v]) continue;
			nei2hop.push_back(v);
			visNei[v]=true;
		}
		int curDeg=pstart[u+1]-pstart[u];
		for (int j = 0; j < curDeg; j++){
			int v=nei2hop[j];
			//find the neighbors of u
			for (int k = pstart[v]; k < pstart[v+1]; k++){
				int w=edges[k];
				if(visNei[w]) continue;
				nei2hop.push_back(w);
				visNei[w]=true;
			}
		}
		deg2hop[u]=nei2hop.size();
		min2hop=min(min2hop, deg2hop[u]), max2hop=max(max2hop, deg2hop[u]);
		for(auto v: nei2hop) visNei[v]=false;
		visNei[u]=false;
		nei2hop.clear();
	}
	// for (int i = 0; i < n; i++) printf("deg2hop: %d\n",deg2hop[i]);
#ifdef _TEST_
	printf("#max2hop=%d\n#min2hop=%d\n#2hopbuildtime=%.2f\n", max2hop, min2hop,double(t.elapsed())/1000000);
#elif
	printf("The max 2hop is: %d, the min 2hop is: %d, the time of 2hop construction is: %.2f\n", max2hop, min2hop, double(t.elapsed())/1000000);
#endif
	
	// printf("out construction\n");
	//2. reorder by 2-hop degenracy
	// printf("enter degeneracy\n");
	if(n>0){
		heap->init(n,n-1,seq+del_n,deg2hop);
		max2hopCore=0;
		for (int i = 0; i < n; i++){
			int u, key; heap->pop_min(u,key);
			if(key > max2hopCore) max2hopCore=key; core2hop[u]=max2hopCore;
			seq[i]=u; vis[u]=1;//vis should be false at start
			// for (int j = pstart2hop[u]; j < pstart2hop[u+1]; j++) heap->decrement(adj2hop[j],1);
			//find the 2hop 
			visNei[u]=true;
			for (int j = pstart[u]; j < pstart[u+1]; j++){
				int v=edges[j];
				// if v is in nei2hop or v is deleted
				if(visNei[v]) continue;
				nei2hop.push_back(v);
				visNei[v]=true;
			}
			int curDeg=nei2hop.size();
			for (int i = 0; i < curDeg; i++){
				int v=nei2hop[i];
				for (int j = pstart[v]; j < pstart[v+1]; j++){
					int w=edges[j];
					if(visNei[w]) continue;
					nei2hop.push_back(w);
					visNei[w]=true;
				}
			}

			for(auto v: nei2hop) {
				if(vis[v]==0) heap->decrement(v,1);
				visNei[v]=false;
			}
			visNei[u]=false;
			nei2hop.clear();
		}
	}
	// printf("out degeneracy\n");

	UB=max2hopCore+1;
#ifdef _TEST_
	if(output) printf("#MaxCore2hop=%d\n#UB=%d\n#Order2hopTime=%.2f\n", max2hopCore, UB, double(t.elapsed())/1000000);
#elif
	if(output) printf("MaxCore2hop: %d, UB: %d, Order Time: %.2f\n", max2hopCore, UB, double(t.elapsed())/1000000);
#endif
	return UB;
}



// in_mapping and out_mapping can be the same array
// note that core is not maintained, and is assumed to not be used anymore
void Graph::shrink_graph(int &n, int &m, int *peel_sequence, int *core, int *out_mapping, int *in_mapping, int *rid, int *pstart, int *edges) {
	int cnt = 0;
	for(int i = 0;i < n;i ++) if(core[i] + K + 1> KDC.size()) {
		rid[i] = cnt;
		if(in_mapping == NULL) out_mapping[cnt] = i;
		else out_mapping[cnt] = in_mapping[i];
		++ cnt;
	}

	if(cnt != n) {
		cnt = 0;
		int pos = 0;
		for(int i = 0;i < n;i ++) if(core[i] + K + 1 > KDC.size()) {
			int t_start = pstart[i]; pstart[cnt] = pos;
			for(int j = t_start;j < pstart[i+1];j ++) if(core[edges[j]] + K + 1> KDC.size()) {
				edges[pos ++] = rid[edges[j]];
			}
			++ cnt;
		}
		pstart[cnt] = pos;

		//printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n-cnt-1]] == 0||core[peel_sequence[n-cnt-1]] + K + 1<= KDC.size());
		assert(cnt == 0||core[peel_sequence[n-cnt]] + K + 1> KDC.size());
		for(int i = 0;i < cnt;i ++) {
			peel_sequence[i] = rid[peel_sequence[n-cnt+i]];
			//core[i] = core[out_mapping[i]];
		}

		n = cnt;
		m = pos;
	}

	printf("*** Core Shrink: n = %d, m = %d \n", n, m/2);
}

// orient graph and triangle counting
void Graph::oriented_triangle_counting(int n, int m, int *peel_sequence, int *pstart, int *pend, int *edges, int *tri_cnt, int *adj) {
	int *rid = adj;
	for(int i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(int i = 0;i < n;i ++) {
		int &end = pend[i] = pstart[i];
		for(int j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for(int i = 0;i < n;i ++) sum += pend[i] - pstart[i];
	// printf("%lld %lld\n", sum, m);
	assert(sum*2 == m);
#endif

	memset(adj, 0, sizeof(int)*n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(int)*m);
	for(int u = 0;u < n;u ++) {
		for(int j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(int j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];
			for(int k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
				++ tri_cnt[j];
				++ tri_cnt[k];
				++ tri_cnt[adj[edges[k]]-1];
				++ cnt;
			}
		}

		for(int j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(int n, int *tri_cnt, int *edge_list, int *pstart, int *pend, int *pend2, int *edges, int *edgelist_pointer, int *buf) {
	for(int i = 0;i < n;i ++) pend2[i] = pend[i];
	int pos = 0;
	for(int i = 0;i < n;i ++) {
		for(int j = pstart[i];j < pend[i];j ++) {
			tri_cnt[pos>>1] = edgelist_pointer[j]; edge_list[pos++] = i; edge_list[pos++] = edges[j];

			int &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos>>1)-1;
			edges[k ++] = i;
		}
	}

#ifndef NDEBUG
	for(int i = 0;i < n;i ++) assert(pend2[i] == pstart[i+1]);
#endif

	for(int i = 0;i < n;i ++) {
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for(int i = 0;i < n;i ++) {
		for(int j = pend2[i];j < pstart[i+1];j ++) {
			int &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k ++] = i;
		}
	}

	int *ids = pend2;
	for(int i = 0;i < n;i ++) {
		if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
		int j = pstart[i], k = pend[i], pos = 0;
		while(j < pend[i]&&k < pstart[i+1]) {
			if(edges[j] < edges[k]) {
				ids[pos] = edges[j];
				buf[pos ++] = edgelist_pointer[j ++];
			}
			else {
				ids[pos] = edges[k];
				buf[pos ++] = edgelist_pointer[k ++];
			}
		}
		while(j < pend[i]) {
			ids[pos] = edges[j];
			buf[pos ++] = edgelist_pointer[j ++];
		}
		while(k < pstart[i+1]) {
			ids[pos] = edges[k];
			buf[pos ++] = edgelist_pointer[k ++];
		}
		for(int j = 0;j < pos;j ++) {
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}


char Graph::find(int u, int w, int &b, int e, char *deleted, int &idx, int *edgelist_pointer, int *edges) {
	if(b >= e) return 0;

	while(b+1 < e) {
		idx = b + (e-b)/2;
		if(edges[idx] > w) e = idx;
		else b = idx;
	}

	if(edges[b] == w) {
		idx = edgelist_pointer[b];
		if(!deleted[idx]) return 1;
	}

	return 0;
}

// return the number of peeled edges
int Graph:: peeling(int critical_vertex, ListLinearHeap *linear_heap, int *Qv, int &Qv_n, int d_threshold, int *Qe, int t_threshold, int *tri_cnt, int *active_edgelist, int &active_edgelist_n, int *edge_list, int *edgelist_pointer, char *deleted, int *degree, int *pstart, int *pend, int *edges, char *exists) {
	int Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	bool initialize_Qe=t_threshold>0;
	if(initialize_Qe) {
		int active_edgelist_newn = 0;
		for(int j = 0;j < active_edgelist_n;j ++) if(!deleted[active_edgelist[j]]) {
			if(tri_cnt[active_edgelist[j]] < t_threshold) Qe[Qe_n++] = active_edgelist[j];
			else active_edgelist[active_edgelist_newn ++] = active_edgelist[j];
		}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	//printf("%lu\n", Qe_n);

	int deleted_edges_n = 0;
	int Qv_idx = 0;
	while(Qv_idx < Qv_n || Qe_n) {
		if(Qe_n == 0) {
			//printf("hit\n");
			int u = Qv[Qv_idx ++]; // delete u from the graph due to have a degree < d_threshold
			int u_n = pstart[u];
			for(int k = pstart[u];k < pend[u];k ++) if(!deleted[edgelist_pointer[k]]) {
				edges[u_n] = edges[k]; edgelist_pointer[u_n++] = edgelist_pointer[k];
				exists[edges[k]] = 1;
			}
			pend[u] = u_n;

			for(int k = pstart[u];k < pend[u];k ++) deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			degree[u] = 0;
			if(linear_heap != NULL) linear_heap->del(u);
			//printf("Removed %u\n", u);

			for(int k= pstart[u];k < pend[u];k ++) {
				int v = edges[k];
				int v_n = pstart[v];
				for(int x = pstart[v];x < pend[v];x ++) if(!deleted[edgelist_pointer[x]]) {
					edges[v_n] = edges[x]; edgelist_pointer[v_n++] = edgelist_pointer[x];
					if(edges[x] > v&&exists[edges[x]]) {
						if( (tri_cnt[edgelist_pointer[x]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[x];
					}
				}
				pend[v] = v_n;

				if( (degree[v]--) == d_threshold) {
					Qv[Qv_n++] = v;
					if(v == critical_vertex) {
						for(int k = pstart[u];k < pend[u];k ++) exists[edges[k]] = 0;
						return 0;
					}
				}
				if(linear_heap != NULL) linear_heap->decrement(v, 1);
			}

			for(int k = pstart[u];k < pend[u];k ++) exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for(int j = 0;j < Qe_n;j ++) {
			int idx = Qe[j];
			int u = edge_list[idx<<1], v = edge_list[(idx<<1)+1];
			int tri_n = tri_cnt[idx];
			//printf("remove %u %u\n", u, v);
			deleted[idx] = 1;
			if( (degree[u] --) == d_threshold) {
				Qv[Qv_n++] = u;
				if(u == critical_vertex) return 0;
			}
			if( (degree[v] --) == d_threshold) {
				Qv[Qv_n++] = v;
				if(v == critical_vertex) return 0;
			}
			//printf("before\n");
			if(linear_heap != NULL) {
				linear_heap->decrement(u, 1);
				linear_heap->decrement(v, 1);
			}
			//printf("after\n");
			deleted_edges_n ++;
			
			if(degree[u] < degree[v]) swap(u,v);
			//printf("here\n");

			if(degree[u] > degree[v]*2) { // binary search
			//if(false) {
				int v_n = pstart[v], start = pstart[u];
				for(int k = pstart[v];k < pend[v];k ++) if(!deleted[edgelist_pointer[k]]) {
					edges[v_n] = edges[k]; edgelist_pointer[v_n++] = edgelist_pointer[k];

					if(tri_n&&find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges)) {
						-- tri_n;
						if( (tri_cnt[idx]--) == t_threshold) Qe[Qe_n++] = idx;
						if( (tri_cnt[edgelist_pointer[k]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[k];
					}
				}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else { // sorted_merge
				int ii = pstart[u], jj = pstart[v];
				int u_n = pstart[u], v_n = pstart[v];

				while(true) {
					while(ii < pend[u]&&deleted[edgelist_pointer[ii]]) ++ ii;
					while(jj < pend[v]&&deleted[edgelist_pointer[jj]]) ++ jj;
					if(ii >= pend[u]||jj >= pend[v]) break;

					if(edges[ii] == edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if( (tri_cnt[edgelist_pointer[ii]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[ii];
						if( (tri_cnt[edgelist_pointer[jj]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[jj];

						++ ii;
						++ jj;
					}
					else if(edges[ii] < edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ ii;
					}
					else {
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++ jj;
					}
				}
				while(ii < pend[u]) {
					if(!deleted[edgelist_pointer[ii]]) {
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ ii;
				}
				while(jj < pend[v]) {
					if(!deleted[edgelist_pointer[jj]]) {
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++ jj;
				}
				pend[u] = u_n; pend[v] = v_n;
			}
		}
		Qe_n = 0;
	}
	return deleted_edges_n;
}


int Graph::subgraph_heuri(ui *ids, ui &_n, vector<pair<int,int> > &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists) {
	ui s_n; ept s_m;
	load_graph_from_edgelist(_n, edge_list, s_n, s_m, s_degree, s_pstart, s_edges);
	return degen(s_n, s_peel_sequence, s_core, s_pstart, s_edges, s_degree, s_vis, s_heap,false);
}

void Graph::print_info(){
	double density=m/(1.0*n*(n-1));
	// printf("#graph=%s\n#gamma=%lf\n",dir, gamma);
	printf("graph=%s\n",this->dir.c_str());
	printf("n=%d, m=%d, max degree=%d, dense=%.3f\n", this->n, this->m, maxDeg, density);
	printf("gamma=%Lf\n ", gamma);
	cout<<"miss edges upper bound: "<<K<<endl;
}

void Graph::setK(){
	K=floor((1-gamma)*n*(n-1)/2.0);
}