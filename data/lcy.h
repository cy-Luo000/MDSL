// #pragma once
// #define MAX_INDUCE_SIZE 10000
// #include<bitset>
// #include"LinearHeap.h"
// #include"MBitSet.h"
// using namespace std;
// //!warning: lcy's code
// typedef bitset<MAX_INDUCE_SIZE> Vset;
// Vset P;//the v in P
// Vset C;//the v in C
// Vset NP;//the v not in P
// Vset NC;//the v not in C
// Vset bitMtx[MAX_INDUCE_SIZE];// the adjacent matrix
// Vset rbitMtx[MAX_INDUCE_SIZE];// the non-adjacent matrix
// Vset bitcolor[MAX_INDUCE_SIZE];//to record the vertices of each color
// int tree_cnt=0;
// typedef struct{
// 	MBitSet64 bitcolor;
// 	int sz;
// }ColorType;
// //!warning: lcy's code

// class HeuriSearcher{
// public:
// 	int n;
// 	int m;
// 	int* pstart;
// 	int* edges;
//     int k; //total missing edges
//     int LB; // the lowerbound
//     int UB; // the upperbound
// 	// vector<int> &KDC;

//     // int* PC;
//     // int* PC_rid;// to record the positon in PC of each vertex v of G
//     // int P_end;// the P is from PC[0] to PC[P_end-1]
//     // int C_end;// the C is from PC[P_end] to PC[C_end-1]

//     // int* neiInP;// to record every vertex's degree in P
//     // int* neiInG;// to record every vertex's degree in C
//     // int MEInP;// the missing edges in P 
//     // int MEInG;// the missing edges in C

// 	// Your code
// 	void search(){}
// };

// //exact_solver(K, ids_n, vp, KDC, UB)
// class ExactSearcher{
// public:
// 	int n;
// 	int m;
// 	int* pstart;
// 	int* edges;
//     int k; //total missing edges
//     int LB; // the lowerbound
//     int UB; // the upperbound
// 	vector<int> &KDC;

//     int* PC;
//     int* PC_rid;// to record the positon in PC of each vertex v of G
//     int P_end;// the P is from PC[0] to PC[P_end-1]
//     int C_end;// the C is from PC[P_end] to PC[C_end-1]

//     int* neiInP;// to record every vertex's degree in P
//     int* neiInG;// to record every vertex's degree in C
// 	bool* matrix;
//     int MEInP;// the missing edges in P 
//     int MEInG;// the missing edges in C

// 	int max_CC;//! lcy's code
// 	int cur_B;
// 	//!warning: lcy's code
// 	int *neiInC;// to record every vertex's degree in C
// 	MBitSet64 **rows;
// 	ColorType *colors;
// 	ListLinearHeap *C_heap_P;// to record the order of vertices based on the num of non-neighbors in P
// 	ListLinearHeap *C_heap_C;// to record the order of vertices based on the num of non-neighbors in C
// 	int *color_lable;//to record the color lable of each vertex
// 	int *color_ME;// to record the allowed missing edges in each colorset(independent set)
// 	int *color_num;// the vertices can be added into C in each color
// 	// to record the reverse graph of the induced graph
// 	int* rpstart;
// 	int* redges;
// 	// to record the reverse graph of the induced graph
// 	//!warning: lcy's code
	
// 	bool& isAdj(int x, int y){
// 		return matrix[x*n+y];
// 	}
// 	bool notAdj(int x,int y){
// 		return rbitMtx[x].test(y);
// 	}
//    ExactSearcher(int _k, int _n, vector<pair<int,int>> &_vp, vector<int> &_KDC, int _UB)
//    :KDC(_KDC)
//    {
// 		n=_n;
// 		m=_vp.size();
// 		k=_k;
// 		LB=_KDC.size();
// 		UB=_UB;

// 		pstart=new int[n+1];
// 		edges=new int[m*2];
//         PC=new int[n];
//         PC_rid=new int[n];
//         neiInP=new int[n];
//         neiInG=new int[n];
// 		matrix=new bool[n*n];
// 		//!warning: lcy's code
// 		rows=new MBitSet64*[n];
// 		for (int i = 0; i < n; i++)
// 			rows[i]=new MBitSet64(n);
// 		colors=new ColorType[n];
// 		for (int i = 0; i < n; i++)
// 		{
// 			colors[i].bitcolor=MBitSet64(n),colors[i].sz=0;
// 		}
		
// 		max_CC=0;//! lcy's code
// 		neiInC=new int[n];
// 		color_lable=new int[n];
// 		color_ME=new int[n+1];// here the size of array color_ME should be the color cnt of C, but we don't know color cnt before coloring
// 		color_num=new int[n+1];
// 		rpstart=new int[n+1];
// 		redges=new int[n*n];
// 		C_heap_P=new ListLinearHeap(n,n-1);
// 		C_heap_C=new ListLinearHeap(n,n-1);
// 		//!warning: lcy's code
// 		// KDC.clear();//! lcy's check
// 		P_end=0; C_end=n;
// 		memset(matrix,0,sizeof(bool)*n*n);
// 		memset(neiInP,0,sizeof(int)*n);
		
// 		//!warning: lcy's code
// 		// NP_heap->init(n,n-1,PC,neiInP);
// 		// NC_heap->init(n,n-1,PC,neiInP);
// 		memset(neiInC,0,sizeof(int)*n);
// 		memset(color_ME,0,sizeof(int)*n);
// 		memset(color_num,0,sizeof(int)*n);
// 		P.reset();//at the beginning ,each vertex is not in P
// 		NP.set();
// 		C.set();// at the beginning, each vertex is in C
// 		NC.reset();
// 		for (int i = 0; i < MAX_INDUCE_SIZE; i++)
// 		{
// 			bitMtx[i].reset();//the adjacent matrix of orginal graph is all set 0
// 			rbitMtx[i].set();
// 			rbitMtx[i].reset(i);// the adjacent matrix of complement graph is all set to 1 except the diagonal elements
// 		}
// 		//!warning: lcy's code

// 		for(auto pr : _vp){
// 			int u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
// 			//!warning: lcy's code
// 			rows[u]->set(v);
// 			rows[v]->set(u);
// 			bitMtx[u].set(v);
// 			bitMtx[v].set(u);
// 			rbitMtx[u].reset(v);
// 			rbitMtx[v].reset(u);
// 			//!warning: lcy's code
// 		}
// 		int idx=0,idx_r=0;//!warning: lcy's code
// 		for(int i=0;i<n;++i){
// 			pstart[i]=idx;	
// 			rpstart[i]=idx_r;//!warning: lcy's code
// 			for(int j=0;j<n;++j){
// 				if(bitMtx[i].test(j)) edges[idx++]=j;  //!warning: lcy's code
// 				if(rbitMtx[i].test(j)) redges[idx_r++]=j; //!warning: lcy's code
// 			}
// 			neiInG[i]=idx-pstart[i];
// 			neiInC[i]=neiInG[i];//!warning: lcy's code
// 			// assert(non_neiInC[i]<=n-1);
// 			PC[i]=i;
// 			PC_rid[i]=i;
// 		}
// 		pstart[n]=idx;
// 		rpstart[n]=idx_r;//!warning: lcy's code
// 		MEInP=0;
// 		MEInG=n*(n-1)/2-m;
// 		//!warning: lcy's code
// 		// printf("checkNP_heapinit\n");
// 		C_heap_P->init(n,n-1,PC,neiInP);
// 		// printf("checkNC_heapinit\n");
// 		C_heap_C->init(n,n-1,PC,neiInC);
// 		if (!checkinit())
// 		{
// 			printf("error init\n");
// 		}
// 		checkinit2();
// 		//!warning: lcy's code

//    }

//    ~ExactSearcher(){
//         delete[] PC;
//         delete[] PC_rid;
//         delete[] neiInP;
//         delete[] neiInG;
// 		delete[] matrix;
// 		delete[] pstart;
// 		delete[] edges;

// 		//!warning: lcy's code
// 		delete[] neiInC;
// 		delete[] color_lable;
// 		delete[] color_ME;
// 		delete[] color_num;
// 		delete[] rpstart;
// 		delete[] redges;
// 		delete[] rows;
// 		delete[] colors;
// 		if(C_heap_C!=NULL){
// 			delete C_heap_C;
// 			C_heap_C=NULL;
// 		}
// 		if (C_heap_P!=NULL){
// 			delete C_heap_P;
// 			C_heap_P=NULL;
// 		}
		
// 		//!warning: lcy's code
//    }

//    bool checkinit(){
// 		for(int i=0;i<n;i++){
// 			for (int j = 0; j < n; j++){
// 				if (rows[i]->test(j)!=bitMtx[i].test(j)){
// 					return false;
// 				}
				
// 			}
// 		}
// 		return true;
//    }
//    void checkinit2(){
// 		for(int i=0;i<n;i++){
// 			for (int j = 0; j < n; j++){
// 				assert(rows[i]->test(j)==bitMtx[i].test(j));
// 			}
// 		}
//    }
// 	//!warning: lcy's code
//    inline int intsecColor(int u, int c) {
// 		return colors[c].bitcolor.intersect(*rows[u]);
// 	}
// 	//!warning: lcy's code
//     void checklcy(){
// 		assert(C_heap_C->size==C_end-P_end);
// 		assert(C_heap_P->size==C_end-P_end);
// 		for (int i = 0; i < n; i++)
// 		{
// 			assert((C.test(PC[i])&P.test(PC[i]))==0);
// 		}
		
// 		for (int i = 0; i < C_end; i++)
// 		{
// 			int ele=PC[i];
// 			assert(neiInG[ele]==neiInC[ele]+neiInP[ele]);
// 		}
// 		for (int i = P_end; i < C_end; i++)
// 		{
// 			int ele=PC[i];
// 			assert(neiInP[ele]==C_heap_P->key_s[ele]);
// 			assert(neiInC[ele]==C_heap_C->key_s[ele]);
// 		}
		
// 	}
// 	int checkAbleAdd(){
// 		for (int i = P_end; i < C_end; i++)
// 		{
// 			if (k-MEInP>=(P_end-neiInP[PC[i]]))
// 			{
// 				return PC[i];
// 			}
// 		}
// 		return -1;
// 	}
// 	bool inC(int id){
// 		return id>=P_end&&id<C_end;
// 	}

// 	void C2P(int id){
// 		// checklcy();
// 		assert(inC(id));
// 		int u = PC[id];
// 		swapID(id,P_end++);
// 		//maintain data structure
// 		int nonadj_inP=P_end-1-neiInP[u];
// 		MEInP+=nonadj_inP;
// 		//!warning: lcy's code
// 		P.set(u);C.reset(u);
// 		NP.reset(u);NC.set(u);
// 		C_heap_C->del(u);
// 		C_heap_P->del(u);
// 		for (int i = pstart[u]; i < pstart[u+1]; i++){
// 			neiInP[edges[i]]++,neiInC[edges[i]]=neiInG[edges[i]]-neiInP[edges[i]];
// 			if (C.test(edges[i])){	//if edges[i] in C
// 				C_heap_P->increment(edges[i],1);
// 				C_heap_C->decrement(edges[i],1); //if edges[i] in C,then alter the heap
// 			}
// 		}
// 		// checklcy();
// 		//!warning: lcy's code
		
// 	}

// 	void P2C(int push_cnt){
// 		while (push_cnt--)
// 		{	
// 			// checklcy();
// 			P_end--;
// 			int u=PC[P_end];
// 			int nonadj_inP=P_end-neiInP[PC[P_end]];
// 			MEInP-=nonadj_inP;
// 			//!warning: lcy's code
// 			P.reset(u);C.set(u);
// 			NP.set(u);NC.reset(u);
			
// 			C_heap_P->insert(u,neiInP[u]);
// 			C_heap_C->insert(u,neiInC[u]);
// 			for (int i = pstart[PC[P_end]]; i < pstart[PC[P_end]+1]; i++){
// 				neiInP[edges[i]]--,neiInC[edges[i]]=neiInG[edges[i]]-neiInP[edges[i]];
// 				if (C.test(edges[i])){  //if edges[i] in C
// 					C_heap_C->increment(edges[i],1);
// 					C_heap_P->decrement(edges[i],1);
// 				}
// 			}
// 			// checklcy();
// 			//!warning: lcy's code
// 		}
// 	}

// 	void P2X(){
// 		// checklcy();
// 		int u=PC[P_end-1];
// 		swapID(--P_end,--C_end);
// 		int nonadj_inP=P_end-neiInP[u];
// 		int nonadj_inG=C_end-neiInG[u];
// 		MEInP-=nonadj_inP;
// 		MEInG-=nonadj_inG;
// 		//!warning: lcy's code
// 		P.reset(u);
// 		NP.set(u);
// 		for (int i = pstart[u]; i < pstart[u+1]; i++){
// 			// checklcy();
// 			neiInP[edges[i]]--, neiInG[edges[i]]--;
// 			if (C.test(edges[i])){  //if edges[i] in C
// 				C_heap_P->decrement(edges[i],1);
// 			}
// 			// checklcy();
// 		}
// 		// checklcy();
// 		//!warning: lcy's code
// 	}

// 	void C2X(int id){
// 		// checklcy();
// 		assert(inC(id));
// 		int u=PC[id];
// 		swapID(id,--C_end);
// 		int nonadj_inC=C_end-neiInG[u];
// 		MEInG-=nonadj_inC;
// 		//!warning: lcy's code
// 		C.reset(u);
// 		NC.set(u);
// 		C_heap_C->del(u);
// 		C_heap_P->del(u);
// 		for (int i = pstart[u]; i < pstart[u+1]; i++){
// 			neiInG[edges[i]]--,neiInC[edges[i]]=neiInG[edges[i]]-neiInP[edges[i]];
// 			if (C.test(edges[i])){  //if edges[i] in C
// 				C_heap_C->decrement(edges[i],1);
// 			}
// 		}
// 		// checklcy();
// 		//!warning: lcy's code
// 	}

// 	void X2C(int pop_cnt){
// 		while (pop_cnt--)
// 		{
// 			// checklcy();
// 			C_end++;
// 			int nonadj_inG=C_end-1-neiInG[PC[C_end-1]];
// 			MEInG+=nonadj_inG;
// 			//MEInP not change
// 			//!warning: lcy's code
// 			int u=PC[C_end-1];
// 			C.set(u);
// 			NC.reset(u);
// 			C_heap_P->insert(u,neiInP[u]);
// 			C_heap_C->insert(u,neiInC[u]);
// 			for (int i = pstart[PC[C_end-1]]; i < pstart[PC[C_end-1]+1]; i++){
// 				neiInG[edges[i]]++,neiInC[edges[i]]=neiInG[edges[i]]-neiInP[edges[i]];
// 				if(C.test(edges[i])){
// 					C_heap_C->increment(edges[i],1);
// 				}
// 			}
// 			// checklcy();
// 			//!warning: lcy's code
// 		}
// 	}

// 	int hiReduce(){
// 		int push_cnt=0;
// 		for (int id = P_end; id < C_end;)
// 			if (neiInG[PC[id]]>=C_end-2 && neiInP[PC[id]]==P_end) C2P(id), push_cnt++;
// 			else id++;
// 		return push_cnt;
// 	}

// 	int lwReduce(){
// 		int pop_cnt=0;
// 		for (int id = P_end; id < C_end; )
// 			if ( P_end - neiInP[PC[id]] + MEInP > k || 1+neiInG[PC[id]]+(k-MEInP)<=LB ) C2X(id), pop_cnt++;
// 			else id++;
// 		return pop_cnt;
// 	}

// 	int bound(int neis){
// 		return min(P_end+neis+(k-MEInP),UB);
// 	}
// 	int bound2(){
// 		return min(P_end+C_heap_P->key_num[P_end]+(k-MEInP),UB);
// 	}
// 	int candibound(){
// 		int totME=0;
// 		int C_B=0;
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--)
// 		{
// 			if ((totME+C_heap_P->key_num[deg]*(P_end-deg))<=k-MEInP){
// 				C_B+=C_heap_P->key_num[deg];
// 				totME+=C_heap_P->key_num[deg]*(P_end-deg);
// 			}
// 			else{
// 				C_B+=(k-MEInP-C_B)/(P_end-deg);
// 				// totME+=C_heap_P->key_num[deg]*((k-MEInP-C_B)/(P_end-deg));
// 				break;
// 			}
// 		}
// 		return C_B+P_end;
// 	}
// 	int candibound2(){
// 		int totME=0;
// 		int C_B2=0;
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for (int id = C_heap_P->head_s[deg]; id != n; id=C_heap_P->next_s[id]){
// 				if (totME+(P_end-C_heap_P->key_s[id])<=k-MEInP){
// 					totME+=(P_end-C_heap_P->key_s[id]);
// 					C_B2++;
// 				}
// 				else
// 					break;
// 			}
// 		}
// 		return P_end+C_B2;
// 	}
// 	//! warning: lcy's code
// 	int color_candibound(){
// 		int cnt_color = 0,k_cur=k-MEInP;
// 		int CC_bound=0;
// 		vector<int> lv,lv2;
// 		C_heap_C->load_heap_inc_order(lv);// put the vertices in decreasing non-degree order of C into lv
// 		C_heap_P->load_heap_dec_order(lv2);// put the vertices in increasing non-degree order of P into lv2
// 		Vset C_copy=C;
// 		while(!lv.empty()) {
// 			Vset cur = C_copy;
// 			++cnt_color;
// 			for (vector<int>::iterator it = lv.begin(); it != lv.end(); ) {
// 				if (cur[*it]) {
// 					cur &= rbitMtx[*it];
// 					C_copy.reset(*it);
// 					color_lable[*it]=cnt_color;
// 					it = lv.erase(it);
// 				} else ++it;
// 			}
//    		}
// 		fill(color_ME+1,color_ME+cnt_color+1,k_cur);
// 		fill(color_num+1,color_num+cnt_color+1,0);
// 		for (int i = 0; i < lv2.size(); i++)
// 		{
// 			if (color_ME[color_lable[lv2[i]]]-color_num[color_lable[lv2[i]]]-(P_end-neiInP[lv2[i]])>=0)
// 			{
// 				if (k_cur-(P_end-neiInP[lv2[i]])>=0)
// 				{
// 					CC_bound++;
// 					color_ME[color_lable[lv2[i]]]-=(color_num[color_lable[lv2[i]]]+(P_end-neiInP[lv2[i]]));
// 					k_cur-=(P_end-neiInP[lv2[i]]);
// 					color_num[color_lable[lv2[i]]]++;
// 				}
// 				else
// 					break;
// 			}
			
// 		}
// 		return P_end+CC_bound;
// 	}
// 	//! warning: lcy's code

// 	//! warning: lcy's code
// 	int color_candibound2(){
// 		int cnt_color=0,CC_bound=0;
// 		int k_cur=k-MEInP;
// 		vector<int> lv,lv2;
// 		C_heap_C->load_heap_inc_order(lv);// put the vertices in decreasing non-degree order of C into lv
// 		C_heap_P->load_heap_dec_order(lv2);// put the vertices in increasing non-degree order of P into lv2
// 		// fill(color_ME+1,color_ME+1+n,k_cur);
// 		// fill(color_num+1,color_num+n+1,0);
// 		for (int i = 0; i < lv.size(); i++){
// 			int cur_color=1;
// 			Vset cur_con;
// 			while (true){
// 				cur_con=bitcolor[cur_color]&bitMtx[lv[i]];
// 				if (cur_con!=0)
// 					cur_color++;
// 				else{
// 					bitcolor[cur_color].set(lv[i]);
// 					color_lable[lv[i]]=cur_color;
// 					cnt_color=max(cnt_color,cur_color);
// 					break;
// 				}
// 			}
// 		}
// 		fill(color_ME+1,color_ME+1+cnt_color,k_cur);
// 		fill(color_num+1,color_num+1+cnt_color,0);
// 		for (int i = 0; i < lv2.size(); i++){
// 			if (color_ME[color_lable[lv2[i]]]-color_num[color_lable[lv2[i]]]-(P_end-neiInP[lv2[i]])>=0)
// 			{
// 				if (k_cur-(P_end-neiInP[lv2[i]])>=0)
// 				{
// 					CC_bound++;
// 					color_ME[color_lable[lv2[i]]]-=(color_num[color_lable[lv2[i]]]+(P_end-neiInP[lv2[i]]));
// 					k_cur-=(P_end-neiInP[lv2[i]]);
// 					color_num[color_lable[lv2[i]]]++;
// 				}
// 				else
// 					break;
// 			}
// 		}
// 		for (int i = 1; i <= cnt_color; i++){
// 			bitcolor[i].reset();
// 		}
// 		return P_end+CC_bound;
// 	}
// 	//! warning: lcy's code

// 	//! warning: lcy's code
// 	int color_candibound3(){
// 		int cnt_color=0,CC_bound=0;
// 		int k_cur=k-MEInP;
// 		fill(color_ME+1,color_ME+n,k_cur);
// 		fill(color_num+1,color_num+n,0);
// 		for (int deg = C_heap_C->min_key; deg <= C_heap_C->max_key; deg++){
// 			for (int i = C_heap_C->head_s[deg]; i != C_heap_C->n; i=C_heap_C->next_s[i]){
// 				int cur_color=1;
// 				Vset cur_con;
// 				while (true){
// 					cur_con=bitcolor[cur_color]&bitMtx[i];
// 					if (cur_con!=0){
// 						cur_color++;
// 					}else{
// 						bitcolor[cur_color].set(i);
// 						color_lable[i]=cur_color;
// 						cnt_color=max(cnt_color,cur_color);
// 						break;
// 					}
// 				}
// 			}
// 		}
// 		// fill(color_ME+1,color_ME+1+cnt_color,k_cur);
// 		// fill(color_num+1,color_num+cnt_color+1,0);

// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for(int id=C_heap_P->head_s[deg];id!=C_heap_P->n;id=C_heap_P->next_s[id]){
// 				if (color_ME[color_lable[id]]-color_num[color_lable[id]]-(P_end-neiInP[id])>=0)
// 				{
// 					if (k_cur-(P_end-neiInP[id])>=0)
// 					{
// 						CC_bound++;
// 						color_ME[color_lable[id]]-=(color_num[color_lable[id]]+(P_end-neiInP[id]));
// 						k_cur-=(P_end-neiInP[id]);
// 						color_num[color_lable[id]]++;
// 					}
// 					else
// 						break;
// 				}
// 			}
// 		}
// 		for (int i = 1; i <= cnt_color; i++){
// 			bitcolor[i].reset();
// 		}
// 		return P_end+CC_bound;
// 	}
// 	//! warning: lcy's code
// 	int color_candibound4(){
// 		int max_color=0,CC_bound=0;
// 		int k_cur=k-MEInP;
// 		fill(color_ME,color_ME+n,k_cur);
// 		fill(color_num,color_num+n,0);
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for (int id = C_heap_P->head_s[deg]; id != C_heap_P->n; id=C_heap_P->next_s[id]){
// 				if (k_cur-(P_end-neiInP[id])<0)
// 					goto RESET;
// 				else{
// 					int cur_color=0;
// 					Vset con;
// 					while (true){//graph coloring
// 						con=bitcolor[cur_color]&bitMtx[id];
// 						if (con!=0)
// 							cur_color++;
// 						else{
// 							bitcolor[cur_color].set(id);
// 							color_lable[id]=cur_color;
// 							max_color=max(max_color,cur_color);
// 							break;
// 						}
// 					}
// 					if (color_ME[color_lable[id]]-color_num[color_lable[id]]-(P_end-neiInP[id])>=0)
// 					{
// 						CC_bound++;
// 						color_ME[color_lable[id]]-=(color_num[color_lable[id]]+(P_end-neiInP[id]));
// 						k_cur-=(P_end-neiInP[id]);
// 						color_num[color_lable[id]]++;
// 					}
// 				}
// 			}
// 		}
// RESET: 	for (int i = 0; i <= max_color; i++){
// 			bitcolor[i].reset();
// 		}
// 		return P_end+CC_bound;
// 	}
// 	int color_candibound5(){
// 		int max_color=0,CC_bound=0;
// 		int k_cur=k-MEInP;
// 		fill(color_ME,color_ME+n,k_cur);
// 		fill(color_num,color_num+n,0);
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for (int id = C_heap_P->head_s[deg]; id != C_heap_P->n; id=C_heap_P->next_s[id]){
// 				if (k_cur-(P_end-neiInP[id])<0)
// 					goto RESET;
// 				else{
// 					int cur_color=0;
// 					int cur_color2=0;
// 					int con;
// 					Vset con2;
// 					while (true){//graph coloring
// 						con=intsecColor(id,cur_color);//! error
// 						con2=bitcolor[cur_color]&bitMtx[id];
// 						if ((con==0 && con2!=0))
// 						{
// 							printf("error\n");
// 						}
						
// 						if (con!=0){//if changing into con2!=0 then it's correct
// 							// assert((con!=0 && con2 != 0)||(con==0 && con2==0));
// 							// assert(con!=0);
// 							cur_color++;
// 						}
// 						else{
// 							// assert(con==0);
// 							bitcolor[cur_color].set(id);
// 							colors[cur_color].bitcolor.set(id);
// 							color_lable[id]=cur_color;
// 							max_color=max(max_color,cur_color);
// 							break;
// 						}
// 					}
// 					// while (true){
// 					// 	con2=bitcolor[cur_color2]&bitMtx[id];
// 					// 	if (con2!=0){
// 					// 		cur_color2++;
// 					// 	}
// 					// 	else{
// 					// 		assert(cur_color2-cur_color==0);
// 					// 		bitcolor[cur_color2].set(id);
// 					// 		// if(cur_color2!=cur_color)
// 					// 		// 	printf("Error!\n");
// 					// 		// else	
// 					// 		// 	printf("Correct!\n");
// 					// 		color_lable[id]=cur_color2;
// 					// 		max_color=max(max_color,cur_color2);
// 					// 		break;
// 					// 	}
// 					// }
// 					if (color_ME[color_lable[id]]-color_num[color_lable[id]]-(P_end-neiInP[id])>=0)
// 					{
// 						CC_bound++;
// 						color_ME[color_lable[id]]-=(color_num[color_lable[id]]+(P_end-neiInP[id]));
// 						k_cur-=(P_end-neiInP[id]);
// 						color_num[color_lable[id]]++;
// 					}
// 				}
// 			}
// 		}
// RESET: 	for (int i = 0; i <= max_color; i++){
// 			bitcolor[i].reset();
// 			colors[i].bitcolor.clear();
// 		}
// 		return P_end+CC_bound;
// 	}
// 	int Facolor_candibound(){
// 		int max_color=1,CC_bound=0;
// 		int k_cur=k-MEInP;
// 		fill(color_ME+1,color_ME+n+1,k_cur);
// 		fill(color_num+1,color_num+n+1,0);
// 		colors[1].bitcolor.clear();
// 		colors[2].bitcolor.clear();
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for (int id = C_heap_P->head_s[deg]; id != C_heap_P->n; id=C_heap_P->next_s[id]){
// 				if (k_cur-(P_end-neiInP[id])<0)
// 					goto RET;
// 				else{
// 					int cur_color=1;
// 					while (1)
// 					{
// 						if(intsecColor(id,cur_color)!=0)
// 							cur_color++;
// 						else
// 							break;
// 					}
// 					colors[cur_color].bitcolor.set(id);
// 					color_lable[id]=cur_color;
// 					if (cur_color>max_color){
// 						max_color=cur_color;
// 						colors[max_color+1].bitcolor.clear();
// 					}
// 					if (color_ME[color_lable[id]]-color_num[color_lable[id]]-(P_end-neiInP[id])>=0){
// 						CC_bound++;
// 						color_ME[color_lable[id]]-=(color_num[color_lable[id]]+(P_end-neiInP[id]));
// 						k_cur-=(P_end-neiInP[id]);
// 						color_num[color_lable[id]]++;
// 					}
// 				}
// 			}
// 		}
// RET: 	return P_end+CC_bound;
// 	}
// 	void store(int newLB){
// 		KDC.resize(LB=newLB);
// 		assert(KDC.size()<=n);
// 		for (int i = 0; i < LB; i++) KDC[i]=PC[i];
// 	}

// 	void showcolors(){
// 		printf("The num of C is: %d; The num of P is: %d\n",C_heap_P->size,P_end);
// 		printf("Missing edges in P is: %d\n",MEInP);
// 		for (int deg = C_heap_P->max_key; deg >= C_heap_P->min_key; deg--){
// 			for (int id = C_heap_P->head_s[deg]; id != C_heap_P->n; id=C_heap_P->next_s[id]){
// 				printf("%d:%d:%d, ",id,neiInP[id],color_lable[id]);
// 			}
// 		}
// 		printf("\n");
// 	}
// 	void branch(){
// 		// if(C_end <= LB) return; if(P_end > LB) store(P_end);
// 		tree_cnt++;
// 		int push_cnt=hiReduce();
// 		int pop_cnt=lwReduce();
// 		int pid=P_end, neis=0;
// 		int min_d_inP=0;
// 		if(C_end <= LB) goto REC;
// 		if (P_end > LB) store(P_end);
// 		if (MEInG<=k) {store(C_end); goto REC;}

		
// 		// for(int id=P_end;id<C_end;id++) (neiInP[PC[pid]]>neiInP[PC[id]])&&(pid=id),neis+=(neiInP[PC[id]]==P_end);
// 		C_heap_P->get_min(pid,min_d_inP);
// 		pid=PC_rid[pid];
// 		// if (bound2()<=LB) goto REC;
// 		// if (candibound()<=LB) goto REC;
// 		// if (candibound2()<=LB) goto REC;
// 		// if(color_candibound()<=LB) goto REC;
// 		// if(color_candibound2()<=LB) goto REC;
// 		// if(color_candibound3()<=LB) goto REC;
// 		// if(color_candibound4()<=LB) goto REC;
// 		if(color_candibound5()<=LB) goto REC;
// 		// if(Facolor_candibound()<=LB) goto REC;
// 		// if (color_candibound4()<=LB)
// 		// {
// 		// 	if (color_candibound4()==30&&LB==30)
// 		// 	{
// 		// 		showcolors();
// 		// 	}
// 		// 	goto REC;
// 		// }
		
// 		// if (P_end>=LB/4)
// 		// {
// 		// 	if(color_candibound()<=LB) goto REC;
// 		// }
// 		C2P(pid); branch(); P2X(); branch(); X2C(1);
// 	REC:
// 		X2C(pop_cnt); P2C(push_cnt);
// 	}

// 	void search(){
// 		// checklcy();
// 		C2P(0); branch(); P2C(0);
// 		// checklcy();
// 	}

// 	void swapID(int i, int j) {
// 		swap(PC[i], PC[j]);
// 		PC_rid[PC[i]] = i;
// 		PC_rid[PC[j]] = j;
// 	}

// };