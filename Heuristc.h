#ifndef _HEURISTC_H_
#define _HEURISTC_H_
#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
using namespace std;
class HeuriSearcher{
public:
	ui n;
	ept m;
	ept* pstart;
	ui* edges;
    // ui k; //total missing edges
	double gamma;
    ui K;
	ui maxDeg;
    ui LB; // the lowerbound
    ui UB; // the upperbound
	
    ui n_s;
    ept m_s;
	std::vector<ui> id_s;
	std::vector<ui> rid_s;
	std::vector<ui> deg_s;
    std::vector<ui> core_s;
    std::vector<ui> seq_s;
	bool* exist_s;
    bool* delete_s;
    char* vis_s;
	std::vector<ui> heuMDS;
	std::vector<std::vector<ui>> adjList;
    ListLinearHeap *s_heap;

	// Your code
	HeuriSearcher(ui _n,ui _m, ept *_pstart, ui *_edges, double _gamma, ui _maxDeg, std::vector<ui>& _QC){
		n=_n, m=_m, gamma=_gamma;
		this->pstart=_pstart;
		this->edges=_edges;
		maxDeg=_maxDeg;
		exist_s=new bool[n];
        delete_s=new bool[n];
        vis_s=new char[1+maxDeg];
		memset(exist_s, false, n*sizeof(bool));
        memset(delete_s, false, n*sizeof(bool));
        LB=ui(_QC.size());
        
	}
    HeuriSearcher(ui _n,ui _m, ept *_pstart, ui *_edges, ui _K, ui _maxDeg, std::vector<ui>& _KDC){
		n=_n, m=_m, gamma=_K;
		this->pstart=_pstart;
		this->edges=_edges;
		maxDeg=_maxDeg;
		exist_s=new bool[n];
        delete_s=new bool[n];
        vis_s=new char[1+maxDeg];
		memset(exist_s, false, n*sizeof(bool));
        memset(delete_s, false, n*sizeof(bool));
        LB=ui(_KDC.size());
	}
 	void induceNei(ui u){
		n_s=0; m_s=0;
		//1. get the max id
		ui maxId=u;
		for (ept i = pstart[u]; i < pstart[u+1]; i++) maxId=max(maxId, edges[i]);
		rid_s.resize(1+maxId);
		//2. get the veritces
		id_s.push_back(u); rid_s[u]=id_s.size()-1; exist_s[u]=true;
		
		for (ept i = pstart[u]; i < pstart[u+1]; i++){
			id_s.push_back(edges[i]);
            exist_s[edges[i]]=true;
			rid_s[edges[i]]=id_s.size()-1;
		}
        n_s=id_s.size(); 
        core_s.resize(n_s);//core_s initialization
        assert(seq_s.empty());
        for (ui i = 0; i < n_s; i++) seq_s.push_back(i);  
		//3. get the edges
		for(auto v: id_s) exist_s[v]=true;
		adjList.resize(ui(id_s.size()));
		for (ui i = 0; i < id_s.size(); i++){
			ui v=id_s[i];
			for (ui j = pstart[v]; j < pstart[v+1]; j++){
				ui w=edges[j];
				if(exist_s[w]) {
					adjList[rid_s[v]].push_back(rid_s[w]);
				}
			}
		}
		//4. construct the degrees
		for (ui i = 0; i < n_s; i++) deg_s.push_back(ui(adjList[i].size())), m_s+=deg_s[i];
        m_s/=2;
        //5. clear the exist
        for (ui i = 0; i < n_s; i++) exist_s[id_s[i]]=false;   
	}
    ui induceNeiBack(ui u){
        n_s=m_s=0;
		//1. get the max id
		ui maxId=u;
		for (ui i = pstart[u]; i < pstart[u+1]; i++) {
            if(delete_s[edges[i]]) continue;
            maxId=max(maxId, edges[i]);
        }
		rid_s.resize(1+maxId);
		//2. get the veritces
		id_s.push_back(u); rid_s[u]=id_s.size()-1; exist_s[u]=true;
		
		for (ept i = pstart[u]; i < pstart[u+1]; i++){
            if(delete_s[edges[i]]) continue;
			id_s.push_back(edges[i]);
            exist_s[edges[i]]=true;
			rid_s[edges[i]]=id_s.size()-1;
		}
        n_s=(ui)id_s.size(); 
        core_s.resize(n_s);//core_s initialization
        assert(seq_s.empty());
        for (ui i = 0; i < n_s; i++) seq_s.push_back(i);  
		//3. get the edges
		for(auto v: id_s) exist_s[v]=true;
		adjList.resize(ui(id_s.size()));
		for (ui i = 0; i < id_s.size(); i++){
			ui v=id_s[i];
			for (ept j = pstart[v]; j < pstart[v+1]; j++){
				ui w=edges[j];
				if(exist_s[w]) {
					adjList[rid_s[v]].push_back(rid_s[w]);
				}
			}
		}
		//4. construct the degrees
		for (ui i = 0; i < n_s; i++) deg_s.push_back(ui(adjList[i].size())), m_s+=deg_s[i];
        m_s/=2;
        //5. clear the exist
        for (ui i = 0; i < n_s; i++) exist_s[id_s[i]]=false;  
        return n_s;
    }
    ui degenHeu(ui u_0, ui u_deg, ListLinearHeap *heap, ui _kind){
        ui subSz=1+u_deg, s_maxcore=0;
        ui edgeCnt_s=m_s, idx=subSz;
        // printf("u0: %d, subSz: %d\n",u_0, subSz);
        memset(vis_s, 0, subSz*sizeof(char));
        heap->init(subSz, subSz-1, seq_s.data(), deg_s.data());
        //delete the vertex by degeneracy order
        for (ui i = 0; i < subSz; i++){
            switch (_kind){
                case 1:
                    if(idx == subSz && edgeCnt_s >= ceil(gamma*(1.0*(subSz - i)*(subSz - i - 1)/2.0)) ) idx = i;
                    break;
                case 2:
                    if(idx == subSz && edgeCnt_s >= (subSz - i)*(subSz - i - 1)/2 - K) idx = i;
                    break;
            }
            // if(idx == subSz && edgeCnt_s >= ceil(gamma*(1.0*(subSz - i)*(subSz - i - 1)/2.0)) ) idx = i;
            ui u, key; heap->pop_min(u,key);
            if(key> s_maxcore) s_maxcore=key; 
            core_s[u]=s_maxcore; vis_s[u]=1; seq_s[i]=u;
            for(auto v: adjList[u]){
                if(vis_s[v]==0) heap->decrement(v,1);
            }
            edgeCnt_s-=key;
        }
        if(subSz-idx>heuMDS.size()){
            printf("renew heuMDS:%d\n", subSz-idx);
            LB=subSz-idx;
            heuMDS.clear();
            for (ui i = idx; i < subSz; i++) heuMDS.push_back(id_s[seq_s[i]]);
            if(!checkExist(u_0, heuMDS)) {printf("error, not include u: %d\n", u_0);}//this should be useless because u_0 is adjacent to all vertices 
        }

        seq_s.clear();
        return n;
    }
	void clear(){
		id_s.clear(), rid_s.clear(); seq_s.clear();
		for(auto edgeVec: adjList) edgeVec.clear();
		adjList.clear(); deg_s.clear();
	}
    bool checkExist(ui u, std::vector<ui> &_QC){
        bool flag=false;
        for(auto v:_QC){
            if(v==u) {flag=true;break;}
        }
        return flag;
    }
	void search(std::vector<ui> &_MDS, ui _kind){
        Timer t;
        ui max_n=1+maxDeg;
        ListLinearHeap* s_heap=new ListLinearHeap(max_n, max_n-1);
        for (ui u = 0; u < n; u++){
            //1. build the subgraph
            // printf("heuristic search subgraph: %d\n",u);
            induceNei(u);
            ui u_deg=adjList[rid_s[u]].size();
            this->UB=degenHeu(u, u_deg, s_heap, _kind);
            this->clear();
        }
        if(LB>ui(_MDS.size())){
            _MDS.resize(LB);
            for (ui i = 0; i < LB; i++){
                _MDS[i]=heuMDS[i];
            }
            // printf("renew the result in heuristic search\n");
        }
        printf("#HeuSize=%u\n#HeuTime=%.2f\n", LB, double(t.elapsed())/1000000);
	}
    void degenSearch(std::vector<ui> &_MDS, ui* seq, ui _kind){
        Timer t;
        ui max_n=1+maxDeg;
        ui subMax=0, subSz=0;
        ListLinearHeap* s_heap=new ListLinearHeap(max_n, max_n-1);
        ui pre_size = (ui)heuMDS.size(), updt_gap = 0;
        for (ui i = 0; i < n; i++){
            //1. build the subgraph
            // printf("heuristic search subgraph: %d\n",u);
            ui u = seq[i];
            subSz=induceNeiBack(u);
            subMax=max(subMax, subSz);
            ui u_deg=adjList[rid_s[u]].size();
            this->UB=degenHeu(u, u_deg, s_heap, _kind);
            this->clear();
            delete_s[u]=true;
            if((ui)heuMDS.size()<=pre_size) updt_gap++;
            else pre_size=(ui)heuMDS.size(), updt_gap=0;
            if(updt_gap>3000) break;
        }
        if(LB>ui(_MDS.size())){
            _MDS.resize(LB);
            for (ui i = 0; i < LB; i++){
                _MDS[i]=heuMDS[i];
            }
            // printf("renew the result in heuristic search\n");
        }
        printf("#HeuSize=%u\n#HeuSubMax=%u\n#HeuTime=%.2f\n", LB, subMax,double(t.elapsed())/1000000);
    }
};



#endif