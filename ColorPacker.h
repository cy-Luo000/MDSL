#ifndef _COLOR_PACKER_H_
#define _COLOR_PACKER_H_

#include "Utility.h"
using namespace std;
class colorPacker
{
private:
	/* data */
public:
	ui n;// the length of the color vector
	ui *colLab;//the color of the element in color vector(color lable)
	vector<vector<ui>> AdjList;//the graph of non mutual pairs
	vector<ui> vChose;
	vector<ui> newWeiB; //the new weights of the chosen vertices in the same colorset
	ui maxWeight;
	bool *colMtx;
	ui colorNum;
	ui maxWsum;

	colorPacker(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz);
	colorPacker(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz,  ui KCur);
	colorPacker(ui n);
	colorPacker();
	void set(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz);
	void set(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz,  ui KCur);
	void set(std::vector<std::pair<ui, ui>>& conflictPairs, ui st, ui end, vector<ui>& colVec, ui ColSz,  ui KCur);
	void reset();
	void resetL();
	void clear();
	void coloring();
	void coloring(ui *colVec, ui* neiInP, ui P_end);
	void coloring(vector<ui>& colVec, ui* neiInP, ui P_end);
	ui colPack(ui *colVec, ui colSz ,ui *neilInP,ui P_end,ui threshold);
	void getNewWeights(ui *colVec, ui ColSz, ui *neiInP, ui P_end, ui kCur);
	void getNewWeights(vector<ui> & colVec, ui ColSz, ui *neiInP, ui P_end, ui kCur);
	~colorPacker();
};
colorPacker::colorPacker(ui n_){// n_ is the upperbound of n
	colLab=new ui[n_];
	memset(colLab, 0, n_*sizeof(ui));
	AdjList.resize(n_);
	for(int i=0; i<AdjList.size(); i++){
		AdjList[i].clear();
	}
	colMtx=new bool[n_*n_];
	memset(colMtx, false, n_*n_*sizeof(bool));
	n=0, colorNum=0, maxWsum=0;maxWeight = 0;
}
colorPacker::colorPacker(){// n_ is the upperbound of n
	n=0, colorNum=0, maxWsum=0;maxWeight = 0;
	AdjList.clear();newWeiB.clear();vChose.clear();
}
void colorPacker::set(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz,  ui KCur){
	n=ColSz;// n is the color size
	maxWsum=0;
	newWeiB.resize(KCur+1); newWeiB.assign(KCur+1, 0);
	for (ui i = 0; i < ColSz; i++){
		ui u=colVec[i];
		for (ui j = i+1; j < ColSz; j++){
			ui v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
}
void colorPacker::set(std::vector<std::pair<ui, ui>>& conflictPairs, ui st, ui end, vector<ui>& colVec, ui ColSz,  ui KCur){
	// printf("get into the colpacker set\n");
	// // printf("lcysddd\n");
	// printf("colPacker: colsz=%u\n",ColSz);
	// printf("the size of colpack: %u\n", n);
	this->n=ColSz;// n is the color size
	// printf("lcysddd\n");
	AdjList.resize(ColSz); //the adjlist has at most n vertices
	for (ui i = 0; i < ColSz; i++) AdjList[i].clear();
	
	// printf("AdjList init complete\n");
	// exit(0);
	maxWsum=0;
	
	newWeiB.resize(KCur+1);  newWeiB.assign(KCur+1, 0);// the weight bucket to sort the new C
	std::unordered_map<ui, ui> rid;
	// std::unordered_set<pair<ui,ui>> confEdges;
	std::unordered_set<string> confEdges2;
	// printf("colpacker set init complete\n");
	for(ui i=0; i<ColSz; i++){ //get the rid of each vertex u
		ui u = colVec[i];
		rid[u]=i;
	}
	for (ui i = st; i < end; i++){
		std::pair<ui, ui> p = conflictPairs[i];
		ui u = p.first, v= p.second;
		ui j = rid[u], k=rid[v];
		// confEdges.insert(make_pair(j,k)), confEdges.insert(make_pair(k,j));
		confEdges2.insert(to_string(j)+","+to_string(k)), confEdges2.insert(to_string(k)+","+to_string(j));
		//
	}
	//construct the graph
	for(ui i=0; i < ColSz; i++){
		for(ui j=i+1; j<ColSz;j++){
			// pair<ui,ui> p = make_pair(i,j);
			// if(confEdges.find(p)!=confEdges.end()) continue;
			string edge = to_string(i)+","+to_string(j);
			if(confEdges2.find(edge)!=confEdges2.end()) continue;
			AdjList[i].push_back(j);
			AdjList[j].push_back(i);
		}
	}
}
void colorPacker::set(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz){
	n=ColSz;// n is the color size
	maxWsum=0;
	for (ui i = 0; i < ColSz; i++){
		ui u=colVec[i];
		for (ui j = i+1; j < ColSz; j++){
			ui v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
}
void colorPacker::reset(){
	//reset the color lable, adjlist and color mtx
	memset(colLab, 0, n*sizeof(ui));
	for (ui i = 0; i < n; i++) AdjList[i].clear();
	memset(colMtx, false, n*n*sizeof(bool));
	vChose.clear(); newWeiB.clear();
	colorNum=0, n=0, maxWsum=0; maxWeight=0;
}
void colorPacker::resetL(){
	AdjList.clear();
	vChose.clear(); newWeiB.clear();
	colorNum=0, n=0, maxWsum=0;
}
colorPacker::colorPacker(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz)
{
	for (ui i = 0; i < ColSz; i++){
		ui u=colVec[i];
		if(u<0){
			printf("u: %d, i: %d, sz: %d\n",u, i, ColSz);
			exit(0);
		}
#ifdef _DBUG_
		assert(u>=0);
		assert(u<n_);
#endif
	}
	
	n=ColSz;// n is the color size
	colorNum=0;
	AdjList.resize(n);
	//construct the graph
	for (ui i = 0; i < ColSz; i++){
		ui u=colVec[i];
		for (ui j = i+1; j < ColSz; j++){
			ui v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
	colMtx=new bool[n*n];
	memset(colMtx,false, n*n*sizeof(bool));
	colLab=new ui[n];
	memset(colLab,0, n*sizeof(ui));
}
colorPacker::colorPacker(std::vector<std::vector<bool>>& MuEx, ui *colVec, ui ColSz,  ui kCur){
	n=ColSz;// n is the color size
	colorNum=0;
	maxWeight = 0;
	AdjList.resize(n);
	newWeiB.resize(kCur+1);
	newWeiB.assign(kCur+1, 0);
	//construct the graph
	for (ui i = 0; i < ColSz; i++){
		ui u=colVec[i];
		for (ui j = i+1; j < ColSz; j++){
			ui v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
}
void colorPacker::clear(){
	n=0;
	if(colLab!=NULL){
		delete[] colLab;
		colLab=NULL;
	}
	if(AdjList.size()>0){
		AdjList.clear();
	}
	if(!vChose.empty()) vChose.clear();
	if(!newWeiB.empty()) newWeiB.clear();
	if(colMtx!=NULL){
		delete[] colMtx;
		colMtx=NULL;
	}
}
void colorPacker::coloring(){
	int maxCol=-1;
	// do the coloring
	//! does the order influence the result???
	for (ui i = 0; i < n; i++){
		int col=0;
		while (colMtx[n*i+col]) col++;
		if(col>maxCol) {
			maxCol=col;
			vChose.push_back(i);
		}
		colLab[i]=col;
		for (auto j:AdjList[i]) colMtx[n*j+col]=true;
		// vChose.push_back(i);
	}
	colorNum=maxCol+1;
	return;
}
void colorPacker::coloring(ui* colVec, ui *neiInP, ui P_end){
	int maxCol=-1;
	// do the coloring
	for (ui i = 0; i < n; i++){
		int col=0;
		while (colMtx[n*i+col]) col++;
		if(col>maxCol) {
			maxCol=col;
			vChose.push_back(i);//get the index in the colorvec
			maxWsum+=(P_end-neiInP[colVec[i]]+vChose.size()-1);
		}
		colLab[i]=col;
		for (auto j:AdjList[i]) colMtx[n*j+col]=true;
	}
	colorNum=maxCol+1;
	return;
}
void colorPacker::coloring(vector<ui>& colVec, ui *neiInP, ui P_end){
	int maxCol=-1;
	vector<ui> colorlabel; //to record the color of each vertex
	vector<bool> colvis; // to record if the color has used
	
	colorlabel.resize(this->n); colvis.resize(this->n);
	colorlabel.assign(n, -1); colvis.assign(n, false);

	//begin coloring
	for(ui i=0; i<n; i++){
		ui col=0;
		for(auto j:AdjList[i]){
			if(colorlabel[j]>=0) colvis[j]=true;
		}
		while(colvis[col]) col++;
		if(col>maxCol){
			maxCol=col;
			vChose.push_back(i);
		}
		colorlabel[i] = col;
		colvis.assign(n, false);
	}


	//recover
	colorlabel.clear(); colvis.clear();
	colorNum=maxCol+1;
	return;

}
ui colorPacker::colPack(ui *colVec, ui colSz ,ui *neiInP,ui P_end,ui threshold){
	ui choseNum=0,mssEdge=0;
	if(threshold>=maxWsum) return vChose.size();
	for (ui i = 0; i < vChose.size(); i++){
		ui u=colVec[vChose[i]];
		if(mssEdge+i+(P_end-neiInP[u])<=threshold) choseNum++;
		else break;
		mssEdge+=(i+P_end-neiInP[u]);
	}
	return choseNum;
}
void colorPacker::getNewWeights(ui *colVec, ui ColSz, ui *neiInP, ui P_end, ui kCur){
	for(ui i = 0; i < vChose.size(); i++){
		ui u = colVec[vChose[i]];
		ui wei = i+P_end-neiInP[u];
		if(wei <= kCur){
			// printf("weight = %u, cur k = %u, bucksize=%d\n", wei, kCur, newWeiB.size());
			maxWeight = max(maxWeight, wei);
			// printf("here2.1\n");
			newWeiB[wei]++;
		}
	}
	// printf("get out the getNewWeights\n");
}
void colorPacker::getNewWeights(vector<ui> & colVec, ui ColSz, ui *neiInP, ui P_end, ui kCur){
	for (ui i = 0; i < vChose.size(); i++){
		ui u = colVec[vChose[i]];
		ui wei = i+P_end-neiInP[u];
		if(wei <= kCur){
			maxWeight = max(maxWeight, wei);
			newWeiB[wei]++;
		}
	}
}
colorPacker::~colorPacker()
{
	n=0;
	if(colLab!=NULL){
		delete[] colLab;
		colLab=NULL;
	}
	if(!AdjList.empty()){
		AdjList.clear();
	}
	if(!vChose.empty()) vChose.clear();
	if(colMtx!=NULL){
		delete[] colMtx;
		colMtx=NULL;
	}
	if(!newWeiB.empty()) newWeiB.clear();
}

#endif