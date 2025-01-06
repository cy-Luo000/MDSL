#ifndef _MIP_H_
#define _MIP_H_
#include <gurobi_c++.h>
#include<vector>
#include "Utility.h"
using namespace std;
int feasible=0, unfeasible=0;

class MQCMIP;
class MKDCMIP;
class MQCMIP
{
private:
    int n,m;
    int lb, ub;
    double gamma;
    vector<pair<int,int>> edges;
    ept * pstart;
    ept* pend;
    ui *edges_ptr;
public:
    MQCMIP(/* args */);
    // MIP(int _n, double _gamma, int _lb);
    ~MQCMIP();
    void load_graph(ui _n, ept *_pstart, ept *_pend, ui* _edges,ui _lb, double _gamma);
    int dist(ui u, ui v);// to decide whether u and v are 2 hop
    int MIPSolve();
    float getUB2hop();
    int MIPSolve2();
    int MIPSolve2hop();
    void load_subgraph(ui _n, vector<pair<ui, ui>>& _vp, ui _lb, double _gamma);
    void printinfo();
    void getCN(vector<ui>& cmm_neighbors, ui u, ui v);
};

MQCMIP::MQCMIP(/* args */){
    n=m=0;
    lb=0, ub=INT32_MAX;
    gamma=1.0;
}
int MQCMIP::dist(ui u, ui v){
    unordered_set<ui> neighbors;
    for (ept i = pstart[u]; i < pend[u]; i++){
        ui u_nei = edges_ptr[i];
        if(u_nei == v) return 1;
        neighbors.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pend[v]; i++){
        ui v_nei = edges_ptr[i];
        if(neighbors.find(v_nei)!=neighbors.end()) return 2;
    }
    return 3;
}
void MQCMIP::load_graph(ui _n, ept *_pstart, ept *_pend, ui* _edges,ui _lb, double _gamma){
    n=_n; lb=_lb; gamma=_gamma;
    pstart=_pstart, pend=_pend, edges_ptr = _edges;
    for (ui u = 0; u < n; u++){
        for (ept j = _pstart[u]; j < _pend[u]; j++){
            int v=_edges[j];
            if(u<v) edges.push_back(make_pair(u,v)), m++;
        }
    }
    ub=floor(0.5+0.50*sqrt(1+8.0*m/gamma));
    ub=min(n, ub);
}
void MQCMIP::load_subgraph(ui _n, vector<pair<ui, ui>>& _vp, ui _lb, double _gamma){
    n=_n; lb=_lb; gamma=_gamma;
    // maxSubSz=max(maxSubSz, n);
    for(auto &p: _vp){
        if(p.first>p.second){
            int a=p.second;
            p.second=p.first, p.first=a;
        }
    }
    sort(_vp.begin(), _vp.end());
    for(auto e: _vp){
        if(!edges.empty()){
            pair<int,int> eback=edges.back();
            if(eback.first==e.first && eback.second==e.second) continue;
        }
        if(e.first>=n || e.second>=n){
            printf("edge error\n");
        }
        edges.push_back(e), m++;
    }
    ub=floor(0.5+0.50*sqrt(1+8.0*m/gamma));
    ub=min(n, ub);
}
void MQCMIP::printinfo(){
    printf("vertex num: %d, edge num: %d, density: %.2f\n", n, m, 2.0*m/(n*(n-1)));
}
void MQCMIP::getCN(vector<ui>& cmm_neighbors, ui u, ui v){
    cmm_neighbors.clear();
    unordered_set<ui> neighbors;
    for (ept i = pstart[u]; i < pend[u]; i++){
        ui u_nei = edges_ptr[i];
        neighbors.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pend[v]; i++){
        ui v_nei = edges_ptr[i];
        if (neighbors.find(v_nei)!=neighbors.end()){
            cmm_neighbors.push_back(v_nei);
        }
    }
}
int MQCMIP::MIPSolve(){
    GRBEnv env=GRBEnv(true);
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    // model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    // model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d\n", lb, ub);
    for (int k = lb; k <= ub; k++){
        eNumLB+= gamma*k*(k-1)*z[k]/2;
    }
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();

    int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Optimal solution found" << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
        cout << "Time limit reached, no optimal solution found within the specified time." << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else {
        cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    }
    cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    return obj_gurobi;
}
int MQCMIP::MIPSolve2(){
    GRBEnv env=GRBEnv(true);
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    // model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    // model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    //add 2hop constraint
    vector<ui> cmm_neighbors;
    for (int u = 0; u < n; u++){
        for (int v = u+1; v < n; v++){
            if(dist(u,v)==2){
                getCN(cmm_neighbors, u, v);
                GRBLinExpr constr=0;
                constr+=(x[u]+x[v]);
                for(auto cmm_nei:cmm_neighbors){
                    constr-=x[cmm_nei];
                }
                model.addConstr(constr<=1, "d_"+to_string(u)+to_string(v));
            }else if(dist(u,v)==3){
                model.addConstr(x[u]+x[v]<=1, "d_"+to_string(u)+to_string(v));
            }
        }
    }
    

    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d\n", lb, ub);
    for (int k = lb; k <= ub; k++){
        eNumLB+= gamma*k*(k-1)*z[k]/2;
    }
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++){
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++){
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");
    
    //optimize
    // Optimize model
    model.optimize();

    int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Optimal solution found" << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
        cout << "Time limit reached, no optimal solution found within the specified time." << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else {
        cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    }
    cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    return obj_gurobi;
}
int MQCMIP::MIPSolve2hop(){
    // lb=1;
    GRBEnv env=GRBEnv(true);
    
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
    // model.getEnv().set(GRB_IntParam_ResultFlag, 1); 
    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+"_"+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        // if(edges.first>)
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d, n: %d\n", lb, ub, n);
    for (int k = lb; k <= ub; k++){
        eNumLB+= gamma*k*(k-1)*z[k]/2;
    }
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if(status==GRB_OPTIMAL){
        int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
        feasible++;
        return obj_gurobi;
    }
    // // Print the status
    // switch (status) {
    //     case GRB_OPTIMAL:
    //         std::cout << "Optimal solution found." << std::endl;
    //         break;
    //     case GRB_INF_OR_UNBD:
    //         std::cout << "Model is infeasible or unbounded." << std::endl;
    //         break;
    //     case GRB_INFEASIBLE:
    //         std::cout << "Model is infeasible." << std::endl;
    //         break;
    //     case GRB_UNBOUNDED:
    //         std::cout << "Model is unbounded." << std::endl;
    //         break;
    //     case GRB_TIME_LIMIT:
    //         std::cout << "Time limit reached." << std::endl;
    //         break;
    //     case GRB_ITERATION_LIMIT:
    //         std::cout << "Iteration limit reached." << std::endl;
    //         break;
    //     case GRB_NODE_LIMIT:
    //         std::cout << "Node limit reached." << std::endl;
    //         break;
    //     // Add more cases as needed
    //     default:
    //         std::cout << "Optimization ended with status " << status << "." << std::endl;
    //         break;
    // }
    // model.computeIIS();
    // model.write("model.ilp");
    // int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    // if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    //     cout << "Optimal solution found" << endl;
    //     cout << "Obj: " << obj_gurobi << endl;
    // } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    //     cout << "Time limit reached, no optimal solution found within the specified time." << endl;
    //     cout << "Obj: " << obj_gurobi << endl;
    // } else {
    //     cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    // }
    // cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    // cout<<11111<<endl;
    // return obj_gurobi;
    return 0;
}
float MQCMIP::getUB2hop(){
    // lb=1;
    GRBEnv env=GRBEnv(true);
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
    // model.getEnv().set(GRB_IntParam_ResultFlag, 1); 
    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+"_"+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        // if(edges.first>)
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d, n: %d\n", lb, ub, n);
    for (int k = lb; k <= ub; k++){
        eNumLB+= gamma*k*(k-1)*z[k]/2;
    }
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if(status==GRB_OPTIMAL){
        int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
        feasible++;
        return obj_gurobi;
    }
    return 0.0;
}
MQCMIP::~MQCMIP()
{
    if(!edges.empty()){
        edges.clear();
    }
}

class MKDCMIP
{
private:
    int n,m;
    int lb, ub;
    int s;
    vector<pair<int,int>> edges;
    ept * pstart;
    ept* pend;
    ui *edges_ptr;
public:
    MKDCMIP(/* args */);
    // MIP(int _n, int _s, int _lb);
    ~MKDCMIP();
    void load_graph(ui _n, ept *_pstart, ept *_pend, ui* _edges,ui _lb, ui _s);
    int dist(ui u, ui v);// to decide whether u and v are 2 hop
    void getCN(vector<ui>& cmm_neighbors, ui u, ui v);
    int MIPSolve();
    int MIPSolve2();
    int MIPSolve2hop();
    float getUB2hop();
    void load_subgraph(ui _n, vector<pair<ui, ui>>& _vp, ui _lb, ui _s);
    void printinfo();
};

MKDCMIP::MKDCMIP(/* args */)
{
    n=m=0;
    lb=0, ub=INT32_MAX;
    s=1;
}
void MKDCMIP::load_graph(ui _n, ept *_pstart, ept *_pend, ui* _edges,ui _lb, ui _s){
    n=_n; lb=_lb; s=_s;
    pstart=_pstart, pend=_pend, edges_ptr = _edges;
    for (int u = 0; u < n; u++){
        for (int j = _pstart[u]; j < _pend[u]; j++){
            int v=_edges[j];
            if(u<v) edges.push_back(make_pair(u,v)), m++;
        }
    }
    ub=n;
}
int MKDCMIP::dist(ui u, ui v){
    unordered_set<ui> neighbors;
    for (ept i = pstart[u]; i < pend[u]; i++){
        ui u_nei = edges_ptr[i];
        if(u_nei == v) return 1;
        neighbors.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pend[v]; i++){
        ui v_nei = edges_ptr[i];
        if(neighbors.find(v_nei)!=neighbors.end()) return 2;
    }
    return 3;
}
void MKDCMIP::getCN(vector<ui>& cmm_neighbors, ui u, ui v){
    cmm_neighbors.clear();
    unordered_set<ui> neighbors;
    for (ept i = pstart[u]; i < pend[u]; i++){
        ui u_nei = edges_ptr[i];
        neighbors.insert(u_nei);
    }
    for (ept i = pstart[v]; i < pend[v]; i++){
        ui v_nei = edges_ptr[i];
        if (neighbors.find(v_nei)!=neighbors.end()){
            cmm_neighbors.push_back(v_nei);
        }
    }
}
void MKDCMIP::load_subgraph(ui _n, vector<pair<ui, ui>>& _vp, ui _lb, ui _s){
    n=_n; lb=_lb; s=_s;
    // maxSubSz=max(maxSubSz, n);
    for(auto &p: _vp){
        if(p.first>p.second){
            int a=p.second;
            p.second=p.first, p.first=a;
        }
    }
    sort(_vp.begin(), _vp.end());
    for(auto e: _vp){
        if(!edges.empty()){
            pair<int,int> eback=edges.back();
            if(eback.first==e.first && eback.second==e.second) continue;
        }
        if(e.first>=n || e.second>=n){
            printf("edge error\n");
        }
        edges.push_back(e), m++;
    }
    ub=n;
}
void MKDCMIP::printinfo(){
    printf("vertex num: %d, edge num: %d, density: %.2f\n", n, m, 2.0*m/(n*(n-1)));
}

int MKDCMIP::MIPSolve(){
    GRBEnv env=GRBEnv(true);
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    // model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    // model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d\n", lb, ub);
    for (int k = lb; k <= ub; k++){
        eNumLB+= k*(k-1)*z[k]/2;
    }
    eNumLB-=s;
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();

    int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Optimal solution found" << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
        cout << "Time limit reached, no optimal solution found within the specified time." << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else {
        cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    }
    cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    return obj_gurobi;
}
int MKDCMIP::MIPSolve2(){
    GRBEnv env=GRBEnv(true);
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    // model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    // model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    //add 2hop constraint
    vector<ui> cmm_neighbors;
    for (int u = 0; u < n; u++){
        for (int v = u+1; v < n; v++){
            if(dist(u,v)==2){
                getCN(cmm_neighbors, u, v);
                GRBLinExpr constr=0;
                constr+=(x[u]+x[v]);
                for(auto cmm_nei:cmm_neighbors){
                    constr-=x[cmm_nei];
                }
                model.addConstr(constr<=1, "d_"+to_string(u)+to_string(v));
            }else if(dist(u,v)==3){
                model.addConstr(x[u]+x[v]<=1, "d_"+to_string(u)+to_string(v));
            }
        }
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+to_string(j));
        }
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d\n", lb, ub);
    for (int k = lb; k <= ub; k++){
        eNumLB+= k*(k-1)*z[k]/2;
    }
    eNumLB-=s;
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();

    int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Optimal solution found" << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
        cout << "Time limit reached, no optimal solution found within the specified time." << endl;
        cout << "Obj: " << obj_gurobi << endl;
    } else {
        cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    }
    cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    return obj_gurobi;
}
int MKDCMIP::MIPSolve2hop(){
    // lb=1;
    GRBEnv env=GRBEnv(true);
    
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
    // model.getEnv().set(GRB_IntParam_ResultFlag, 1); 
    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+"_"+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        // if(edges.first>)
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d, n: %d\n", lb, ub, n);
    for (int k = lb; k <= ub; k++){
        eNumLB+= k*(k-1)*z[k]/2;
    }
    eNumLB-=s;
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if(status==GRB_OPTIMAL){
        int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
        feasible++;
        return obj_gurobi;
    }
    // // Print the status
    // switch (status) {
    //     case GRB_OPTIMAL:
    //         std::cout << "Optimal solution found." << std::endl;
    //         break;
    //     case GRB_INF_OR_UNBD:
    //         std::cout << "Model is infeasible or unbounded." << std::endl;
    //         break;
    //     case GRB_INFEASIBLE:
    //         std::cout << "Model is infeasible." << std::endl;
    //         break;
    //     case GRB_UNBOUNDED:
    //         std::cout << "Model is unbounded." << std::endl;
    //         break;
    //     case GRB_TIME_LIMIT:
    //         std::cout << "Time limit reached." << std::endl;
    //         break;
    //     case GRB_ITERATION_LIMIT:
    //         std::cout << "Iteration limit reached." << std::endl;
    //         break;
    //     case GRB_NODE_LIMIT:
    //         std::cout << "Node limit reached." << std::endl;
    //         break;
    //     // Add more cases as needed
    //     default:
    //         std::cout << "Optimization ended with status " << status << "." << std::endl;
    //         break;
    // }
    // model.computeIIS();
    // model.write("model.ilp");
    // int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
    // if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    //     cout << "Optimal solution found" << endl;
    //     cout << "Obj: " << obj_gurobi << endl;
    // } else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    //     cout << "Time limit reached, no optimal solution found within the specified time." << endl;
    //     cout << "Obj: " << obj_gurobi << endl;
    // } else {
    //     cout << "Optimization was terminated with status " << model.get(GRB_IntAttr_Status) << endl;
    // }
    // cout << "Gurobi Time: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        // cout << "Total Weight: " << obj_gurobi + weight_fixed_ << endl;
    // cout<<11111<<endl;
    // return obj_gurobi;
    return 0;
}
float MKDCMIP::getUB2hop(){
    // lb=1;
    GRBEnv env=GRBEnv(true);
    
    env.set("LogFile", "opt.log");
    env.start();
    GRBModel model=GRBEnv(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    // set number of threads
    model.set(GRB_IntParam_Threads, 1);

    // Set time limit
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
    // model.getEnv().set(GRB_IntParam_ResultFlag, 1); 
    // set tolerance
    model.set(GRB_DoubleParam_MIPGap, 0);
    
    //2. init the vars
    GRBVar *x=new GRBVar[n];
    GRBVar **e=NULL;
    e=new GRBVar*[n];
    for (int i = 0; i < n; i++){
        e[i]=new GRBVar[n];
    }
    GRBVar *z=new GRBVar[n+1];
    // GRBVar *k=new GRBVar[n];
    GRBLinExpr objective = 0;
    //init x[i]
    for (int v = 0; v < n; v++) {
        x[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x_" + to_string(v));
        objective += x[v];
    }
    model.setObjective(objective, GRB_MAXIMIZE);
    //init e[i,j]
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            e[i][j]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_"+to_string(i)+"_"+to_string(j));
        }
        
    }
    //init z[i]
    for (int i = 0; i < n+1; i++){
        z[i]=model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_"+to_string(i));
    }

    //add constraints
    // add constraints (6b)
    GRBLinExpr eNum=0; GRBLinExpr eNumLB=0;
    for(auto edge: edges){
        // if(edges.first>)
        eNum+=e[edge.first][edge.second];
    }
    // printf("lb: %d, ub: %d, n: %d\n", lb, ub, n);
    for (int k = lb; k <= ub; k++){
        eNumLB+= k*(k-1)*z[k]/2;
    }
    eNumLB-=s;
    model.addConstr(eNum>=eNumLB, "c0");

    // add constraints (6c)
    for(auto edge: edges){
        int i=edge.first, j=edge.second;
        model.addConstr(e[i][j]<= x[i],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(i));
        model.addConstr(e[i][j]<= x[j],"c_"+to_string(i)+"_"+to_string(j)+"_"+to_string(j));
    }

    //add constraints (6d)
    GRBLinExpr vNum1=0, vNum2=0; GRBLinExpr zSum=0;
    for (int v = 0; v < n; v++)
    {
        vNum1+=x[v];
    }
    for (int k = lb; k <= ub; k++)
    {
        vNum2+=k*z[k];
        zSum+=z[k];
    }
    model.addConstr(vNum1==vNum2, "c1");
    model.addConstr(zSum==1, "c2");

    //optimize
    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if(status==GRB_OPTIMAL){
        int obj_gurobi = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal));
        feasible++;
        return obj_gurobi;
    }
    return 0.0;
}
MKDCMIP::~MKDCMIP()
{
    if(!edges.empty()){
        edges.clear();
    }
}





#endif