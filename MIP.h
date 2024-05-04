#ifndef _MIP_H_
#define _MIP_H_
#include <gurobi_c++.h>
#include<vector>
#include "Utility.h"
using namespace std;
int feasible=0, unfeasible=0;

class MIP;

class MIP
{
private:
    int n,m;
    int lb, ub;
    double gamma;
    vector<pair<int,int>> edges;
public:
    MIP(/* args */);
    // MIP(int _n, double _gamma, int _lb);
    ~MIP();
    void load_graph(int _n, int *_pstart, int *_pend, int* _edges,int _lb, double _gamma);
    int MIPSolve();
    int MIPSolve2hop();
    void load_subgraph(int _n, vector<pair<int, int>>& _vp, int _lb, double _gamma);
    void printinfo();
};

MIP::MIP(/* args */)
{
    n=m=0;
    lb=0, ub=INT32_MAX;
    gamma=1.0;
}
void MIP::load_graph(int _n, int *_pstart, int *_pend, int *_edges,int _lb, double _gamma){
    n=_n; lb=_lb; gamma=_gamma;
    for (int u = 0; u < n; u++){
        for (int j = _pstart[u]; j < _pend[u]; j++){
            int v=_edges[j];
            if(u<v) edges.push_back(make_pair(u,v)), m++;
        }
    }
    ub=floor(0.5+0.50*sqrt(1+8.0*m/gamma));
    ub=min(n, ub);
}
void MIP::load_subgraph(int _n, vector<pair<int, int>>& _vp, int _lb, double _gamma){
    n=_n; lb=_lb; gamma=_gamma;
    maxSubSz=max(maxSubSz, n);
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
void MIP::printinfo(){
    printf("vertex num: %d, edge num: %d, density: %.2f\n", n, m, 2.0*m/(n*(n-1)));
}

int MIP::MIPSolve(){
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
int MIP::MIPSolve2hop(){
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
MIP::~MIP()
{
    if(!edges.empty()){
        edges.clear();
    }
}







#endif