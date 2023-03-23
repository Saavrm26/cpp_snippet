#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;
vvi adj;
void input(){
    int n,m;
    adj.resize(n+1);
    for(int i=1;i<=m;i++){
        int u,v;
        cin>>u>>v;
        adj[u].emplace_back(v);
    }
}
void dfs(int u,int p){
    for(auto &v :adj[u]){
        if(v==p) continue;
        dfs(v,u);
    }
}