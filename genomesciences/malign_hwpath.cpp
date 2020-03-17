// Uniform Highest Weight Path for Directed Acyclical Graphs
// (n) Seth Temple (e) sdtemple at uw dot edu

#include<stdlib.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<sstream>
#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

struct node
{
    string id;
    double wt = 0;
    node* before;
    string back = "";
};
struct edge
{
    string id;
    double wt;
    node* parent;
    node* child;
};
struct path
{
    vector<string> seq;
    string start;
    string end;
    double score;
};

void update_wt(edge e);
node get_hwpath(node ns[], int size);
path return_hwpath(node n, edge es[], int size);

int main()
{
    // I/O
    ifstream vertices("../../vertices.txt", ios::in);
    ifstream edges("../../edges.txt", ios::in);
    string line;

    // count the nodes and edges
    int nlen = 0, elen = 0;
    while (getline(vertices, line)) nlen++;
    vertices.clear(); vertices.seekg(0, ios::beg);
    while (getline(edges, line)) elen++;
    edges.clear(); edges.seekg(0, ios::beg);

    // populate nodes and edges
    node* ns = new node[nlen]; edge* es = new edge[elen];
    nlen = elen = 0;

    istringstream iss;
    string a, b, c, d, e;
    unordered_map<string, node*> umap;
    while (getline(vertices, line)) {
        istringstream iss(line);
        iss >> a >> b;
        ns[nlen].id = b;
        umap[b] = &ns[nlen];
        ns[nlen].before = &ns[nlen];
        nlen++;
    };
    while (getline(edges, line)) {
        istringstream iss(line);
        iss >> a >> b >> c >> d >> e;
        es[elen].id = b;
        es[elen].wt = stod(e); // string to float
        es[elen].parent = umap.at(c);
        es[elen].child = umap.at(d);
        elen++;
    };
    vertices.close(); edges.close();

    // dynamic programming, unrestricted
    for (int i = 0; i < elen; i++) {
        update_wt(es[i]);
    }
    path up = return_hwpath(get_hwpath(ns, nlen), es, elen);

    delete[] ns; delete[] es;
    return(0);
}

void update_wt(edge e)
{
    double m = (e.parent)->wt + e.wt;
    if (m > (e.child)->wt)
    {
        (e.child)->wt = m;
        (e.child)->back = e.id;
        (e.child)->before = e.parent;
    }
}
node get_hwpath(node ns[], int size)
{
    node n = ns[0];
    for (int i = 1; i < size; i++)
    {
        if (ns[i].wt > n.wt) n = ns[i];
    }
    return(n);
}
path return_hwpath(node n, edge es[], int size)
{
    node* cn = &n;
    path p;
    p.end = cn->id;
    p.score = cn->wt;
    while (cn != cn->before)
    {
        p.seq.push_back(cn->back);
        cn = cn->before;
    }
    p.start = cn->id;
    return(p);
}