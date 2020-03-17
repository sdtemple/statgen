// Rectangular Prism Edit Graph for Multiple Sequence Alignment
// (n) Seth Temple (e) sdtemple at uw dot edu

#include<sstream>
#include<fstream>
#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
using namespace std;

// substitution matrix
size_t amino_id(char a);
string amino_iter(string line);
int amino_wt(string line);

// edit graph

// utility
int pw_score(char a, char b, int** m);
int trp_score(string& e, int** m);
int incube(char a, char b, char c, int*** cube);

// axes
void xcoord(string& e, int i,
    ofstream& vs, ofstream& es, int** m);
void ycoord(string& e, int j,
    ofstream& vs, ofstream& es, int** m);
void zcoord(string& e, int k, ofstream& vs, ofstream& es, int** m);
void xaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);
void yaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);
void zaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);

// planes
void xycoord(string& e1, string& e2, string& e3,
    int i, int j,
    ofstream& vs, ofstream& es,
    int** m);
void xzcoord(string& e1, string& e2, string& e3,
    int i, int k,
    ofstream& vs, ofstream& es,
    int** m);
void yzcoord(string& e1, string& e2, string& e3,
    int j, int k,
    ofstream& vs, ofstream& es,
    int** m);
void xyplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);
void xzplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);
void yzplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube);

// 3D
void xyzcoord(string& e1, string& e2, string& e3,
    string& e4, string& e5, string& e6, string& e7,
    int i, int j, int k,
    ofstream& vs, ofstream& es,
    int** m);
void xyzspace(string& l1, string& l2, string& l3,
    vector<string>& edges, ofstream& vs, ofstream& es,
    int** m, int*** cube);

int main()
{
    // I/O
    ifstream prot1("../../../D/P02818.txt", ios::in);
    ifstream prot2("../../../D/A2D670.txt", ios::in);
    ifstream prot3("../../../D/XP_029364794.1.txt", ios::in);
    ofstream vs("../../vertices.txt", ios::out);
    ofstream es("../../edges.txt", ios::out);
    string line1, line2, line3;
    getline(prot1, line1); getline(prot2, line2); getline(prot3, line3);
    getline(prot1, line1); getline(prot2, line2); getline(prot3, line3);
    prot1.close(); prot2.close(); prot3.close();

    // substitution matrix

    // initialization
    int** m = new int* [21];
    for (int i = 0; i < 21; i++) {
        m[i] = new int[21];
    }
    int gap = -6;
    string line;
    ifstream blosum("../../../D/blosum62.txt", ios::in);

    // data population
    getline(blosum, line);
    for (int i = 0; i < 20; i++) {
        getline(blosum, line);
        for (int j = 0; j < 20; j++) {
            line = amino_iter(line);
            m[i][j] = amino_wt(line);
        }
        m[i][20] = gap;
    }
    for (int k = 0; k < 20; k++) {
        m[20][k] = gap;
    }
    m[20][20] = 0;
    blosum.close();

    // edit graph

    // initialization
    vector<string> edges;
    int*** cube = new int** [21];
    for (int i = 0; i < 21; i++) {
        cube[i] = new int* [21];
        for (int j = 0; j < 21; j++) {
            cube[i][j] = new int[21];
            for (int k = 0; k < 21; k++) {
                cube[i][j][k] = 0;
            }
        }
    }

    // data population
    vs << "V " << 100020003000 << endl; // origin
    xaxis(line1, edges, vs, es, m, cube);
    yaxis(line2, edges, vs, es, m, cube);
    zaxis(line3, edges, vs, es, m, cube);
    xyplane(line1, line2, edges, vs, es, m, cube);
    xzplane(line1, line3, edges, vs, es, m, cube);
    yzplane(line2, line3, edges, vs, es, m, cube);
    xyzspace(line1, line2, line3, edges, vs, es, m, cube);

    // closing
    vs.close();
    es.close();

    // output
    ofstream output("../../ancillary.txt", ios::out);

    sort(edges.begin(), edges.end());
    output << "Edge weights:" << endl;
    for (auto e : edges) {
        output << e << " = " << trp_score(e, m) << endl;
    }
    output << endl;
    output << "Edge counts:" << endl;
    for (auto e : edges) {
        output << e << " = ";
        output << cube[amino_id(e[0])][amino_id(e[1])][amino_id(e[2])] << endl;
    }
    output.close();

    // deallocating

    // substitution matrix
    for (int i = 0; i < 21; i++) {
        delete[] m[i];
    }
    delete[] m;

    // cube
    for (int i = 0; i < 21; i++) {
        for (int j = 0; j < 21; j++) {
            delete[] cube[i][j];
        }
        delete[] cube[i];
    }
    delete[] cube;

    return(0);
}

// substitution matrix
size_t amino_id(char a) {
    const string aa = "ARNDCQEGHILKMFPSTWYV-";
    size_t i = aa.find(a);
    return(i);
}
string amino_iter(string line) {
    size_t i, j;
    i = line.find(" ");
    j = line.find("  ");
    if (i == j) return(line.substr(i + 2));
    else return(line.substr(i + 1));
}
int amino_wt(string line) {
    size_t i = line.find(" ");
    return(stoi(line.substr(0, i)));
}

// edit graph

// utility
int pw_score(char a, char b, int** m) {
    size_t i = amino_id(a);
    size_t j = amino_id(b);
    return(m[i][j]);
}
int trp_score(string& e, int** m) {
    int p1, p2, p3;
    p1 = pw_score(e[0], e[1], m);
    p2 = pw_score(e[0], e[2], m);
    p3 = pw_score(e[1], e[2], m);
    return(p1 + p2 + p3);
}
int incube(char a, char b, char c, int*** cube) {
    size_t A = amino_id(a);
    size_t B = amino_id(b);
    size_t C = amino_id(c);
    return(cube[A][B][C]);
}

// axes
void xcoord(string& e, int i,
    ofstream& vs, ofstream& es, int** m) {
    vs << "V " << (1000 + i) << 2000 << 3000 << endl;
    es << "E " << e << " "; // name
    es << (1000 + i - 1) << 2000 << 3000 << " ";
    es << (1000 + i) << 2000 << 3000 << " ";
    es << trp_score(e, m) << endl;
}
void ycoord(string& e, int j,
    ofstream& vs, ofstream& es, int** m) {
    vs << "V " << 1000 << (2000 + j) << 3000 << endl;
    es << "E " << e << " "; // name
    es << 1000 << (2000 + j - 1) << 3000 << " ";
    es << 1000 << (2000 + j) << 3000 << " ";
    es << trp_score(e, m) << endl;
}
void zcoord(string& e, int k,
    ofstream& vs, ofstream& es, int** m) {
    vs << "V " << 1000 << 2000 << (3000 + k) << endl;
    es << "E " << e << " "; // name
    es << 1000 << 2000 << (3000 + k - 1) << " ";
    es << 1000 << 2000 << (3000 + k) << " ";
    es << trp_score(e, m) << endl;
}
void xaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e;
    for (int i = 0; i < l.length(); i++) {
        // count
        e = l[i]; e.append("-"); e.append("-");
        val = incube(l[i], '-', '-', cube);
        cube[amino_id(l[i])][20][20]++;
        if (!val) edges.push_back(e);
        // write vertex and edge to graph
        xcoord(e, i + 1, vs, es, m);
    }
}
void yaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e;
    for (int j = 0; j < l.length(); j++) {
        // count
        e = l[j]; e.append("-"); e.insert(0, "-");
        val = incube('-', l[j], '-', cube);
        cube[20][amino_id(l[j])][20]++;
        if (!val) edges.push_back(e);
        // write vertex and edge to graph
        ycoord(e, j + 1, vs, es, m);
    }
}
void zaxis(string& l, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e;
    for (int k = 0; k < l.length(); k++) {
        // count
        e = l[k]; e.insert(0, "-"); e.insert(0, "-");
        val = incube('-', '-', l[k], cube);
        cube[20][20][amino_id(l[k])]++;
        if (!val) edges.push_back(e);
        // write vertex and edge to graph
        zcoord(e, k + 1, vs, es, m);
    }
}

// planes
void xycoord(string& e1, string& e2, string& e3,
    int i, int j,
    ofstream& vs, ofstream& es,
    int** m) {
    vs << "V " << (1000 + i) << (2000 + j) << 3000 << endl;
    // first edge
    es << "E " << e1 << " ";
    es << (1000 + i - 1) << (2000 + j) << 3000 << " ";
    es << (1000 + i) << (2000 + j) << 3000 << " ";
    es << trp_score(e1, m) << endl;
    // second edge
    es << "E " << e2 << " ";
    es << (1000 + i - 1) << (2000 + j - 1) << 3000 << " ";
    es << (1000 + i) << (2000 + j) << 3000 << " ";
    es << trp_score(e2, m) << endl;
    // third edge
    es << "E " << e3 << " ";
    es << (1000 + i) << (2000 + j - 1) << 3000 << " ";
    es << (1000 + i) << (2000 + j) << 3000 << " ";
    es << trp_score(e3, m) << endl;
}
void xzcoord(string& e1, string& e2, string& e3,
    int i, int k,
    ofstream& vs, ofstream& es,
    int** m) {
    vs << "V " << (1000 + i) << 2000 << (3000 + k) << endl;
    // first edge
    es << "E " << e1 << " ";
    es << (1000 + i - 1) << 2000 << (3000 + k) << " ";
    es << (1000 + i) << 2000 << (3000 + k) << " ";
    es << trp_score(e1, m) << endl;
    // second edge
    es << "E " << e2 << " ";
    es << (1000 + i - 1) << 2000 << (3000 + k - 1) << " ";
    es << (1000 + i) << 2000 << (3000 + k) << " ";
    es << trp_score(e2, m) << endl;
    // third edge
    es << "E " << e3 << " ";
    es << (1000 + i) << 2000 << (3000 + k - 1) << " ";
    es << (1000 + i) << 2000 << (3000 + k) << " ";
    es << trp_score(e3, m) << endl;
}
void yzcoord(string& e1, string& e2, string& e3,
    int j, int k,
    ofstream& vs, ofstream& es,
    int** m) {
    vs << "V " << 1000 << (2000 + j) << (3000 + k) << endl;
    // first edge
    es << "E " << e1 << " ";
    es << 1000 << (2000 + j - 1) << (3000 + k) << " ";
    es << 1000 << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e1, m) << endl;
    // second edge
    es << "E " << e2 << " ";
    es << 1000 << (2000 + j - 1) << (3000 + k - 1) << " ";
    es << 1000 << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e2, m) << endl;
    // third edge
    es << "E " << e3 << " ";
    es << 1000 << (2000 + j) << (3000 + k - 1) << " ";
    es << 1000 << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e3, m) << endl;
}
void xyplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e1, e2, e3;
    for (int i = 0; i < l1.length(); i++) {
        for (int j = 0; j < l2.length(); j++) {
            // count
            // third edge
            e3 = l2[j]; e3.append("-"); e3.insert(0, "-");
            cube[20][amino_id(l2[j])][20]++;
            // second edge
            e2 = l1[i]; e2 = e2 + l2[j]; e2.append("-");
            val = incube(l1[i], l2[j], '-', cube);
            cube[amino_id(l1[i])][amino_id(l2[j])][20]++;
            if (!val) edges.push_back(e2);
            // first edge
            e1 = l1[i]; e1.append("-"); e1.append("-");
            cube[amino_id(l1[i])][20][20]++;
            // write vertices and edges to graph
            xycoord(e1, e2, e3, i + 1, j + 1, vs, es, m);
        }
    }
}
void xzplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e1, e2, e3;
    for (int i = 0; i < l1.length(); i++) {
        for (int k = 0; k < l2.length(); k++) {
            // count
            // third edge
            e3 = l2[k]; e3.insert(0, "-"); e3.insert(0, "-");
            cube[20][20][amino_id(l2[k])]++;
            // second edge
            e2 = l1[i]; e2 = e2 + l2[k]; e2.insert(1, "-");
            val = incube(l1[i], '-', l2[k], cube);
            cube[amino_id(l1[i])][20][amino_id(l2[k])]++;
            if (!val) edges.push_back(e2);
            // first edge
            e1 = l1[i]; e1.append("-"); e1.append("-");
            cube[amino_id(l1[i])][20][20]++;
            // write vertices and edges to graph
            xzcoord(e1, e2, e3, i + 1, k + 1, vs, es, m);
        }
    }
}
void yzplane(string& l1, string& l2, vector<string>& edges,
    ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e1, e2, e3;
    for (int j = 0; j < l1.length(); j++) {
        for (int k = 0; k < l2.length(); k++) {
            // count
            // third edge
            e3 = l2[k]; e3.insert(0, "-"); e3.insert(0, "-");
            cube[20][20][amino_id(l2[k])]++;
            // second edge
            e2 = l1[j]; e2 = e2 + l2[k]; e2.insert(0, "-");
            val = incube('-', l1[j], l2[k], cube);
            cube[20][amino_id(l1[j])][amino_id(l2[k])]++;
            if (!val) edges.push_back(e2);
            // first edge
            e1 = l1[j]; e1.append("-"); e1.insert(0, "-");
            cube[20][amino_id(l1[j])][20]++;
            // write vertices and edges to graph
            yzcoord(e1, e2, e3, j + 1, k + 1, vs, es, m);
        }
    }
}

// 3D
void xyzcoord(string& e1, string& e2, string& e3,
    string& e4, string& e5, string& e6, string& e7,
    int i, int j, int k,
    ofstream& vs, ofstream& es,
    int** m) {
    vs << "V " << (1000 + i) << (2000 + j) << (3000 + k) << endl;
    // first edge
    es << "E " << e1 << " ";
    es << (1000 + i - 1) << (2000 + j) << (3000 + k - 1) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e1, m) << endl;
    // second edge
    es << "E " << e2 << " ";
    es << (1000 + i - 1) << (2000 + j) << (3000 + k) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e2, m) << endl;
    // seventh edge
    es << "E " << e7 << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k - 1) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e7, m) << endl;
    // third edge
    es << "E " << e3 << " ";
    es << (1000 + i - 1) << (2000 + j - 1) << (3000 + k - 1) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e3, m) << endl;
    // fourth edge
    es << "E " << e4 << " ";
    es << (1000 + i - 1) << (2000 + j - 1) << (3000 + k) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e4, m) << endl;
    // fifth edge
    es << "E " << e5 << " ";
    es << (1000 + i) << (2000 + j - 1) << (3000 + k - 1) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e5, m) << endl;
    // sixth edge
    es << "E " << e6 << " ";
    es << (1000 + i) << (2000 + j - 1) << (3000 + k) << " ";
    es << (1000 + i) << (2000 + j) << (3000 + k) << " ";
    es << trp_score(e6, m) << endl;
}
void xyzspace(string& l1, string& l2, string& l3,
    vector<string>& edges, ofstream& vs, ofstream& es,
    int** m, int*** cube) {
    int val;
    string e1, e2, e3, e4, e5, e6, e7;
    for (int i = 0; i < l1.length(); i++) {
        for (int j = 0; j < l2.length(); j++) {
            for (int k = 0; k < l3.length(); k++) {
                // count
                // sixth edge
                e6 = l2[j]; e6.append("-"); e6.insert(0, "-");
                cube[20][amino_id(l2[j])][20]++;
                // fifth edge
                e5 = l2[j]; e5 = e5 + l3[k]; e5.insert(0, "-");
                cube[20][amino_id(l2[j])][amino_id(l3[k])]++;
                // fourth edge
                e4 = l1[i]; e4 = e4 + l2[j]; e4.append("-");
                cube[amino_id(l1[i])][amino_id(l2[j])][20]++;
                // third edge
                e3 = l1[i]; e3 = e3 + l2[j]; e3 = e3 + l3[k];
                val = incube(l1[i], l2[j], l3[k], cube);
                cube[amino_id(l1[i])][amino_id(l2[j])][amino_id(l3[k])]++;
                if (!val) edges.push_back(e3);
                // seventh edge
                e7 = l3[k]; e7.insert(0, "-"); e7.insert(0, "-");
                cube[20][20][amino_id(l3[k])]++;
                // second edge
                e2 = l1[i]; e2.append("-"); e2.append("-");
                cube[amino_id(l1[i])][20][20]++;
                // first edge
                e1 = l1[i]; e1 = e1 + l3[k]; e1.insert(1, "-");
                cube[amino_id(l1[i])][20][amino_id(l3[k])]++;
                // write vertices and edges to graph
                xyzcoord(e1, e2, e3, e4, e5, e6, e7,
                    i + 1, j + 1, k + 1, vs, es, m);
            }
        }
    }
}