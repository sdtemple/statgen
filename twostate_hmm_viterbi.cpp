// Viterbi Training for 2-state Hidden Markov Model
// (n) Seth Temple (e) sdtemple at uw dot edu

#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<algorithm>
#include<vector>
#include<math.h>
#include<limits>
#include<unordered_map>
using namespace std;

struct vertex {
    int state;
    double wt = -numeric_limits<double>::infinity();
    vertex* back;
};

struct edge {
    vertex* from;
    vertex* to;
    string nucl;
};

void update_wt(edge e, double** A, double** B,
    unordered_map<string, int> b);
string viterbi_parse(vertex v);
void update_model(string parse, string obs, double** A, double** B);
void reset_wts(vertex vs[], int vlen);
void viterbi_train(vertex vs[], edge es[], string obs,
    double** A, double** B, unordered_map<string, int> b,
    int vlen, int elen);
void cpg_segment(string parse);

int main()
{
    // initialization
    double** A = new double* [4];
    for (int i = 0; i < 4; i++) A[i] = new double[4];

    // log transition matrix
    A[0][0] = log10(.999);
    A[0][1] = log10(.001);
    A[1][0] = log10(.01);
    A[1][1] = log10(.99);

    // log first and last step transition probs
    A[2][0] = log10(.996);
    A[2][1] = log10(.004);
    A[0][3] = 0;
    A[1][3] = 0;

    unordered_map<string, int> b;
    b["A"] = 0;
    b["C"] = 1;
    b["G"] = 2;
    b["T"] = 3;
    b["END"] = 4;

    double** B = new double* [5];
    for (int i = 0; i < 4; i++) B[i] = new double[5];
    B[0][0] = log10(.3); B[0][3] = log10(.3);
    B[0][1] = log10(.2); B[0][2] = log10(.2);
    B[0][4] = 0;
    B[1][0] = log10(.15); B[1][3] = log10(.15);
    B[1][1] = log10(.35); B[1][2] = log10(.35);
    B[1][4] = 0;
    B[3][4] = 0;

    // I/O
    ifstream fasta("../../../D/Pyrococcus_horikoshii.fasta", ios::in);
    string title;
    getline(fasta, title);

    string obs;
    getline(fasta, obs); fasta.close();
    int len = obs.length();

    // build HMM
    int vlen = 2 * len + 2;
    int elen = 4 * len;
    vertex* vs = new vertex[vlen];
    edge* es = new edge[elen];

    vs[0].state = 2;
    vs[0].back = &vs[0];
    vs[0].wt = 0;
    vs[1].state = 0;
    vs[1].back = &vs[1];
    vs[2].state = 1;
    vs[2].back = &vs[2];

    es[0].from = &vs[0];
    es[0].to = &vs[1];
    es[0].nucl = obs[0];
    es[1].from = &vs[0];
    es[1].to = &vs[2];
    es[1].nucl = obs[0];

    for (int i = 1; i < len; i++) {
        // vertices
        vs[2 * i + 1].state = 0;
        vs[2 * i + 1].back = &vs[2 * i + 1];

        vs[2 * i + 2].state = 1;
        vs[2 * i + 2].back = &vs[2 * i + 2];

        // edges
        es[4 * i - 2].from = &vs[2 * (i - 1) + 1];
        es[4 * i - 2].to = &vs[2 * i + 1];
        es[4 * i - 2].nucl = obs[i];

        es[4 * i - 1].from = &vs[2 * (i - 1) + 1];
        es[4 * i - 1].to = &vs[2 * i + 2];
        es[4 * i - 1].nucl = obs[i];

        es[4 * i].from = &vs[2 * i];
        es[4 * i].to = &vs[2 * i + 1];
        es[4 * i].nucl = obs[i];

        es[4 * i + 1].from = &vs[2 * i];
        es[4 * i + 1].to = &vs[2 * i + 2];
        es[4 * i + 1].nucl = obs[i];
    }

    vs[vlen - 1].state = 3;
    vs[vlen - 1].back = &vs[vlen - 1];

    es[elen - 2].from = &vs[vlen - 3];
    es[elen - 2].to = &vs[vlen - 1];
    es[elen - 2].nucl = "END";
    es[elen - 1].from = &vs[vlen - 2];
    es[elen - 1].to = &vs[vlen - 1];
    es[elen - 1].nucl = "END";

    for (int i = 1; i < 11; i++) {
        viterbi_train(vs, es, obs, A, B, b, vlen, elen);
    }

    for (int i = 0; i < elen; i++) update_wt(es[i], A, B, b);
    string parse = viterbi_parse(vs[vlen - 1]);
    cpg_segment(parse);

    // deallocate
    for (int i = 0; i < 4; i++)
        delete[] B[i];
    delete[] B;

    for (int i = 0; i < 4; i++)
        delete[] A[i];
    delete[] A;

    return(0);
}

void update_wt(edge e, double** A, double** B,
    unordered_map<string, int> b) {
    int i = (e.from)->state;
    int j = (e.to)->state;
    int k = b.at(e.nucl);
    double wt = A[i][j] + B[j][k] + (e.from)->wt;
    if (wt > (e.to)->wt) {
        (e.to)->wt = wt;
        (e.to)->back = e.from;
    }
}

string viterbi_parse(vertex v) {
    string parse = "";
    vertex* currentv = &v;
    while (currentv != currentv->back) {
        parse.append(to_string((currentv->back)->state));
        currentv = currentv->back;
    }
    reverse(parse.begin(), parse.end());
    parse = parse.substr(1, parse.length());
    return(parse);
};


void update_model(string parse, string obs, double** A, double** B) {
    double a00 = 0, a01 = 0, a10 = 0, a11 = 0;
    double b0AT = 0, b1AT = 0, b0CG = 0, b1CG = 0;
    double s0 = 0, s1 = 0;
    int eol = obs.length() - 1;
    for (int i = 0; i < eol; i++) {
        // emissions
        if (parse[i] == '0') {
            s0++;
            if (obs[i] == 'A' || obs[i] == 'T') b0AT++;
            else b0CG++;
        }
        else {
            s1++;
            if (obs[i] == 'A' || obs[i] == 'T') b1AT++;
            else b1CG++;
        }
        // transitions
        if (parse[i] == '0' && parse[i + 1] == '0') a00++;
        else if (parse[i] == '0' && parse[i + 1] == '1') a01++;
        else if (parse[i] == '1' && parse[i + 1] == '0') a10++;
        else a11++;
    }

    // transitions
    double A00 = a00 / s0;
    double A01 = a01 / s0;
    double A10 = a10 / s1;
    double A11 = a11 / s1;

    A[0][0] = log10(A00);
    A[0][1] = log10(A01);
    A[1][0] = log10(A10);
    A[1][1] = log10(A11);

    // emissions
    if (parse[eol] == '0') {
        s0++;
        if (obs[eol] == 'A' || obs[eol] == 'T') b0AT++;
        else b0CG++;
    }
    else {
        s1++;
        if (obs[eol] == 'A' || obs[eol] == 'T') b1AT++;
        else b1CG++;
    }
    double B0AT = b0AT / s0;
    double B0CG = b0CG / s0;
    double B1AT = b1AT / s1;
    double B1CG = b1CG / s1;

    // zeroeth state emissions
    B[0][0] = B[0][3] = log10(B0AT / 2);
    B[0][1] = B[0][2] = log10(B0CG / 2);

    // first state emissions
    B[1][0] = B[1][3] = log10(B1AT / 2);
    B[1][1] = B[1][2] = log10(B1CG / 2);
}

void reset_wts(vertex vs[], int vlen) {
    for (int i = 1; i < vlen; i++)
        vs[i].wt = -numeric_limits<double>::infinity();
}

void viterbi_train(vertex vs[], edge es[], string obs,
    double** A, double** B, unordered_map<string, int> b,
    int vlen, int elen) {
    // dynamic programming
    for (int i = 0; i < elen; i++) update_wt(es[i], A, B, b);
    string parse = viterbi_parse(vs[vlen - 1]);
    update_model(parse, obs, A, B);
    reset_wts(vs, vlen);
}

void cpg_segment(string parse) {
    cout << "Segment List:" << endl;
    int eol = parse.length();
    int ct = 0;
    while (ct != eol) {
        if (parse[ct] == '1') {
            cout << ct + 1 << " ";
            while ((parse[ct] == '1') && (ct != eol - 1)) ct++;
            cout << ct << endl;
        }
        ct++;
    }
}


