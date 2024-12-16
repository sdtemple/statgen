// Baum-Welch Training with Rescaling for 2-state Hidden Markov Model
// (n) Seth Temple (e) sdtemple at uw dot edu

#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<algorithm>
#include<math.h>
using namespace std;

int emit_state(char o);
double arr_sum(double* arr, int len);
double arr_sum_some(double* arr, double* some, int len);
void fwd_alg(string obs, int len,
    double** f, double* c, double** t, double** e, double* i);
void bck_alg(string obs, int len,
    double** b, double* c, double** t, double** e);
void bw_update(string obs, int len,
    double** gamma, double** psi,
    double* aa, double* cc, double* gg, double* tt,
    double** f, double** b, double** t, double** e, double* i);
double log_likelihood(double* c, int len);

int main()
{
    // I/O
    ifstream fasta("../../../D/Pyrococcus_horikoshii.fasta", ios::in);

    string title;
    getline(fasta, title);

    string obs;
    getline(fasta, obs);
    int len = obs.length();
    fasta.close();

    double* aa = new double[len];
    double* cc = new double[len];
    double* gg = new double[len];
    double* tt = new double[len];
    for (int j = 0; j < len; j++) {
        aa[j] = 0;
        cc[j] = 0;
        gg[j] = 0;
        tt[j] = 0;
    }
    for (int j = 0; j < len; j++) {
        if (obs[j] == 'A') aa[j] = 1;
        else if (obs[j] == 'C') cc[j] = 1;
        else if (obs[j] == 'G') gg[j] = 1;
        else tt[j] = 1;
    }

    // initialization
    double i[2] = { .996, .004 };

    // transition probs
    double** t = new double* [2];
    t[0] = new double[2]; t[1] = new double[2];

    t[0][0] = .999;
    t[0][1] = .001;
    t[1][0] = .01;
    t[1][1] = .99;

    // emission probs
    double** e = new double* [2];
    e[0] = new double[4]; e[1] = new double[4];

    e[0][0] = .3;
    e[0][1] = .2;
    e[0][2] = .2;
    e[0][3] = .3;
    e[1][0] = .15;
    e[1][1] = .35;
    e[1][2] = .35;
    e[1][3] = .15;

    // forward probs
    double** f = new double* [2];
    f[0] = new double[len];
    f[1] = new double[len];

    // backward probs
    double** b = new double* [2];
    b[0] = new double[len];
    b[1] = new double[len];

    // scaling factors
    double* c = new double[len];

    double** gamma = new double* [2];
    gamma[0] = new double[len];
    gamma[1] = new double[len];

    double** psi = new double* [4];
    psi[0] = new double[len - 1];
    psi[1] = new double[len - 1];
    psi[2] = new double[len - 1];
    psi[3] = new double[len - 1];

    // iterate
    fwd_alg(obs, len, f, c, t, e, i);
    bck_alg(obs, len, b, c, t, e);
    bw_update(obs, len, gamma, psi, aa, cc, gg, tt,
        f, b, t, e, i);
    double ll = log_likelihood(c, len), LL;
    int ct = 1;
    do {
        ct += 1;
        LL = ll;
        fwd_alg(obs, len, f, c, t, e, i);
        bck_alg(obs, len, b, c, t, e);
        bw_update(obs, len, gamma, psi, aa, cc, gg, tt,
            f, b, t, e, i);
        ll = log_likelihood(c, len);
    } while ((LL - ll <= -.1) && ct >= 2);

    // deallocating
    delete[] e[0]; delete[] e[1]; delete[] e;
    delete[] t[0]; delete[] t[1]; delete[] t;
    delete[] f[0]; delete[] f[1]; delete[] f;
    delete[] b[0]; delete[] b[1]; delete[] b;
    delete[] psi[0]; delete[] psi[1]; delete[] psi[2]; delete[] psi[3]; delete[] psi;
    delete[] gamma[0]; delete[] gamma[1]; delete[] gamma;
    delete[] aa; delete[] cc; delete[] gg; delete[] tt;
    delete[] c;

    return(0);
}

int emit_state(char o) {
    if (o == 'A') return(0);
    else if (o == 'C') return(1);
    else if (o == 'G') return(2);
    else if (o == 'T') return(3);
    else return(-1); // error
}
double arr_sum(double* arr, int len) {
    double s = 0;
    for (int j = 0; j < len; j++) {
        s += arr[j];
    }
    return(s);
}
double arr_sum_some(double* arr, double* some, int len) {
    double s = 0;
    for (int j = 0; j < len; j++) {
        s += (arr[j] * some[j]);
    }
    return(s);
}
void fwd_alg(string obs, int len,
    double** f, double* c, double** t, double** e, double* i) {
    int s = emit_state(obs[0]);
    f[0][0] = i[0] * e[0][s];
    f[1][0] = i[1] * e[1][s];
    c[0] = 1 / (f[0][0] + f[1][0]);
    f[0][0] = c[0] * f[0][0];
    f[1][0] = c[0] * f[1][0];
    for (int j = 1; j < len; j++) {
        s = emit_state(obs[j]);

        f[0][j] = f[0][j - 1] * t[0][0] * e[0][s] +
            f[1][j - 1] * t[1][0] * e[0][s];
        f[1][j] = f[0][j - 1] * t[0][1] * e[1][s] +
            f[1][j - 1] * t[1][1] * e[1][s];
        c[j] = 1 / (f[0][j] + f[1][j]);
        f[0][j] = c[j] * f[0][j];
        f[1][j] = c[j] * f[1][j];
    }
}
void bck_alg(string obs, int len,
    double** b, double* c, double** t, double** e) {
    int s;
    b[0][len - 1] = c[len - 1];
    b[1][len - 1] = c[len - 1];
    for (int j = len - 2; j >= 0; j--) {
        s = emit_state(obs[j + 1]);

        b[0][j] = c[j] * (b[0][j + 1] * t[0][0] * e[0][s] +
            b[1][j + 1] * t[0][1] * e[1][s]);
        b[1][j] = c[j] * (b[0][j + 1] * t[1][0] * e[0][s] +
            b[1][j + 1] * t[1][1] * e[1][s]);
    }
}
void bw_update(string obs, int len,
    double** gamma, double** psi,
    double* aa, double* cc, double* gg, double* tt,
    double** f, double** b, double** t, double** e, double* i) {

    double g, g0, g1;
    double p, p00, p01, p10, p11;
    int s;
    for (int j = 0; j < len - 1; j++) { // len - 1
        s = emit_state(obs[j + 1]);

        // update gamma
        g0 = b[0][j] * f[0][j];
        g1 = b[1][j] * f[1][j];
        g = 1 / (g0 + g1);
        gamma[0][j] = g * g0;
        gamma[1][j] = g * g1;

        // update psi
        p00 = b[0][j + 1] * e[0][s] * t[0][0] * f[0][j];
        p01 = b[1][j + 1] * e[1][s] * t[0][1] * f[0][j];
        p10 = b[0][j + 1] * e[0][s] * t[1][0] * f[1][j];
        p11 = b[1][j + 1] * e[1][s] * t[1][1] * f[1][j];
        p = 1 / (p00 + p01 + p10 + p11);
        psi[0][j] = p * p00;
        psi[1][j] = p * p01;
        psi[2][j] = p * p10;
        psi[3][j] = p * p11;
    }

    g0 = b[0][len - 1] * f[0][len - 1];
    g1 = b[0][len - 1] * f[1][len - 1];
    g = 1 / (g0 + g1);
    gamma[0][len - 1] = g0 * g;
    gamma[1][len - 1] = g1 * g;


    // update i
    i[0] = gamma[0][0];
    i[1] = gamma[1][0];

    // update t
    double t00, t01, t10, t11, e0, e1;
    e0 = arr_sum(gamma[0], len - 1);
    e1 = arr_sum(gamma[1], len - 1);
    t00 = arr_sum(psi[0], len - 1) / e0;
    t01 = arr_sum(psi[1], len - 1) / e0;
    t10 = arr_sum(psi[2], len - 1) / e1;
    t11 = arr_sum(psi[3], len - 1) / e1;

    // adjustment
    /*
    double T;
    T = t00 + t01;
    t00 = t00 / T;
    t01 = t01 / T;
    T = t10 + t11;
    t10 = t10 / T;
    t11 = t11 / T;
    */

    t[0][0] = t00;
    t[0][1] = t01;
    t[1][0] = t10;
    t[1][1] = t11;

    // update e
    e0 += gamma[0][len - 1];
    e1 += gamma[1][len - 1];
    double aa0, aa1, cc0, cc1, gg0, gg1, tt0, tt1;
    aa0 = arr_sum_some(gamma[0], aa, len) / e0;
    aa1 = arr_sum_some(gamma[1], aa, len) / e1;
    cc0 = arr_sum_some(gamma[0], cc, len) / e0;
    cc1 = arr_sum_some(gamma[1], cc, len) / e1;
    gg0 = arr_sum_some(gamma[0], gg, len) / e0;
    gg1 = arr_sum_some(gamma[1], gg, len) / e1;
    tt0 = arr_sum_some(gamma[0], tt, len) / e0;
    tt1 = arr_sum_some(gamma[1], tt, len) / e1;

    /*
    // adjustment
    double Z = (aa0 + cc0 + gg0 + tt0), O = (aa1 + cc1 + gg1 + tt1);
    aa0 = (aa0 / Z); cout << aa0 << endl;
    cc0 = (cc0 / Z); cout << cc0 << endl;
    gg0 = (gg0 / Z); cout << gg0 << endl;
    tt0 = (tt0 / Z); cout << tt0 << endl;
    aa1 = (aa1 / O); cout << aa1 << endl;
    cc1 = (cc1 / O); cout << cc1 << endl;
    gg1 = (gg1 / O); cout << gg1 << endl;
    tt1 = (tt1 / O); cout << tt1 << endl;
    */

    e[0][0] = aa0; e[0][3] = tt0;
    e[0][1] = cc0; e[0][2] = gg0;
    e[1][0] = aa1; e[1][3] = tt1;
    e[1][1] = cc1; e[1][2] = gg1;

}
double log_likelihood(double* c, int len) {
    double l = 0;
    for (int j = 0; j < len; j++) {
        l -= log2(c[j]);
    }
    return(l);
}
