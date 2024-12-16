// Genome Matching via Suffix Arrays
// (n) Seth Temple (e) sdtemple at uw dot edu

#define _CRT_SECURE_NO_WARNINGS
#include<algorithm>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<cstring>
#include<string.h>
using namespace std;

// nucleotide class
struct nucl
{
    char* sfxarry; // suffix array
    int mstrlen; // length of string match
    int genome;
    const char* strand;
    int pos;
};

// user-defined functions
int sfxarry_compare(const void* a, const void* b);
int mstrlen_compare(const void* a, const void* b);
int mstrlen(const char* a, const char* b);
int nmaxlen(nucl a[]);
char ncomplement(const char a);

int main()
{
    time_t start, end;
    time(&start);

    // I/O
    ifstream g1("BA000007.3.fasta", ios::in); // first genome
    ifstream g2("AL935263.2.fasta", ios::in); // second genome

    // allocate enough space to heap
    char* h1 = new char[12]; char* h2 = new char[12];
    char* x = new char[17613706]; char* y = new char[3308274];
    char* w = new char[5498579]; char* z = new char[3308274]; // reverse complement strand

    // reading
    g1.getline(h1, 12); g2.getline(h2, 12);
    g1.getline(x, 5498579); g2.getline(y, 3308274);
    g1.close(); g2.close();

    int n = strlen(x);
    int m = strlen(y);

    // reverse complements
    int ct = 0;
    for (int i = (m - 1); i >= 0; i--)
    {
        if (ncomplement(y[i]))
        {
            z[ct] = ncomplement(y[i]);
            ct++;
        }
    }
    z[ct] = '\0';

    ct = 0;
    for (int i = (n - 1); i >= 0; i--)
    {
        if (ncomplement(x[i]))
        {
            w[ct] = ncomplement(x[i]);
            ct++;
        }
    }
    w[ct] = '\0';

    // nucleotide summary for first genome
    int A, C, G, T, NN;
    A = 0; C = 0; G = 0; T = 0; NN = 0;
    char f;
    for (int i = 0; i < n; i++)
    {
        f = x[i];
        switch (f)
        {
        case 'A':
            A++; break;
        case 'C':
            C++;  break;
        case 'G':
            G++; break;
        case 'T':
            T++; break;
        default:
            NN++; break;
        }
    }

    // nucleotide summary for second genome
    int a, c, g, t, nn;
    a = 0; c = 0; g = 0; t = 0; nn = 0;
    for (int i = 0; i < m; i++)
    {
        f = y[i];
        switch (f)
        {
        case 'A':
            a++; break;
        case 'C':
            c++; break;
        case 'G':
            g++; break;
        case 'T':
            t++; break;
        default:
            nn++; break;
        }
    }

    // concatenations
    strcat(x, "*"); strcat(x, w);
    strcat(x, "*"); strcat(x, y);
    strcat(x, "*"); strcat(x, z);

    // initialize
    nucl* X = new nucl[17613706];
    int N = 2 * n + 2 * m + 4;
    for (int i = 0; i < N; i++)
    {
        X[i].sfxarry = &x[i];
        if (i < (n + 1))
        {
            X[i].genome = 1;
            X[i].strand = "forward";
            X[i].pos = i + 1;
        }
        else if (i < (2 * n + 2))
        {
            X[i].genome = 1;
            X[i].strand = "reverse";
            X[i].pos = 2 * (n + 1) - i;
        }
        else if (i < (2 * n + m + 3))
        {
            X[i].genome = 2;
            X[i].strand = "forward";
            X[i].pos = i - 2 * n - 1;
        }
        else
        {
            X[i].genome = 2;
            X[i].strand = "reverse";
            X[i].pos = 2 * (n + m + 2) - i;
        }
    }

    // sort suffix arrays lexicographically
    qsort(X, N, sizeof(nucl), sfxarry_compare);

    // calculate lengths of string matches
    int l, u, k;
    bool b, B;
    for (int i = 0; i < N; i++)
    {
        // find best neighbor down the list
        k = 0;
        B = 1;
        b = 1;
        while (B && b)
        {
            if (!(k + i - N + 1)) b = 0; // end reached
            else
            {
                k++;
                B = !(X[i].genome - X[i + k].genome);
            }
        }
        if (!b)
        {
            u = 0;
        }
        else u = mstrlen(X[i].sfxarry, X[i + k].sfxarry);

        // find best neighbor up the list
        k = 0;
        B = 1;
        b = 1;
        while (B && b)
        {
            if (!(k + i)) b = 0; // end reached
            else
            {
                k--;
                B = !(X[i].genome - X[i + k].genome);
            }
        }
        if (!b)
        {
            l = 0;
        }

        else l = mstrlen(X[i].sfxarry, X[i + k].sfxarry);

        X[i].mstrlen = max(u, l);
    }

    // sort based on length of string match
    qsort(X, N, sizeof(nucl), mstrlen_compare);
    int maxlen = nmaxlen(X);
    qsort(X, maxlen, sizeof(nucl), sfxarry_compare);

    for (int i = 0; i < N; i++)
    {
        if (X[i].strand == "reverse") X[i].pos = X[i].pos - X[i].mstrlen;
    }
    
    // deallocate heap variables
    delete[] w; delete[] x; delete[] y; delete[] z;
    delete[] h1; delete[] h2;

    return(0);
}

// suffix array comparator
int sfxarry_compare(const void* a, const void* b)
{
    char* p = ((struct nucl*)a)->sfxarry; // access the suffix array
    char* q = ((struct nucl*)b)->sfxarry;
    return(strcmp(p, q)); // compare suffix arrays
}

// match length comparator
int mstrlen_compare(const void* a, const void* b)
{
    int i = ((struct nucl*)a)->mstrlen; // access the match length
    int j = ((struct nucl*)b)->mstrlen;
    return(j - i); // compare based on match length
}

// length of string match
int mstrlen(const char* a, const char* b)
{
    int i = -1; // iterator
    bool l; // logic
    do {
        i++;
        if ((a[i] == '*') & (b[i] == '*'))
        {
            i--;
            break;
        }
        else
        {
            l = !(a[i] - b[i]);
        }
    } while (l);
    return(i);
}

// # of matches of max length
int nmaxlen(nucl a[])
{
    int i = 1;
    int maxlen = a[0].mstrlen;
    while (maxlen == a[i].mstrlen) i++;
    return(i);
}

// nucleotide complement
char ncomplement(char a)
{
    switch (a)
    {
    case 'A':
        return('T');
    case 'T':
        return('A');
    case 'G':
        return('C');
    case 'C':
        return('G');
    default:
        cout << "Sequence error" << endl;
        return(0);
    }
}