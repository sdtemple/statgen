#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
using namespace std;

size_t amino_id(char a);
string amino_iter(string line);
int amino_wt(string line);

int main()
{
    int matrix[25][25];
    int gap = -6;
    string line;
    ifstream blosum("blosum62.txt", ios::in);
    getline(blosum, line);
    for (int l = 0; l < 24; l++) {
        getline(blosum, line);
        for (int k = 0; k < 24; k++) {
            line = amino_iter(line);
            matrix[l][k] = amino_wt(line);
        }
        matrix[l][24] = gap;
    }

    for (int k = 0; k < 25; k++) {
        matrix[24][k] = gap;
    }

    return(0);
}

size_t amino_id(char a) {
    const string aa = "ARNDCQEGHILKMFPSTWYVBZX*-";
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