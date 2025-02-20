#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <stdlib.h>
#include <valarray>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

class coo_matrix{
  // COOrdinate sparse matrix
  //Member variables declarations
  public:
    int n;      //No. of gates in our netlist
    int nnz;    //No. of non zero elements in n*n sparse matrix 
    valarray<int> row;      //Valarray is similar to array.. row is a valarray of integer data types
    valarray<int> col;      //col is a valarray of integer data types
    valarray<double> dat;   //dat is a valarray of double data type 


    //Function declarations below. The actual function definitions in solver.cpp ...

    void read_coo_matrix(const char *fname);     //function to read non zero elements of our sparse matrix from file with name under "fname" 
    void matvec(const valarray<double> &x, valarray<double> &y);   //Perfroms matrix multiplication of matrix represented by object of coo_matrix (say a) and x and stores in matrix y ...... y=a*x
    void solve(const valarray<double> &b, valarray<double> &x);    //Solves for x in a*x=b (a,b and c are matrices), where a is represented as an object of coo_matrix class
};


//Print valarray function
template<typename T>
void print_valarray(valarray<T>& v) {
  for (size_t i = 0; i < v.size(); ++i) {
    cout << v[i] << endl;
  }
}

#endif
