// Demonstrate how to use the sparse matrix class and solve the linear system
// Compilation:
//   g++ -o demo demo.cpp solver.cpp
//
// Solves three matrices: 3x3, 20x20, and a big one 400x400


#include <cstdio>
#include <cstdlib>
#include <valarray>
using namespace std;

// this is the header file of solver class
#include "solver.h"

// Solves x for a small 3*3 sparse matrix a, write the matrix a and b directly in code
void solve_small() {
  cout << "** small demonstration" << endl;
  coo_matrix A;                                 // A ----> Object of coo_matrix type.
  int R[]    = {0, 0, 1, 1, 1, 2, 2};          //Temporarily stores row elements of A
  int C[]    = {0, 1, 0, 1, 2, 1, 2};          //Temporarily stores col elements of A
  double V[] = {4.0, -1.0, -1.0,  4.0, -1.0, -1.0, 4.0};  //Temporarily stores dat elements of A..... NOTE: Here, in this example all the pad nets have a weight of 3. (This can be deduced by the numbers in the sparse matrix)
  A.n = 3;              
  A.nnz = sizeof(R) / sizeof(int);     //Number of non zero elements is the number of elements in array R
  A.row.resize(A.nnz);                 //Fixing size of valarrays row, col and dat as nnz
  A.col.resize(A.nnz);
  A.dat.resize(A.nnz);
  
  A.row = valarray<int>(R, A.nnz);         //Copying row, col, dat of A from temporary R,C and V
  A.col = valarray<int>(C, A.nnz);
  A.dat = valarray<double>(V, A.nnz);

  // initialize as [1, 1, 1]
  valarray<double> x(1.0, A.n);         //Assume all 3 gates at coordinate x=1 (simply for demonstration purpose)
  valarray<double> b(A.n);              //Size of b column matrix is n--->number of gates 
  A.matvec(x, b); // b = A*x is done

  cout << "b should equal [3,2,3]" << endl;
  print_valarray(b);

  // solve for x...   Note: Since we assumed all gates at x=1 and then found B, now if we do the reverse process of taking this B and then solving for x in Ax=B , we will get x=[1,1,1]
  cout << "x = " << endl;
  A.solve(b, x);
  print_valarray(x);
}

// Solves x for a 20*20 sparse matrix a, read matrix a and b from txt files
void solve_psd() {
  cout << "** PSD demonstration" << endl;
  coo_matrix psd;
  psd.read_coo_matrix("psd.txt");   //read sparse matrix non zero elements from "psd.txt"
  valarray<double> x(psd.n);        //Create an array of length n to implement x matrix
  valarray<double> b(psd.n);        //Create an array of length n to implement b matrix

  
  ifstream ifs("b.txt");           //read sparse matrix-->vector ( b ) elements from "b.txt"
  for (int i = 0; i < psd.n; ++i) {
    ifs >> b[i];
  }
  ifs.close();

  cout << "b = " << endl;
  print_valarray(b);

  psd.solve(b, x);
  cout << "x = " << endl;
  print_valarray(x);
}

/*
// // Solves x for a 100*100 sparse matrix a, read matrix a from txt file, and for demonstration purpose we have set b matrix to be all 1s
void solve_big() {
  cout << "** big demonstration" << endl;       //Working similar to solve_psd
  coo_matrix big;
  big.read_coo_matrix("../data/mat_helmholtz.txt");
  valarray<double> x(big.n);
  
  // set the matrix b as ones;
  valarray<double> b(1.0, big.n);
  big.solve(b, x);
  
  cout << "x = " << endl;
  print_valarray(x);
}

*/

int main(int argc, char *argv[]) {
  solve_small();
  solve_psd();
  //solve_big();*/
  return 0;
}
