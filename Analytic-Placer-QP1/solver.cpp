#include "solver.h"

void coo_matrix::read_coo_matrix(const char *fname){    //Function to read from a txt file and store row, col and dat variables of object of this class.
  // matrix format:
  //  n nnz                NOTE: n is total number of gates in netlist .. nnz is number of non zero elements in our sparse matrix. 
  //  i1 j1 A[i1,j1]       NOTE: One line for each non zero element. i1 is row number, j1 is column number, A[i1, j1] is the number present in the sparse matrix at row i1 and column j1.
  //  i2 j2 A[i2,j2]       This functions reads line by line according to the above format from the txt file and stores in object in 3 valarrays(one for row, one for column and one for non zero element in that respective row and column)
  //    .  .
  //    .  .
  ifstream fin(fname);
  fin >> n >> nnz;
  row.resize(nnz);
  col.resize(nnz);
  dat.resize(nnz);
  for(int i=0; i<nnz; ++i){
    fin >> row[i] >> col[i] >> dat[i];
  }
}

double dot(const valarray<double> &x, const valarray<double> &y){   //Function used in some computation in solve function. Not important to understand.
  return (x * y).sum();
}       

void coo_matrix::matvec(const valarray<double> &x, valarray<double> &y){
  y = 0.0; // need to reset to 0 first.
  // y = A * x
  for(int i = 0; i < nnz; ++i){
    y[row[i]] += dat[i] * x[col[i]];
  }
}

void coo_matrix::solve(const valarray<double> &b, valarray<double> &x){
  // x = A^{-1} b with CG
  // https://en.wikipedia.org/wiki/Conjugate_gradient#Example_code_in_Matlab

  int maxit = 1000;
  valarray<double> Ax(n);
  valarray<double> Ap(n);
  valarray<double> r(n);
  valarray<double> p(n);
  double rnormold, alpha, rnorm;
  double error, errorold = 1.0;

  for(size_t i=0; i<x.size(); ++i){
  	x[i] = (double)rand()/(double)RAND_MAX;
  }
  
  matvec(x,Ax);
  r = b - Ax;
  p = r;
  rnormold = dot(r, r);

  int i;
  for(i=0; i<maxit; ++i){
    matvec(p,Ap);
    alpha = rnormold / dot(p, Ap);

    // p *= alpha;
    x += alpha * p;

    // Ap *= alpha;
    r -= alpha * Ap;

    rnorm = dot(r,r);
    if(sqrt(rnorm) < 1e-8){ break; }
    else { 
      error = abs( dot( r, x) );
      // clog << "||e[" << i << "]||_A = " << error 
            // << "     ratio = " << error/errorold << endl;
      errorold = error;
    }

    p *= (rnorm/rnormold);
    p += r;

    rnormold = rnorm;
  }

  if ( i == maxit )
    cerr << "Warning: reaches maximum iteration." << endl;
}