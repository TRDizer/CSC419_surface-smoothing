#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = G;
  // Mt Vt = (Mt - lambda * Lt) Vt+1
  Eigen::MatrixXd edge_length;
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> mass;
  Eigen::SparseMatrix<double> cots;

  igl::edge_lengths(V, F, edge_length);
  cotmatrix(edge_length, F, cots);
  massmatrix(edge_length, F, mass);

  // Eigen::MatrixXd A = -lambda * cots;
  // A.diagonal() += mass;   // No go

  Eigen::SparseMatrix<double> A;
  A = -lambda * cots;
  for (int i = 0; i < V.rows(); i++){
  	A.coeffRef(i, i) += mass.diagonal()(i);
  }

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solve;
  solver.compute(A);
  // solver.compute(A.sparseView());   // No go
  U = solver.solve(mass * G);
}
