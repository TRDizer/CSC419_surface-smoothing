#include "massmatrix.h"
#include <cmath> 

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  // New to find all faces that contain vertex i so that M_ii can be assembled properly with the 
  // correct summation value. Can be done  getting collecting Area(t_i) on each face visit individually
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(F.maxCoeff() + 1);
  double ops, adj1, adj2, area;
  
  for (int f = 0; f < F.rows(); f++){
  	ops = l(f, 0);
  	adj1 = l(f, 1);
  	adj2 = l(f, 2);
    area = std::sqrt((ops+adj1+adj2)*(-ops+adj1+adj2)*(ops-adj1+adj2)*(ops+adj1-adj2)) / 4.0;
  	
    temp(F(f, 0)) += area;
    temp(F(f, 1)) += area;
    temp(F(f, 2)) += area;
  }

  temp /= 3.0;
  M = temp.asDiagonal();
}

