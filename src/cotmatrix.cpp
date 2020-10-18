#include "cotmatrix.h"
#include <cmath>
#include <iostream>

#define sqr(x)    (x)*(x)
#define get_ops_v1(v)    ((v + 1) % 3)
#define get_ops_v2(v)    ((v + 2) % 3)

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  //        a
  //       /|\ 
  //      / | \  L(a,c) gives (cot(α_ac) + cot(β_ac))/2
  //     /  |  \ 
  //    b\fα|fβ/d
  //      \ | /
  //       \|/  
  //        c 
  // gonna get cot(α_ac)/2 when visting fα add cot(β_ac)/2 ontop when visting fβ
  auto get_cot = [](double ops, double adj1, double adj2) {
    double cos_target = (sqr(adj1) + sqr(adj2) - sqr(ops)) / (2.0 * adj1 * adj2);
    // 2A / (adj1 * adj2)
    double double_area = std::sqrt((ops+adj1+adj2)*(-ops+adj1+adj2)*(ops-adj1+adj2)*(ops+adj1-adj2)) / 2.0;
    double sin_target = double_area / (adj1 * adj2);
    return cos_target / sin_target;
  };

  int num_V = F.maxCoeff() + 1;
  L.resize(num_V, num_V);
  L.setIdentity();

  double ops, adj1, adj2, phi;
  int ops_v1, ops_v2;
  for (int f = 0; f < F.rows(); f++) {
    for (int v = 0; v < 3; v++) {
      ops_v1 = F(f,get_ops_v1(v));
      ops_v2 = F(f,get_ops_v2(v));

      ops = l(f,v);
      adj1 = l(f,get_ops_v1(v));
      adj2 = l(f,get_ops_v2(v));

      phi = get_cot(ops, adj1, adj2) / 2.0;
      L.coeffRef(ops_v1, ops_v2) += phi;
      L.coeffRef(ops_v2, ops_v1) += phi;
    }
  }

  // Diagonal entry does not exist beforehand
  L.diagonal().array() = -L * Eigen::VectorXd::Ones(num_V);
  L.diagonal().array() += 1;

  std::cout << "cot matrix update completed!" << std::endl;
}
