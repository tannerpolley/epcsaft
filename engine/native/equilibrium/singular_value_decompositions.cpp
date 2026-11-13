// The SVD classes the equilibrium and regression units use, compiled once; those units declare them extern.
#include <Eigen/Dense>

template class Eigen::JacobiSVD<Eigen::MatrixXd>;
template class Eigen::BDCSVD<Eigen::MatrixXd>;
template class Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
