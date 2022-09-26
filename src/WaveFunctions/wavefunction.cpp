#include <numeric>
#include "wavefunction.h"

WaveFunction::WaveFunction(System *system)
{
    m_system = system;
}

Eigen::Map<Eigen::VectorXd> WaveFunction::flatten(Eigen::MatrixXd A)
{
    return Eigen::Map<Eigen::VectorXd>(A.data(), A.size());
}

Eigen::Map<Eigen::MatrixXd> WaveFunction::reshape(Eigen::VectorXd A,
                                                  const Eigen::Index m,
                                                  const Eigen::Index n)
{
    return Eigen::Map<Eigen::MatrixXd>(A.data(), m, n);
}

Eigen::Map<Eigen::MatrixXd> WaveFunction::square(Eigen::VectorXd A)
{
    Eigen::Index size = A.size();
    return Eigen::Map<Eigen::MatrixXd>(A.data(), Eigen::Index(sqrt(size)), Eigen::Index(sqrt(size)));
}

Eigen::MatrixXd WaveFunction::stack(Eigen::VectorXd A, int n)
{
    Eigen::MatrixXd stacked = Eigen::MatrixXd::Zero(A.size(), n);
    for (int i = 0; i < n; i++) {
        stacked.row(i) = A;
    }
    return stacked;
}

/* ----------------------------------------------------------------------------
  Sort vector v and also return the mapping
---------------------------------------------------------------------------- */

void WaveFunction::argsort(Eigen::VectorXd &v, Eigen::VectorXi &map)
{

    // reset sorting map
    map.resizeLike(v);
    std::iota(map.data(), map.data() + map.size(), 0);

    // sort v
    std::sort(v.data(),v.data()+v.size());
    
    // sort map by value
    std::sort(map.data(),
              map.data() + map.size(),
              [&](int i, int j){ return v(i) < v(j); });
}


WaveFunction::~WaveFunction() {}
