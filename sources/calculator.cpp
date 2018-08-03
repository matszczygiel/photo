#include "calculator.h"
#include "harmonics.h"
#include "gamess.h"

void Calculator::initialize(const Job_control &controler, const Disk_reader &reader)
{
    assert(reader.is_ready());
    job = std::make_shared<Job_control>(controler);
    read = std::make_shared<Disk_reader>(reader);
    ready = true;
}

double Calculator::energy() const
{
    assert(ready);
    auto enI = job->get_energy_f();
    auto k = read->get_kval();
    auto iop = job->get_potential();

    return enI + k * k * 0.5 + iop;
}

Eigen::VectorXcd Calculator::continuum_vec() const
{
    assert(ready);

    auto norms = read->load_norms();
    auto lmax = read->get_lmax();
    auto kvec = read->get_kvec();

    std::vector<std::vector<std::vector<double>>> Dfact(lmax + 1);

    for (int l = 0; l <= lmax; l++)
        Dfact[l].resize(l + 1);
    for (int l = 0; l <= lmax; l++)
        for (int p = 0; p <= l; p++)
            Dfact[l][p].resize(l - p + 1, 0.0);

    using std::pow;
    using std::sqrt;

    for (int l = 0; l <= lmax; l++)
        for (int p = 0; p <= l; p++)
            for (int q = 0; q <= l - p; q++)
            {
                for (int m = -l; m <= l; m++)
                    Dfact[l][p][q] += Harmonics::NoNormCalcClmR(l, m, p, q, l - p - q) * Harmonics::Y(l, m, kvec);
                Dfact[l][p][q] *= pow(M_PI, 0.25) * sqrt(Const_arrays::dfact[p] * Const_arrays::dfact[q] * Const_arrays::dfact[l - p - q]) / pow(2.0, 0.25 + l);
            };

    std::vector<std::complex<double>> vec_cont =
        Gamess::order_set<std::complex<double>, double>(Dfact);

    Eigen::VectorXcd vec = Eigen::Map<Eigen::VectorXcd>(vec_cont.data(), vec_cont.size());
    for (int i = 0; i < vec.size(); i++)
        vec(i) /= norms.tail(vec.size())(i);

    return vec;
}

Eigen::VectorXd Calculator::bound_vec() const
{
    assert(ready);

    auto H = read->load_H();
    auto S = read->load_S();
    auto bnkl = read->get_basis_length_nk();

    Eigen::MatrixXd Hnk = H.topLeftCorner(bnkl, bnkl).real();
    Eigen::MatrixXd Snk = S.topLeftCorner(bnkl, bnkl).real();

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(Hnk, Snk);

    return es.eigenvectors().col(0);
}