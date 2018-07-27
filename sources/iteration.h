#ifndef ITERATION_H
#define ITERATION_H

#include "iteration_controler.h"

#include "two_electron_integrals.h"
#include "job_control.h"

#include <complex>
#include <memory>

#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

class Iteration: public Iteration_controler {
public:
	Iteration(Eigen::Tensor<std::complex<double>, 4>& two_el_) :
			Iteration_controler(), two_el(two_el_) {
	}

	status one_step();
	void initialize(const JobControl& jc);
	Eigen::MatrixXcd get_ion() {
		return Ion;
	}
	Eigen::MatrixXcd get_cont() {
		return Cont;
	}

	Eigen::MatrixXcd U_matrix_1eq, H_matrix, S_matrix;
	Eigen::VectorXcd Corr_coef, Coeff_ion;

private:
	unsigned b_l, b_nk_l;
	bool force_ort, eval_1, eval_2, select_m;
	double En, H_pp, H_kk, norm_k_sqrt, norm_b_sqrt;
	std::complex<double> H_pk;
	Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 4>> two_el;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es;
	Eigen::MatrixXcd k_R_k, kk_R, pp_R, p_R_p, H_non_k, S_non_k, A_matrix_2eq, A_matrix_1eq_dum, A_matrix_1eq, S_matrix_2eq, S_matrix_1eq_dum,
			S_matrix_1eq, S_transformed, U_matrix_1eq_H;
	Eigen::VectorXcd Coeff_ion_beginning, Ion, Cont, Corr_coef_dum, Coeff_ion_dum;
};

#endif // ITERATION_H
