#ifndef ITERATION_H
#define ITERATION_H

#include <complex>
#include <memory>
#include <cassert>

#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

#include "iteration_controler.h"
#include "two_electron_integrals.h"
#include "job_control.h"
#include "data_holder.h"
#include "fun.h"

class Iteration : public Iteration_controler
{
  public:
	Iteration() = default;

	status one_step();
	void initialize(const Job_control &controler, const Disk_reader &reader);
	void load_ints();

	Eigen::MatrixXcd get_ion()
	{
		return Ion;
	}
	Eigen::MatrixXcd get_cont()
	{
		return Cont;
	}

	Eigen::MatrixXcd U_matrix_1eq, H_matrix, S_matrix;
	Eigen::VectorXcd Corr_coef, Coeff_ion;

  private:
	std::shared_ptr<Job_control> job = nullptr;
	//std::shared_ptr<Data_holder> hold = nullptr;
	std::shared_ptr<Disk_reader> read = nullptr;

	bool evalI, evalC;
	double En, H_pp, H_kk, norm_k_sqrt, norm_b_sqrt;
	std::complex<double> H_pk;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es;

	Eigen::MatrixXcd k_R_k, kk_R, pp_R, p_R_p, H_non_k, S_non_k, A_matrix_2eq, A_matrix_1eq_dum, A_matrix_1eq, S_matrix_2eq, S_matrix_1eq_dum,
		S_matrix_1eq, S_transformed, U_matrix_1eq_H;
	Eigen::VectorXcd Coeff_ion_beginning, Corr_coef_dum, Coeff_ion_dum;

	Eigen::MatrixXcd U;


	Eigen::VectorXcd vecI, vecC;

	bool ints_loaded = false;
	Tensor_2Ecd Rints;
	Eigen::MatrixXcd H, S;
	Eigen::MatrixXcd Hnk, Snk;
};

#endif // ITERATION_H
