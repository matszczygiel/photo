#ifndef ITERATION_H
#define ITERATION_H

#include <complex>
#include <memory>
#include <cassert>

#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>

#include "iteration_controler.h"
#include "two_electron_integrals.h"
#include "job_control.h"
#include "functions.h"
#include "calculator.h"

class Iteration : public Iteration_controler
{
  public:
	Iteration() = default;
	status one_step();

	void initialize(const Job_control &controler, const Disk_reader &reader);
	void load_matrices();
	void set_energy(const double &energy);
	void set_starting_vecs(const Eigen::VectorXcd &vec_ion,
						   const Eigen::VectorXcd &vec_cont);

	inline Eigen::VectorXcd get_vecI() const { return vecI; }
	inline Eigen::VectorXcd get_vecC() const { return vecC; }

	void iterate();
	void free_ints();

  private:
	std::shared_ptr<Job_control> job = nullptr;
	std::shared_ptr<Disk_reader> read = nullptr;

	bool evalI, evalC;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es;

	bool starting_vec_loaded = false;
	Eigen::VectorXcd vecI, vecC;
	Eigen::VectorXcd vecIr, vecCr;
	Eigen::VectorXcd vecIrs, vecCrs;

	bool matrices_loaded = false;
	Tensor_2Ecd Rints;
	Eigen::MatrixXcd H, S;
	Eigen::MatrixXd Hnk, Snk;
	Eigen::MatrixXcd Str;
	Eigen::MatrixXcd U;

	bool energy_loaded = false;
	double En;
};

#endif // ITERATION_H
