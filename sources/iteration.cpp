#include "iteration.h"

using namespace std;
using namespace Eigen;

Iteration::status Iteration::one_step()
{
	// R matrices preparation

	if (!(evalI || evalC))
		return finished;

	auto kRk = Rints.contract(vecC, vecC, 0, 3);
	auto kkR = Rints.contract(vecC, vecC, 0, 1);
	auto ppR = Rints.contract(vecI, vecI, 0, 1);
	auto pRp = Rints.contract(vecI, vecI, 0, 3);

	cout << " R matrices preparation done.\n\n";

	// H matrix elements

	auto Hpp = real(scalar_prod(vecI, H, vecI));
	auto Hkk = real(scalar_prod(vecC, H, vecC));
	auto Hpk = scalar_prod(vecI, H, vecC);

	// Print H matrix elements

	cout << " Energy of the system is: " << En << "\n";
	cout << "\n"
		 << " H matrix elements: "
		 << "\n";
	cout << " H_pp = " << Hpp << "\n"
		 << " H_kk = " << Hkk << "\n"
		 << " H_pk = " << Hpk << "\n\n";

	auto normC = real(scalar_prod(vecC, S, vecC));
	auto normI = real(scalar_prod(vecI, S, vecI));

	cout << "\n"
		 << " Norm squared of the continuum state is: " << normC << "\n";
	cout << "\n"
		 << " Norm squared of the bound state is: " << normI << "\n";

	/// Compose the set of equation of form Ax=ESx

	const int bl = read->get_basis_length();
	const int bnkl = read->get_basis_length_nk();
	//A matrices prep

	MatrixXcd AmatI =
		(normC * Hnk) + (H.topRows(bnkl) * vecC) * (vecC.adjoint() * S.leftCols(bnkl)) +
		(S.topRows(bnkl) * vecC) * (vecC.adjoint() * H.leftCols(bnkl)) +
		Hkk * Snk + kRk + kkR;

	MatrixXcd AmatC =
		U.adjoint() *
		(H + (H * vecI) * (vecI.adjoint() * S) + (S * vecI) * (vecI.adjoint() * H) + Hpp * S + pRp + ppR) *
		U;

	if (job->get_force_orth())
		assert(AmatC.rows() == bnkl);
	else
		assert(AmatC.rows() == (bnkl + 1));

	cout << " A matrices preparation done.\n\n";

	//S matrices prep

	MatrixXcd SmatI =
		normI * Snk + (S.topRows(bnkl) * vecC) * (vecC.adjoint() * S.leftCols(bnkl));

	MatrixXcd SmatC =
		U.adjoint() *
		(S + (S * vecI) * (vecI.adjoint() * S)) *
		U;

	cout << " S matrices preparation done. \n\n";
	cout << " Solving eigenvalue problem. \n\n";

	double sel_paramC, seal_paramI;

	auto method = job->get_selection_mth();

	VectorXcd vecIr_dum, vecCr_dum;

	if (evalC)
	{
		es.compute(AmatC, SmatC);
		VectorXd temp;
		int itC;
		int sizeC = AmatC.cols();

		switch (method)
		{
		case Job_control::selection_mth_t::energy:
			temp = es.eigenvalues() - En * VectorXd::Ones(es.eigenvalues().size());
			temp = temp.cwiseAbs2();
			temp.minCoeff(&itC);
			cout << " Iterator to the matching sel_param of C: " << itC
				 << " \n and it's enegry: " << es.eigenvalues()[itC] << "\n\n";

			break;

		case Job_control::selection_mth_t::norm:
			temp.resize(sizeC);
			for (int i = 0; i < sizeC; ++i)
			{
				vecCr_dum = es.eigenvectors().col(i) / es.eigenvectors()(0, i);
				temp(i) = real(
					vecCr_dum.tail(sizeC - 1).dot(Str.bottomRightCorner(sizeC - 1, sizeC - 1) * vecCr_dum.tail(sizeC - 1)));
			}
			temp.minCoeff(&itC);
			cout << " Iterator to the matching norm of 1eq: " << itC;
			cout << "\n and it's energy: " << es.eigenvalues()[itC] << "\n\n";
			break;
		}
		vecCr_dum = es.eigenvectors().col(itC) / es.eigenvectors()(0, itC);
	}

	if (evalI)
	{
		es.compute(AmatI, SmatI);
		VectorXd temp;
		int itI;

		switch (method)
		{
		case Job_control::selection_mth_t::energy:
			temp = es.eigenvalues() - En * VectorXd::Ones(es.eigenvalues().size());
			temp = temp.cwiseAbs2();
			temp.minCoeff(&itI);
			cout << " Iterator to the matching Energy of 2eq: " << itI
				 << " \n and it's enegry: " << es.eigenvalues()[itI] << "\n\n";

			break;

		case Job_control::selection_mth_t::norm:
			temp.resize(bnkl);
			for (int i = 0; i < bnkl; ++i)
			{
				vecIr_dum = es.eigenvectors().col(i);
				vecIr_dum /= sqrt(real(vecIr_dum.dot(Snk * vecIr_dum)));
				temp(i) = norm(vecIrs.dot(Snk * vecIr_dum));
			}
			temp.maxCoeff(&itI);
			cout << " Iterator to the matching norm of 1eq: " << itI;
			cout << "\n and it's energy: " << es.eigenvalues()[itI] << "\n\n";

			break;
		}

		vecIr_dum = es.eigenvectors().col(itI);
		vecIr_dum /= sqrt(real(vecIr_dum.dot(Snk * vecIr_dum)));
	}

	//Check for self consistency

	double conv_paramC = 0.0;
	double conv_paramI = 0.0;

	if (evalC)
	{
		conv_paramC = (vecCr_dum.cwiseAbs() - vecCr.cwiseAbs()).norm();
		cout << "\n Norm of difference for eq1: " << conv_paramC << "\n";
		vecCr = vecCr_dum;
	}
	if (evalI)
	{
		conv_paramI = (vecIr_dum.cwiseAbs() - vecIr.cwiseAbs()).norm();
		cout << "\n Norm of difference for eq2: " << conv_paramI << "\n";
		vecIr = vecIr_dum;
	}

	vecC = U * vecCr;
	vecI << vecIr, VectorXcd::Zero(bl - bnkl);

	if (!(evalI && evalC))
		return finished;

	if (conv_paramC < treshold && conv_paramI < treshold)
		return self_consistent;
	else
		return running;
}

void Iteration::initialize(const Job_control &controler, const Disk_reader &reader)
{
	assert(reader.is_ready());

	job = std::make_shared<Job_control>(controler);
	read = std::make_shared<Disk_reader>(reader);

	evalI = job->get_computeI();
	evalC = job->get_computeC();

	info = ready;
}

void Iteration::set_energy(const double &energy)
{
	En = energy;
	energy_loaded = true;
}

void Iteration::set_starting_vecs(const Eigen::VectorXcd &vec_ion,
								  const Eigen::VectorXcd &vec_cont)
{
	assert(read->is_ready());

	int bnkl = read->get_basis_length_nk();
	int bl = read->get_basis_length();
	int bkl = read->get_basis_length_k();

	vecI = VectorXcd::Zero(bl);
	vecC = VectorXcd::Zero(bl);

	vecI.head(bnkl) = vec_ion;
	vecC.tail(bkl) = vec_cont;
	starting_vec_loaded = true;
}

void Iteration::load_matrices()
{
	assert(read->is_ready());
	assert(starting_vec_loaded);

	if (evalC || evalI)
		Rints = read->load_Rints();

	H = read->load_H();
	S = read->load_S();

	int bnkl = read->get_basis_length_nk();
	int bl = read->get_basis_length();
	int bkl = read->get_basis_length_k();

	Hnk = H.topLeftCorner(bnkl, bnkl).real();
	Snk = S.topLeftCorner(bnkl, bnkl).real();

	MatrixXd HF = read->load_HFv().real();

	if (job->get_force_orth())
	{
		vecCrs.resize(bnkl);
		vecCrs << 1.0, VectorXcd::Zero(bnkl - 1);
		U = MatrixXcd::Zero(bl, bnkl);
		U.col(0) << -(HF.col(0)).dot(S.topRightCorner(bnkl, bkl) * vecC.tail(bkl)) * HF.col(0),
			vecC.tail(bkl);
		U.block(0, 1, bnkl, bnkl - 1) << HF.rightCols(bnkl - 1);
	}
	else
	{
		vecCrs.resize(bnkl + 1);
		vecCrs << 1.0, VectorXcd::Zero(bnkl);
		U = MatrixXcd::Zero(bl, bnkl + 1);
		U.col(0) << VectorXcd::Zero(bnkl), vecC.tail(bkl);
		U.block(0, 1, bnkl, bnkl) << MatrixXd::Identity(bnkl, bnkl);
	}

	Str = U.adjoint() * S * U;

	vecCr = vecCrs;
	vecC = U * vecCr;

	vecIrs = vecI.head(bnkl);
	vecIr = vecIrs;

	info = ready;
	matrices_loaded = true;
}

void Iteration::iterate()
{
	assert(starting_vec_loaded);
	assert(matrices_loaded);
	assert(energy_loaded);
	Iteration_controler::iterate();
}

void Iteration::free_ints()
{
	Rints.resize(0);
	H.resize(0, 0);
	S.resize(0, 0);
	Hnk.resize(0, 0);
	Snk.resize(0, 0);
	Str.resize(0, 0);
	U.resize(0, 0);
}