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

	int low_itC, low_itI;
	double EnergyC, EnergyI;

	int sizeC = AmatC.cols();

	if (evalC)
	{
		es.compute(AmatC, SmatC);
		cout << " 1 eq (continuum) energies: \n\n"
			 << es.eigenvalues() << "\n\n";

		if (select_m == 0)
		{
			low_itC = 0;
			EnergyC = es.eigenvalues()[0];
			for (int i = 0; i < sizeC; i++)
				if (abs(EnergyC - En) > abs(es.eigenvalues()[i] - En))
				{
					low_itC = i;
					EnergyC = es.eigenvalues()[i];
				}
			cout << " Iterator to the matching Energy of 1eq: " << low_itC << " \n and it's enegry: " << EnergyC << "\n\n";

			Corr_coef_dum = es.eigenvectors().col(low_itC) / es.eigenvectors()(0, low_itC);
		}
		else if (select_m == 1)
		{
			low_itI = 0;

			Corr_coef_dum = es.eigenvectors().col(0) / es.eigenvectors()(0, 0);

			EnergyC = real(
				(Corr_coef_dum.tail(sizeC - 1)).dot(S_transformed.bottomRightCorner(sizeC - 1, sizeC - 1) * Corr_coef_dum.tail(sizeC - 1)));
			double norm_temp;

			for (int i = 1; i < sizeC; i++)
			{
				Corr_coef_dum = es.eigenvectors().col(i) / es.eigenvectors()(0, i);
				norm_temp = real(
					(Corr_coef_dum.tail(sizeC - 1)).dot(S_transformed.bottomRightCorner(sizeC - 1, sizeC - 1) * Corr_coef_dum.tail(sizeC - 1)));

				if (abs(EnergyC) > abs(norm_temp))
				{
					low_itC = i;
					EnergyC = norm_temp;
				}
			}
			cout << " Iterator to the matching norm of 1eq: " << low_itC << " \n and it's norm: " << EnergyC;
			cout << "\n and it's energy: " << es.eigenvalues()[low_itC] << "\n\n";

			Corr_coef_dum = es.eigenvectors().col(low_itC) / es.eigenvectors()(0, low_itC);
		}
	}

	if (eval_2)
	{
		es.compute(AmatI, SmatI);
		cout << " 2 eq (bounded) energies: \n\n"
			 << es.eigenvalues() << "\n\n";

		if (select_m == 0)
		{

			low_itI = 0;

			EnergyI = es.eigenvalues()[0];
			for (int i = 0; i < b_nk_l; i++)
				if (abs(EnergyI - En) > abs(es.eigenvalues()[i] - En))
				{
					low_itI = i;
					EnergyI = es.eigenvalues()[i];
				}
			cout << " Iterator to the matching Energy of 2eq: " << low_itI << " \n and it's enegry: " << EnergyI << "\n\n";

			Coeff_ion_dum = es.eigenvectors().col(low_itI);
			Coeff_ion_dum /= sqrt(real(Coeff_ion_dum.dot(S_non_k * Coeff_ion_dum)));
		}
		else if (select_m == 1)
		{
			low_itI = 0;

			Coeff_ion_dum = es.eigenvectors().col(0);
			Coeff_ion_dum /= sqrt(real(Coeff_ion_dum.dot(S_non_k * Coeff_ion_dum)));
			EnergyI = norm(Coeff_ion_beginning.dot(S_non_k * Coeff_ion_dum));
			double norm_temp;

			for (int i = 1; i < b_nk_l; i++)
			{
				Coeff_ion_dum = es.eigenvectors().col(i);
				Coeff_ion_dum /= sqrt(real(Coeff_ion_dum.dot(S_non_k * Coeff_ion_dum)));
				norm_temp = norm(Coeff_ion_beginning.dot(S_non_k * Coeff_ion_dum));

				if (norm_temp > EnergyI)
				{
					low_itI = i;
					EnergyI = norm_temp;
				}
			}
			cout << " Iterator to the matching solution of 2eq: " << low_itI << " \n and it's coefficient: " << EnergyI;
			cout << "\n and it's energy: " << es.eigenvalues()[low_itI] << "\n\n";

			Coeff_ion_dum = es.eigenvectors().col(low_itI);
			Coeff_ion_dum /= sqrt(real(Coeff_ion_dum.dot(S_non_k * Coeff_ion_dum)));
		}
	}

	//Check for self consistency

	double diff_norm_eq1 = 0.0;
	double diff_norm_eq2 = 0.0;

	if (eval_1)
	{
		diff_norm_eq1 = (Corr_coef_dum.cwiseAbs() - Corr_coef.cwiseAbs()).norm();
		cout << "\n Norm of difference for eq1: " << diff_norm_eq1 << "\n";
		Corr_coef = Corr_coef_dum;
	}
	if (eval_2)
	{
		diff_norm_eq2 = (Coeff_ion_dum.cwiseAbs() - Coeff_ion.cwiseAbs()).norm();
		cout << "\n Norm of difference for eq2: " << diff_norm_eq2 << "\n";
		Coeff_ion = Coeff_ion_dum;
	}

	cout << "New continuum: \n";
	cout << Corr_coef;
	cout << "\n\n";

	cout << "New bound: \n";
	cout << Coeff_ion;
	cout << "\n\n";

	Cont = U_matrix_1eq * Corr_coef;
	Ion << Coeff_ion, VectorXcd::Zero(b_l - b_nk_l);

	if (!(eval_1 && eval_2))
		return finished;

	if (diff_norm_eq1 < treshold && diff_norm_eq2 < treshold)
		return self_consistent;
	else
		return running;
}

void Iteration::initialize(const Job_control &controler, const Disk_reader &reader)
{
	assert(reader.is_ready());

	job = std::make_shared<Job_control>(controler);
	read = std::make_shared<Disk_reader>(reader);

	eval_1 = jc.getEval1eq();
	eval_2 = jc.getEval2eq();
	select_m = jc.getSelectionMethod();
	En = Energy(jc.getK(), jc.getIB());

	H_non_k = H_matrix.topLeftCorner(b_nk_l, b_nk_l);
	S_non_k = S_matrix.topLeftCorner(b_nk_l, b_nk_l);

	k_R_k.resize(b_nk_l, b_nk_l);
	kk_R.resize(b_nk_l, b_nk_l);
	p_R_p.resize(b_l, b_l);
	pp_R.resize(b_l, b_l);

	Cont = U_matrix_1eq * Corr_coef;
	Ion.resize(b_l);
	Ion << Coeff_ion, VectorXcd::Zero(b_l - b_nk_l);
	Coeff_ion_beginning = Coeff_ion;

	U_matrix_1eq_H = U_matrix_1eq.adjoint();
	S_transformed = U_matrix_1eq_H * S_matrix * U_matrix_1eq;

	info = ready;
}

void Iteration::load_ints()
{
	Rints = read->load_Rints();
	H = read->load_H();
	S = read->load_S();

	const bnkl = read->get_basis_length_nk();
	Hnk = H.topLeftCorner(bnkl, bnkl);
	Snk = S.topLeftCorner(bnkl, bnkl);
}
/*
	const int bnkl = read->get_basis_length_nk();
	const int bkl = read->get_basis_length_k();
	const int bl = read->get_basis_length();

	auto HF = read->load_HFv();
	MatrixXd Utemp = HF.leftCols(bnkl - 1);

	if (job->get_force_orth())
	{
		Corr_coef.resize(bnkl);
		Corr_coef << 1.0, VectorXcd::Zero(b_nk_l - 1);

		U = MatrixXcd::Zero(bl, bnkl);
		U.col(0) << - scalar_prod(HF.col(0), S.topRightCorner(bnkl, bkl), )
		-(HF.col(0)).dot(S_matrix.topRightCorner(b_nk_l, b_k_l) * Cont_coef) * HF_ground, Cont_coef;
		U_.block(0, 1, bnkl, bnkl - 1) << Utemp;
	}
	else
	{
		Corr_coef.resize(b_nk_l + 1);
		Corr_coef << 1.0, VectorXcd::Zero(b_nk_l);

		U_matrix_1eq = MatrixXcd::Zero(b_l, b_nk_l + 1);

		U_matrix_1eq.col(0) << VectorXcd::Zero(b_nk_l), Cont_coef;
		U_matrix_1eq.block(0, 1, b_nk_l, b_nk_l) << MatrixXd::Identity(b_nk_l, b_nk_l);
	}





	vecI = 

}
*/