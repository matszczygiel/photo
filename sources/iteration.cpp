#include "iteration.h"

#include "fun.h"


using namespace std;
using namespace Eigen;

Iteration::status Iteration::one_step() {
    // R matrices preparation

    if(!(eval_1 || eval_2)) return finished;

    k_R_k.setZero();

    TensorMap<Tensor<complex<double>, 2>> k_R_k_tens(k_R_k.data(), b_nk_l, b_nk_l);
    TensorMap<Tensor<complex<double>, 2>> kk_R_tens(kk_R.data(), b_nk_l, b_nk_l);
    TensorMap<Tensor<complex<double>, 2>> pp_R_tens(pp_R.data(), b_l, b_l);
    TensorMap<Tensor<complex<double>, 2>> p_R_p_tens(p_R_p.data(), b_l, b_l);



    for(unsigned i=0; i<b_nk_l; i++) for(unsigned j=0; j<b_nk_l; j++)
        for(unsigned a =0; a<b_l; a++) for(unsigned b =0; b<b_l; b++)
            k_R_k(i,j) += conj(Cont(a)) * two_el( a, j, i, b ) * Cont(b);

    kk_R.setZero();

    for(unsigned i=0; i<b_nk_l; i++) for(unsigned j=0; j<b_nk_l; j++)
        for(unsigned a =0; a<b_l; a++) for(unsigned b =0; b<b_l; b++)
            kk_R(i,j) += conj(Cont(a)) * Cont(b) * two_el( a ,b ,i ,j );

    pp_R.setZero();

    for(unsigned n =0; n< b_l; n++) for(unsigned m = 0; m< b_l; m++)
        for(unsigned i =0; i< b_nk_l; i++) for(unsigned j =0; j<b_nk_l; j++)
            pp_R(n,m) += conj(Ion(i)) * Ion(j) * two_el( i, j, n, m );

    p_R_p.setZero();

    for(unsigned n =0; n< b_l; n++) for(unsigned m = 0; m< b_l; m++)
        for(unsigned i =0; i< b_nk_l; i++) for(unsigned j =0; j<b_nk_l; j++)
            p_R_p(i,j) += conj(Ion(i)) * Ion(j) * two_el( i, m, n, j );

    cout << " R matrices preparation done.\n\n";

    // H matrix elements

    H_pp = real( Coeff_ion.dot( H_non_k * Coeff_ion) );
    H_kk = real( Cont.dot(H_matrix * Cont));
    H_pk = Ion.dot( H_matrix * Cont );

    // Print H matrix elements

    cout << " Energy of the system is: " << En << "\n";
    cout << "\n" << " H matrix elements: " << "\n";
    cout << " H_pp = " << H_pp << "\n" << " H_kk = " << H_kk
         << "\n" << " H_pk = " << H_pk << "\n\n";

    // Norm of |k>

    norm_k_sqrt = real( Cont.dot(S_matrix * Cont));
    cout << "\n" << " Norm squared of the continuum state is: " << norm_k_sqrt << "\n";

    // norm of |b> state

    norm_b_sqrt = real( Coeff_ion.dot(S_non_k * Coeff_ion) );

    cout << "\n" << " Norm squared of the bound state is: " << norm_b_sqrt << "\n";

    /// Compose the set of equation of form Ax=ESx

    //A matrices prep

    A_matrix_2eq = ( norm_k_sqrt * H_non_k ) +
            ( H_matrix.topRows(b_nk_l) * Cont ) * ( Cont.adjoint() * S_matrix.leftCols(b_nk_l) ) +
            ( S_matrix.topRows(b_nk_l) * Cont ) * ( Cont.adjoint() * H_matrix.leftCols(b_nk_l) ) +
            ( H_kk * S_non_k ) +  k_R_k + kk_R;

    A_matrix_1eq_dum = H_matrix +
            ( H_matrix * Ion ) * ( Ion.adjoint() * S_matrix ) +
            ( S_matrix * Ion ) * ( Ion.adjoint() * H_matrix ) +
            H_pp * S_matrix + p_R_p + pp_R;

    A_matrix_1eq =  U_matrix_1eq_H * A_matrix_1eq_dum * U_matrix_1eq;

    if ( force_ort == true )
        assert( A_matrix_1eq.rows() == b_nk_l );
    else
        assert( A_matrix_1eq.rows() == (b_nk_l + 1) );

    cout << " A matrices preparation done.\n\n";

    //S matrices prep

    S_matrix_2eq = norm_k_sqrt * S_non_k  +
            (S_matrix.topRows(b_nk_l) * Cont)*
            (Cont.adjoint()* S_matrix.leftCols(b_nk_l));

    S_matrix_1eq_dum = S_matrix + (S_matrix * Ion) * (Ion.adjoint() * S_matrix);

    S_matrix_1eq = U_matrix_1eq_H * S_matrix_1eq_dum * U_matrix_1eq;

    cout << " S matrices preparation done. \n\n";
    cout << "\n Solving eigenvalue problem. \n\n";

    int low_it_1eq, low_it_2eq;
    double Energy_1eq, Energy_2eq;

    int size_1eq = A_matrix_1eq.cols();

    if(eval_1) {
        es.compute(A_matrix_1eq, S_matrix_1eq);
        cout <<" 1 eq (continuum) energies: \n\n" << es.eigenvalues() << "\n\n";

        if (select_m == 0) {
            low_it_1eq =0;
            Energy_1eq = es.eigenvalues()[0];
            for(int i=0; i<size_1eq; i++)
                if( abs(Energy_1eq- En) > abs(es.eigenvalues()[i] - En) ) {
                    low_it_1eq = i;
                    Energy_1eq = es.eigenvalues()[i];
                }
            cout << " Iterator to the matching Energy of 1eq: " << low_it_1eq << " \n and it's enegry: " << Energy_1eq << "\n\n";

            Corr_coef_dum = es.eigenvectors().col(low_it_1eq) / es.eigenvectors()(0, low_it_1eq);
        }
        else if (select_m == 1) {
            low_it_1eq = 0;

            Corr_coef_dum = es.eigenvectors().col(0) / es.eigenvectors()(0, 0);

            Energy_1eq = real(
                        (Corr_coef_dum.tail(size_1eq-1)).dot(
                            S_transformed.bottomRightCorner(size_1eq-1, size_1eq-1) * Corr_coef_dum.tail(size_1eq-1) ));
            double norm_temp;

            for(int i = 1; i<size_1eq; i++) {
                Corr_coef_dum = es.eigenvectors().col(i) / es.eigenvectors()(0, i);
                norm_temp = real(
                            (Corr_coef_dum.tail(size_1eq-1)).dot(
                                S_transformed.bottomRightCorner(size_1eq-1, size_1eq-1) * Corr_coef_dum.tail(size_1eq-1) ));

                if( abs(Energy_1eq) > abs(norm_temp) ) {
                    low_it_1eq = i;
                    Energy_1eq = norm_temp;
                }
            }
            cout << " Iterator to the matching norm of 1eq: " << low_it_1eq << " \n and it's norm: " << Energy_1eq;
            cout<< "\n and it's energy: " << es.eigenvalues()[low_it_1eq] << "\n\n";

            Corr_coef_dum = es.eigenvectors().col(low_it_1eq) / es.eigenvectors()(0, low_it_1eq);
        }
    }

    if(eval_2) {
        es.compute(A_matrix_2eq, S_matrix_2eq);
        cout <<" 2 eq (bounded) energies: \n\n" << es.eigenvalues() << "\n\n";

        if (select_m == 0)
        {

            low_it_2eq = 0;

            Energy_2eq = es.eigenvalues()[0];
            for(int i=0; i<b_nk_l; i++)
                if( abs(Energy_2eq- En) > abs(es.eigenvalues()[i] - En) ) {
                    low_it_2eq = i;
                    Energy_2eq = es.eigenvalues()[i];

                }
            cout << " Iterator to the matching Energy of 2eq: " << low_it_2eq << " \n and it's enegry: " << Energy_2eq << "\n\n";

            Coeff_ion_dum = es.eigenvectors().col(low_it_2eq);
            Coeff_ion_dum /= sqrt( real( Coeff_ion_dum.dot( S_non_k * Coeff_ion_dum)) );
        }
        else if (select_m == 1) {
            low_it_2eq = 0;

            Coeff_ion_dum = es.eigenvectors().col(0);
            Coeff_ion_dum /= sqrt( real( Coeff_ion_dum.dot( S_non_k * Coeff_ion_dum) ) );
            Energy_2eq = norm( Coeff_ion_beginning.dot( S_non_k * Coeff_ion_dum ));
            double norm_temp;

            for(int i = 1; i < b_nk_l; i++)
            {
                Coeff_ion_dum = es.eigenvectors().col(i);
                Coeff_ion_dum /= sqrt( real( Coeff_ion_dum.dot( S_non_k * Coeff_ion_dum) ) );
                norm_temp = norm( Coeff_ion_beginning.dot( S_non_k * Coeff_ion_dum ));

                if( norm_temp > Energy_2eq ) {
                    low_it_2eq = i;
                    Energy_2eq = norm_temp;
                }
            }
            cout << " Iterator to the matching solution of 2eq: " << low_it_2eq << " \n and it's coefficient: " << Energy_2eq;
            cout<< "\n and it's energy: " << es.eigenvalues()[low_it_2eq] << "\n\n";

            Coeff_ion_dum = es.eigenvectors().col(low_it_2eq);
            Coeff_ion_dum /= sqrt( real( Coeff_ion_dum.dot( S_non_k * Coeff_ion_dum)) );
        }
    }

    if ( !(eval_1 && eval_2) )
        return finished;

    //Check for self consistency

    double diff_norm_eq1 = 0.0;
    double diff_norm_eq2 = 0.0;

    if(eval_1) {
        diff_norm_eq1 = (Corr_coef_dum.cwiseAbs() - Corr_coef.cwiseAbs()).norm();
        cout << "\n Norm of difference for eq1: " << diff_norm_eq1 << "\n";
        Corr_coef = Corr_coef_dum;
    }
    if(eval_2) {
        diff_norm_eq2 = (Coeff_ion_dum.cwiseAbs() - Coeff_ion.cwiseAbs()).norm();
        cout << "\n Norm of difference for eq2: " << diff_norm_eq2 << "\n";
        Coeff_ion = Coeff_ion_dum;
    }

    Cont = U_matrix_1eq * Corr_coef;
    Ion << Coeff_ion, VectorXcd::Zero(b_l - b_nk_l);

    if(diff_norm_eq1 < treshold && diff_norm_eq2 < treshold)
        return self_consistent;
    else
        return running;
}



void Iteration::initialize(const JobControl &jc) {
    b_l = jc.getBL();
    b_nk_l = jc.getBNKL();

    assert(two_el.dimension(0) == b_l);
    assert(H_matrix.rows() == b_l);
    assert(S_matrix.rows() == b_l);

    force_ort = jc.getForceOrthogonality();
    eval_1 = jc.getEval1eq();
    eval_2 = jc.getEval2eq();
    select_m = jc.getSelectionMethod();
    En = Energy(jc.getK(), jc.getIB());

    H_non_k = H_matrix.topLeftCorner(b_nk_l, b_nk_l);
    S_non_k = S_matrix.topLeftCorner(b_nk_l, b_nk_l);

    k_R_k.resize( b_nk_l, b_nk_l);
    kk_R.resize( b_nk_l, b_nk_l);
    p_R_p.resize( b_l, b_l);
    pp_R.resize( b_l, b_l);

    Cont = U_matrix_1eq * Corr_coef;
    Ion.resize(b_l);
    Ion << Coeff_ion, VectorXcd::Zero(b_l - b_nk_l);
    Coeff_ion_beginning = Coeff_ion;

    U_matrix_1eq_H = U_matrix_1eq.adjoint();
    S_transformed = U_matrix_1eq_H * S_matrix * U_matrix_1eq;

    info = ready;
}
























