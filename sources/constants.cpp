#include "constants.h"

using namespace std;

const vector<int> Const_arrays::belt([]() {
	vector<int> belt(belt_s);
	belt[0] = 1;
	for(unsigned i=1; i< belt.size(); i++)
	belt[i] = -belt[i-1];
	return belt;
}());

const vector<double> Const_arrays::fact([]() {
	vector<double> fact(fact_s);
	fact[0] = 1.0;
	fact[1] = 1.0;
	for(int i=2; i< fact.size(); i++)
	fact[i] = fact[i-1] * i;
	return fact;
}());

const vector<double> Const_arrays::dfact([]() {
	vector<double> dfact(dfact_s);
	dfact[0] = 1.0;
	dfact[1] = 1.0;
	for(unsigned i = 1; i < dfact.size(); i++)
	dfact[i+1] = dfact[i] * ( 2*i + 1 );
	return dfact;
}());

const vector<vector<double>> Const_arrays::binom([]() {
	vector<vector<double>> binom(binom_s);
	for(unsigned i = 0; i < binom.size(); i++) {
		binom[i].resize(binom.size());
		for(unsigned j = 0; j <= i; j++) {
			binom[i][j] = Const_arrays::fact[i];
			binom[i][j] /= Const_arrays::fact[j] * Const_arrays::fact[i-j];
			binom[j][i] = binom[i][j];
		}
	}
	return binom;
}());

const vector<vector<double>> Const_arrays::omega([]() {
	vector<vector<double>> omega(omega_s);
	for(unsigned i = 0; i < omega.size(); i++) {
		omega[i].resize(omega.size());
		for(unsigned j = 0; j <= i; j++)
		omega[i][j] = std::sqrt( ( 2.0*i + 1 ) * Const_arrays::fact[i-j] / ( 2.0 * Const_arrays::fact[i+j] ) );
	}
	return omega;
}());

const vector<unsigned> Const_arrays::crt_siz([]() {
	vector<unsigned> crt_siz(crt_siz_s);
	for(unsigned i = 0; i < crt_siz.size(); i++)
	crt_siz[i] = (i*i + 3*i +2) / 2;
	return crt_siz;
}());

