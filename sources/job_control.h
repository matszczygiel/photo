#ifndef JOB_CONTROL_H
#define JOB_CONTROL_H

#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "fun.h"


class JobControl {

private:
    bool check;
    bool iter;
    bool eval_1eq;
    bool eval_2eq;
    bool forceOrthogonality;
    bool writeResult;
    bool solveIon;

    std::string data_path;
    std::string file1E;
    std::string file2E;
    std::string norms_file;
    std::string HF_file = "c.dat";
    std::string HF_en_file = "e.dat";
    std::string basis_name;
    std::string CI_file = "ci.dat";
    std::string res_file = "res.dat";

    double ionizatoin_pot;
    double k;
    double E_he;
    double E_ion = -2.0;

    std::vector<unsigned> gto_number;

    int ion_pot_int;
    int eq_int;
    int selectionMethod;

    unsigned b_l;
    unsigned b_k_l;
    unsigned b_nk_l;
    unsigned b_l_sqrt;
    unsigned b_k_l_sqrt;
    unsigned b_nk_l_sqrt;
    unsigned l_max;

    constexpr static unsigned crt_siz[11] = { 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 };

public:
    JobControl() {}

    void setInFiles( double kk );

    std::string getFile1E() const       { return data_path + file1E; }
    std::string getFile2E() const       { return data_path + file2E; }
    std::string getFileNorm() const     { return data_path + norms_file; }
    std::string getFileHF() const       { return data_path + HF_file; }
    std::string getFileHFEnergy() const { return data_path + HF_en_file; }
    std::string getFileCI() const       { return data_path + CI_file; }
    unsigned    getBNKL() const         { return b_nk_l; }
    unsigned    getBKL() const          { return b_k_l; }
    unsigned    getBL() const           { return b_l; }
    unsigned    getBNKLsqrt() const     { return b_nk_l_sqrt; }
    unsigned    getBKLsqrt() const      { return b_k_l_sqrt; }
    unsigned    getBLsqrt() const       { return b_l_sqrt; }
    bool        getIter() const         { return iter; }
    bool        getEval1eq() const      { return eval_1eq; }
    bool        getEval2eq() const      { return eval_2eq; }
    unsigned    getLmax() const         { return l_max; }
    double      getK() const            { return k; }
    double      getIB() const           { return ionizatoin_pot; }
    bool getForceOrthogonality() const  { return forceOrthogonality; }
    bool getIfWrite() const             { return writeResult; }
    bool getSolveIon() const            { return solveIon; }
    int getSelectionMethod() const      { return selectionMethod; }


    void readInput( std::string in_file );
    void calcualteBasisFunctNumber();
    void print() const;
    void writeRes(double sig) const;

};

#endif // JOB_CONTROL_H
