#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>
#include <cassert>
#include <vector>
#include <cmath>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include "job_control.h"
#include "disk_reader.h"

class Calculator
{
public:
void initialize(const Job_control& controler, const Disk_reader& reader);

double energy() const;

Eigen::VectorXcd continuum_vec() const;
Eigen::VectorXd bound_vec() const;

private:
    std::shared_ptr<Job_control> job = nullptr;
    std::shared_ptr<Disk_reader> read = nullptr;
    bool ready = false;
};


#endif