#ifndef HARMONICS_H
#define HARMONICS_H

#include "constants.h"

#include <cmath>
#include <algorithm>

#include <eigen3/Eigen/Dense>

class Harmonics
{
  public:
	/**
	 * Returns the coefficents of cartesian monomial in r^l * Y_{lm}, where Y_{lm} is REAL spherical harmonic. 
	 */
	static double NoNormCalcClmR(const int &l, const int &m, const int &lx, const int &ly, const int &lz);

	/**
	 * Return value of r^l * Y_{lm} (r), where Y_{lm} is REAL spherical harmonic. 
	 */
	static double NoNormYrl(const int &l, const int &m, const Eigen::Vector3d &r);

	/**
	 * Return value of Y_{lm} (r), where Y_{lm} is REAL spherical harmonic. 
	 */
	static double NoNormY(const int &l, const int &m, const Eigen::Vector3d &r);
};

#endif // HARMONICS_H
