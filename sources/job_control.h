#ifndef JOB_CONTROL_H
#define JOB_CONTROL_H

#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <iomanip>

#include "input_data.h"

class Job_control {
public:
	enum selection_m_type {
		energy, norm
	};

	enum gauge_type {
		velocity, dipole
	};

	Job_control(const Input_data& in);

private:
	bool compute_bound;
	bool compute_cont;
	bool compute_ion_state;

	double k_val;
	double i_energy;
	double f_energy;
	double potential;

	std::string data_path;
	std::string basis_file;
	std::string E1_file;
	std::string E2_file;
	std::string norm_file;
	std::string HFv_file;
	std::string CI_file;
	std::string res_file;

	selection_m_type selection_m;
	gauge_type gauge;
	bool force_orth;

	bool write;
	std::string res_file;
};

#endif // JOB_CONTROL_H
