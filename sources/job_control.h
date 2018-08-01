#ifndef JOB_CONTROL_H
#define JOB_CONTROL_H

#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <iomanip>

#include "input_data.h"

class Job_control
{
  public:
	enum selection_mth_t
	{
		energy,
		norm
	};

	enum gauge_t
	{
		velocity,
		dipole
	};

	Job_control();

	void read(Input_data &in);

	inline std::string get_file_basis() const { return path_in + file_basis; }
	inline std::string get_file_1E() const { return path_in + file_1E; }
	inline std::string get_file_2E() const { return path_in + file_2E; }
	inline std::string get_file_norm() const { return path_in + file_norm; }
	inline std::string get_file_HFv() const { return path_in + file_HFv; }
	inline std::string get_file_CI() const { return path_in + file_CI; }
	inline std::string get_continuum_id() const { return continuum_id; }

	inline auto get_gauge() const { return gauge; }

  private:
	bool compute_bound;
	bool compute_cont;
	bool compute_ion_state;

	double k_val;
	double energy_i;
	double energy_f;
	double potential;

	std::string path_in;
	std::string file_basis;
	std::string file_1E;
	std::string file_2E;
	std::string file_norm;
	std::string file_HFv;
	std::string file_CI;

	std::string continuum_id;

	selection_mth_t selection_m;
	gauge_t gauge;
	bool force_orth;

	bool write;
	std::string path_out;
	std::string file_out;
};

#endif // JOB_CONTROL_H
