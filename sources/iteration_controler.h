#ifndef ITERATION_CONTROLER_H
#define ITERATION_CONTROLER_H

#include <iostream>
#include <cassert>

class Iteration_controler
{
public:
	enum status
	{
		not_started,
		ready,
		running,
		self_consistent,
		finished,
		iterations_limit
	};

	Iteration_controler() = default;

	virtual status one_step() = 0;
	void iterate();

	inline status get_info() const { return info; }
	inline int get_iter_count() const { return iter_count; }
	inline int get_max_iter_count() const { return max_iter_count; }
	inline int get_self_sc_cout() const { return self_sc_cout; }
	inline double get_treshold() const { return treshold; }
	inline int get_max_sc_count() const { return max_sc_count; }

	inline void set_max_iter_count(const int &max_iter_count = 50) { this->max_iter_count = max_iter_count; }
	inline void set_max_sc_count(const int &max_sc_count = 5) { this->max_sc_count = max_sc_count; }
	inline void set_ready() { info = ready; }
	inline void set_treshold(const double &treshold = 0.00001) { this->treshold = treshold; }

protected:
	double treshold = 0.00001;
	status info = not_started;

private:
	int self_sc_cout = 0;
	int max_sc_count = 5;
	int iter_count = 0;
	int max_iter_count = 50;
};

#endif // ITERATION_CONTROLER_H
