#include "iteration_controler.h"

void Iteration_controler::iterate() {
	assert(info == ready);

	std::cout << "\n" << " Starting iteration.\n";

	info = running;
	while (info == running) {
		iter_count++;
		std::cout << "\n";
		std::cout << "========= Iteration " << iter_count << " =============" << "\n\n";

		info = one_step();

		if (info == self_consistent) {
			self_sc_cout++;
			info = running;
		} else
			self_sc_cout = 0;

		if (self_sc_cout == max_sc_count) {
			info = finished;
			break;
		}
		if (iter_count == max_iter_count) {
			info = iterations_limit;
			std::cout << " Iterations limit reached. \n";
			break;
		}
		if (self_sc_cout != 0)
			std::cout << "\n" << " Self consistency counter:" << self_sc_cout << "\n\n";
	}
	std::cout << "\n";
	std::cout << " ======= End of iteration " << " ===========" << "\n\n";
}

