#include "iteration_controler.h"

using namespace std;


void Iteration_controler::iterate() {
    if (info != ready) {
        cout << " Iteration not ready!\n";
        return;
    }

    cout << "\n" << " Starting iteration.\n";

    info = running;
    while (info == running) {
        iter_count++;
        cout << "\n";
        cout << "========= Iteration " << iter_count << " =============" << "\n\n";

        info = this->one_step();

        if (info == self_consistent) {
            self_sc_cout++;
            info = running;
        }
        else
            self_sc_cout = 0;

        if (self_sc_cout == max_sc_count) {
            info = finished;
            break;
        }
        if (iter_count == max_iter_count) {
            info = iterations_limit;
            cout << " Iterations limit reached. \n";
            break;
        }
        if (self_sc_cout != 0)
            cout << "\n" << " Self consistency counter:" << self_sc_cout << "\n\n";
    }
    cout << "\n";
    cout << " ======= End of iteration " << " ===========" << "\n\n";
}


