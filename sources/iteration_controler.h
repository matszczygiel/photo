#ifndef ITERATION_CONTROLER_H
#define ITERATION_CONTROLER_H

#include <iostream>

class Iteration_controler
{
public:
    enum status {not_started,
                 ready,
                 running,
                 self_consistent,
                 finished,
                 iterations_limit};

    Iteration_controler() = default;

    virtual status one_step() = 0;

    void iterate();
    void set_treshold( double treshold_ ) { treshold = treshold_; }
    void set_max_iterations( int max ) { iter_count = max; }
    void set_max_selfconsistency_count( int max ) { max_sc_count = max; }

    int get_iterations_count() const { return iter_count; }
    status get_info() const { return info; }
    void set_ready() { info = ready; }

protected:
    double  treshold        = 0.00001;
    status  info = not_started;

private:
    int     self_sc_cout    = 0;
    int     max_sc_count    = 5;
    int     iter_count      = 0;
    int     max_iter_count  = 50;
};

#endif // ITERATION_CONTROLER_H
