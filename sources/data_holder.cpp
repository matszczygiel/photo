#include "data_holder.h"

void Data_holder::free()
{
    H.resize(0, 0);
    S.resize(0, 0);
    Gaugex.resize(0, 0);
    Gaugey.resize(0, 0);
    Gaugez.resize(0, 0);
    HF.resize(0, 0);

    Rints.resize(0);
}

void Data_holder::load(const Disk_reader &r)
{
    assert(r.is_ready());

    H = r.load_H;
    S = r.load_S;
    Gaugex = r.load_Gaugex;
    Gaugey = r.load_Gaugey;
    Gaugez = r.load_Gaugez;
    HF;

    Rints = r.load_Rints;
}
