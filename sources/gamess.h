#ifndef GAMESS_H
#define GAMESS_H

#include <stdexcept>
#include <vector>
#include <complex>

#include "constants.h"

class Gamess
{
  public:
    static int pos_change_gamess(const int &l, const int &pos);

    template <class type1, class type2>
    static std::vector<type1> order(const std::vector<std::vector<type2>> &shl_crt,
                                    const int &l)
    {
        /* transform to the Gamess shell indexing */
        int pp_gam, pos;
        std::vector<type1> shl_crt_dum(Const_arrays::crt_siz.at(l));
        for (int p = 0; p <= l; p++)
            for (int q = 0; q <= l - p; q++)
            {
                pos = (l + 1) * p - p * (p - 1) / 2 + q;
                p_gam = pos_change_gamess(l, pos);
                shl_crt_dum[p_gam] = shl_crt[p][q];
            }

        return shl_crt_dum;
    }

    template <class type1, class type2>
    static std::vector<type1> order_set(const std::vector<std::vector<std::vector<type2>>> shl_crt, const int &l_max)
    {
        std::vector<type1> shl_crt_dum;

        for (int l = 0; l <= l_max; l++)
        {
            std::vector<type1> mem = order(shl_crt.at(l), l);

            for (int i = 0; i < crt_siz[l]; i++)
                shl_crt_dum.push_back(mem[i]);
        }
        return shl_crt_dum;
    }
};

#endif