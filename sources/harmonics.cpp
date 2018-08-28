#include "harmonics.h"

double Harmonics::NoNormCalcClmR(const int &l, const int &m, const int &lx, const int &ly, const int &lz) {
    if (lx + ly + lz != l)
        return 0.0;

    if (m > 0) {
        int bs, am, lima, lim1b, lim2b, twol, twoa, twob;
        double val, fac, faci;
        am   = std::abs(m);
        bs   = l - lz - am;
        twol = 2 * l;
        if (bs % 2 == 1 || bs < 0)
            return 0.0;
        bs  = bs / 2;
        fac = M_SQRT2 * Const_arrays::omega[l][m] * Const_arrays::belt[am] / (std::pow(2.0, l) * Const_arrays::fact[l]);

        lima  = l - am;
        lima  = (lima % 2 == 0) ? lima / 2 : (lima - 1) / 2;
        lim1b = lx - am;
        lim1b = (lim1b % 2 == 0) ? lim1b / 2 : (lim1b + 1) / 2;
        lim1b = std::max(lim1b, 0);
        lim2b = (lx % 2 == 0) ? lx / 2 : (lx - 1) / 2;
        lim2b = std::min(lim2b, bs);

        val = 0.0;
        for (int a = std::max(bs, 0); a <= lima; a++) {
            twoa = 2 * a;
            faci = Const_arrays::binom[l][a] * Const_arrays::belt[a] * Const_arrays::fact[twol - twoa] * Const_arrays::binom[a][bs] / Const_arrays::fact[l - am - twoa];
            for (int b = lim1b; b <= lim2b; b++) {
                twob = 2 * b + am - lx;
                if (twob % 2 == 1)
                    continue;
                val += Const_arrays::binom[bs][b] * Const_arrays::binom[am][twob] * Const_arrays::belt[twob / 2] * faci;
            }
        }
        return Const_arrays::belt[am] * val * fac / std::sqrt(2.0 * M_PI);
    }
    if (m == 0) {
        int lxx, lxy, lim1a, lim2a, twoa, twol;
        lxx = lx;
        lxy = lx + ly;
        if (lxx % 2 == 1)
            return 0.0;
        if (lxy % 2 == 1)
            return 0.0;
        twol = 2 * l;
        lxx /= 2;
        lxy /= 2;
        double val, fac;
        fac   = std::sqrt((2 * l + 1) / 2.) * Const_arrays::binom[lxy][lxx] / (std::pow(2.0, l) * Const_arrays::fact[l]);
        lim1a = lxy;
        lim2a = (l % 2 == 0) ? l / 2 : (l - 1) / 2;

        val = 0.0;
        for (int a = lim1a; a <= lim2a; a++) {
            twoa = 2 * a;
            val += Const_arrays::binom[l][a] * Const_arrays::belt[a] * Const_arrays::fact[twol - twoa] * Const_arrays::binom[a][lxy] / Const_arrays::fact[l - twoa];
        }
        return val * fac / std::sqrt(2.0 * M_PI);
    }
    if (m < 0) {
        int bs, am, lima, lim1b, lim2b, twol, twoa, twob;
        double val, fac, faci;
        am   = std::abs(m);
        bs   = l - lz - am;
        twol = 2 * l;
        if (bs % 2 == 1 || bs < 0)
            return 0.0;
        bs  = bs / 2;
        fac = -M_SQRT2 * Const_arrays::omega[l][am] * Const_arrays::belt[am] / (std::pow(2.0, l) * Const_arrays::fact[l]);

        lima  = l - am;
        lima  = (lima % 2 == 0) ? lima / 2 : (lima - 1) / 2;
        lim1b = lx - am;
        lim1b = (lim1b % 2 == 0) ? lim1b / 2 : (lim1b + 1) / 2;
        lim1b = std::max(lim1b, 0);
        lim2b = (lx % 2 == 0) ? lx / 2 : (lx - 1) / 2;
        lim2b = std::min(lim2b, bs);

        val = 0.0;
        for (int a = std::max(bs, 0); a <= lima; a++) {
            twoa = 2 * a;
            faci = Const_arrays::binom[l][a] * Const_arrays::belt[a] * Const_arrays::fact[twol - twoa] * Const_arrays::binom[a][bs] / Const_arrays::fact[l - am - twoa];
            for (int b = lim1b; b <= lim2b; b++) {
                twob = 2 * b + am - lx + 1;
                if (twob % 2 == 1)
                    continue;
                val += Const_arrays::binom[bs][b] * Const_arrays::binom[am][twob - 1] * Const_arrays::belt[twob / 2] * faci;
            }
        }
        return Const_arrays::belt[am] * val * fac / std::sqrt(2.0 * M_PI);
    }
    return 0;
}

double Harmonics::NoNormYrl(const int &l, const int &m, const Eigen::Vector3d &r) {
    double res = 0;

    for (int lx = 0; lx <= l; ++lx)
        for (int ly = 0; ly <= l - lx; ++ly) {
            int lz = l - lx - ly;
            res += NoNormCalcClmR(l, m, lx, ly, lz) * std::pow(r(0), lx) * std::pow(r(1), ly) * std::pow(r(2), lz);
        }
    return res;
}

double Harmonics::Y(const int &l, const int &m, const Eigen::Vector3d &r) {
    return NoNormYrl(l, m, r) / std::pow(r.norm(), l);
}