#ifndef HARMONICS_H
#define HARMONICS_H

#include "constants.h"

#include <cmath>
#include <algorithm>

class Harmonics
{
public:
    static double NoNormCalcClmR( const int l, const int m, const int lx, const int ly, const int lz );
};

#endif // HARMONICS_H
