#ifndef __COPPERSMITH_H__
#define __COPPERSMITH_H__

#include "lll.h"
#include "polynomial.h"

class CopperSmith {

public:
    CopperSmith(Polynomial __input, const int __m, const int __n);

    std::vector<double> run_attack();

private:
    Polynomial _input;
    const int _m;
    const int _n;
    const int _d;

    std::vector<std::vector<Polynomial>> _polys;

    void generate_polys();
    std::vector<std::vector<long>> generate_lattice();
};

#endif