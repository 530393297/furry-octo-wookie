#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <vector>
#include <iostream>
#include <cmath>

class Polynomial {

public:
    Polynomial(int _degree);
    Polynomial(std::vector<long> coefficients);
    Polynomial(const Polynomial &poly);

    Polynomial operator+ (const Polynomial& rhs);
    Polynomial operator* (const Polynomial& rhs);

    Polynomial exp (int exp);
    Polynomial mod (int n);

    int get_degree();
    std::vector<long> get_coefficients() const;

    Polynomial scalar_mult (const int __x);
    Polynomial partial_evaluate (const int __b);
    long evaluate(const int __x);
    void pad (const int __x);

private:
    int degree;
    std::vector<long> coefficients;

};

#endif