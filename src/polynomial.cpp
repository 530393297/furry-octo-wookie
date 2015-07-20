#include "polynomial.h"

Polynomial::Polynomial(int _degree) : degree(_degree), coefficients(_degree + 1)
{
}

Polynomial::Polynomial(std::vector<long> _coefficients) : degree(_coefficients.size() - 1), coefficients(_coefficients)
{
}

Polynomial::Polynomial(const Polynomial &poly) : Polynomial(poly.get_coefficients())
{
}

std::vector<long> Polynomial::get_coefficients() const
{
    return coefficients;
}

Polynomial Polynomial::operator+ (const Polynomial& rhs)
{
    int coeff = 0;

    if(degree < rhs.degree) {
        std::vector<long> res(rhs.degree);

        for (; coeff <= degree; coeff++)
            res[coeff] = coefficients[coeff] + rhs.coefficients[coeff];
        for (; coeff <= rhs.degree; coeff++)
            res[coeff] = rhs.coefficients[coeff];
        return Polynomial(res);
    } else {
        std::vector<long> res(degree);
        for (; coeff <= rhs.degree; coeff++)
            res[coeff] = coefficients[coeff] + rhs.coefficients[coeff];
        for (; coeff <= degree; coeff++)
            res[coeff] = coefficients[coeff];
        return Polynomial(res);
    }
}

Polynomial Polynomial::operator* (const Polynomial& rhs)
{
    std::vector<long> res(degree + rhs.degree + 1, 0);

    for(int i = 0; i <= degree; i++) {
        for(int j = 0; j <= rhs.degree; j++) {
            res[i + j] += coefficients[i] * rhs.coefficients[j];
        }
    }

    return Polynomial(res);
}

Polynomial Polynomial::exp (int exp)
{
    Polynomial y (std::vector<long> (1, 1));
    Polynomial x (coefficients);

    if(exp == 0) {
        return y;
    }
    while(exp > 1) {
        if(exp % 2 == 0) {
            x = x * x;
            exp = exp >> 1;
        } else {
            y = x * y;
            x = x * x;
            exp = (exp - 1) >> 1;
        }
    }
    return x * y;
}

Polynomial Polynomial::mod (int n)
{
    std::vector<long> res(degree + 1, 0);

    for(int i = 0; i <= degree; i++)
        res[i] = coefficients[i] % n;

    return Polynomial(res);
}

Polynomial Polynomial::scalar_mult(const int __x)
{
    std::vector<long> res(degree + 1, 0);

    for(int i = 0; i <= degree; i++)
        res[i] = coefficients[i] * __x;

    return Polynomial(res);
}

int Polynomial::get_degree()
{
    return degree;
}

Polynomial Polynomial::partial_evaluate(const int __b)
{
    std::vector<long> res(degree + 1, 0);

    for(int i = 0; i <= degree; i++)
        res[i] = coefficients[i] * ((long)std::pow(__b, i))%33;

    return Polynomial(res);
}

void Polynomial::pad(const int __x)
{
    coefficients.resize(__x, 0);
    degree = __x;
}

long Polynomial::evaluate(const int __x) 
{
    long res = 0;

    for(int i = 0; i <= degree; i++)
        res += (coefficients[i] * (long)std::pow(__x, i))%33;

    return res;
}