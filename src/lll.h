#ifndef __LLL_H__
#define __LLL_H__

#include <vector>
#include <cmath>

class LLL {

public:
    LLL(std::vector<std::vector<double>> __b);

    std::vector<std::vector<double>> solve(const double __delta);

private:
    const int _m;
    const int _n;

    const std::vector<std::vector<double>> _b;

    void gram_schmidt(std::vector<std::vector<double>>::const_iterator __first1,
                      std::vector<std::vector<double>> &__out);
    void update_mu(std::vector<std::vector<double>>::const_iterator __first1,
                   std::vector<std::vector<double>>::const_iterator __last1,
                   std::vector<std::vector<double>>::const_iterator __first2,
                   std::vector<std::vector<double>>::const_iterator __last2,
                   std::vector<std::vector<double>> &__out);
    double inner_product(std::vector<double>::const_iterator __first1,
                         std::vector<double>::const_iterator __last1,
                         std::vector<double>::const_iterator __first2);
    std::vector<double> vector_sub(std::vector<double>::const_iterator __first1,
                                   std::vector<double>::const_iterator __last1,
                                   std::vector<double>::const_iterator __first2);
    std::vector<double> scalar_mult(const double __x,
                                    std::vector<double>::const_iterator __first1,
                                    std::vector<double>::const_iterator __last1);
};

#endif
