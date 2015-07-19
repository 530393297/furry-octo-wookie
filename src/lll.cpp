#include "lll.h"

LLL::LLL(std::vector<std::vector<double>> __b) : _m(__b[0].size()), _n(__b.size()), _b(__b)
{
}

std::vector<std::vector<double>> LLL::solve(const double __delta)
{
    std::vector<std::vector<double>> B(_n, std::vector<double>(_m));
    std::vector<std::vector<double>> mu(_n, std::vector<double>(_m));
    std::vector<std::vector<double>> b = _b;

    gram_schmidt(std::begin(b), B);
    update_mu(std::begin(b), std::end(b), std::begin(B), std::end(B), mu);

    auto k = 1;
    while (k < _n) {
        for(auto j = k - 1; j >= 0; j--) {
            if(std::fabs(mu[k][j]) > 0.5) {
                b[k] = vector_sub(std::begin(b[k]), std::end(b[k]),
                                  std::begin(scalar_mult(
                                                 std::round(mu[k][j]),
                                                 std::begin(b[j]),
                                                 std::end(b[j]))));
                gram_schmidt(std::begin(b), B);
                update_mu(std::begin(b), std::end(b), std::begin(B),
                          std::end(B), mu);
            }
        }
        if(inner_product(
                    std::begin(B[k]), std::end(B[k]), std::begin(B[k])) >=
                (__delta - mu[k][k - 1] * mu[k][k - 1]) *
                inner_product(std::begin(B[k - 1]), std::end(B[k - 1]),
                              std::begin(B[k - 1]))) {
            k++;
        } else {
            b[k].swap(b[k - 1]);
            gram_schmidt(std::begin(b), B);
            update_mu(std::begin(b), std::end(b), std::begin(B), std::end(B),
                      mu);
            k = std::max(k - 1, 1);
        }
    }
    return b;
}

void LLL::gram_schmidt(std::vector<std::vector<double>>::const_iterator __first1,
                       std::vector<std::vector<double>> &__out)
{
    int n = __out.size();
    auto first = __first1;

    for(auto i = 0; i < n; i++, ++__first1) {
        __out[i] = *__first1;
        for(auto inner = first; inner != __first1; ++inner) {
            auto begin1 = std::begin(*__first1);
            auto begin2 = std::begin(*inner);
            auto end1 = std::end(*__first1);
            auto end2 = std::end(*inner);
            __out[i] = vector_sub(begin1, end1, std::begin(
                                      scalar_mult(
                                          inner_product(begin1, end1, begin2) /
                                          inner_product(begin2, end2, begin2),
                                          begin2, end2)));
        }
    }
}

void LLL::update_mu(std::vector<std::vector<double>>::const_iterator __first1,
                    std::vector<std::vector<double>>::const_iterator __last1,
                    std::vector<std::vector<double>>::const_iterator __first2,
                    std::vector<std::vector<double>>::const_iterator __last2,
                    std::vector<std::vector<double>> &__out)
{
    for(auto i = 0; __first1 != __last1; i++, ++__first1) {
        auto j = 0;
        for(auto inner = __first2; inner != __last2; j++, ++inner) {
            auto begin = std::begin(*inner);
            auto end = std::end(*inner);
            __out[i][j] = inner_product(std::begin(*__first1),
                                        std::end(*__first1), begin) /
                          inner_product(begin, end, begin);
        }
    }
}

double LLL::inner_product(std::vector<double>::const_iterator __first1,
                          std::vector<double>::const_iterator __last1,
                          std::vector<double>::const_iterator __first2)
{
    auto init = 0.0;

    for (; __first1 != __last1; ++__first1, ++__first2)
        init = init + (*__first1 **__first2);

    return init;
}

std::vector<double> LLL::vector_sub(std::vector<double>::const_iterator __first1,
                                    std::vector<double>::const_iterator __last1,
                                    std::vector<double>::const_iterator __first2)
{
    std::vector<double> out(std::distance(__first1, __last1), 0);

    for(auto begin = std::begin(out);  __first1 != __last1;
            ++__first1, ++__first2, ++begin)
        *begin = *__first1 - *__first2;

    return out;
}

std::vector<double> LLL::scalar_mult(const double __x,
                                     std::vector<double>::const_iterator __first1,
                                     std::vector<double>::const_iterator __last1)
{
    std::vector<double> out(std::distance(__first1, __last1), 0);

    for(auto begin = std::begin(out); __first1 != __last1; ++__first1, ++begin)
        *begin = __x **__first1;

    return out;
}