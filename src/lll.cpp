#include "lll.h"
#include "polynomial.h"

int main()
{
    std::vector<std::vector<double>> input = {	{1, 1, 1},
        { -1, 0, 2},
        {3, 5, 6}
    };

    LLL(input, 0.75);
    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << input[j][i] << ",";
        }
        std::cout << std::endl;
    }


    std::mt19937 mt(1000000000000);
    std::uniform_real_distribution<double> dist(1, 100);
    std::vector<std::vector<double>> q(100, std::vector<double>(100));
    for (int i = 0 ; i < 100; ++i)
        for(int j = 0; j < 100; ++j)
            q[i][j] = dist(mt);
    LLL(q, 0.75);
}


void LLL(std::vector<std::vector<double>> &b, double delta)
{
    const int n = b.size();
    const int m = b[0].size();

    std::vector<std::vector<double>> B(n, std::vector<double>(m));
    std::vector<std::vector<double>> mu(n, std::vector<double>(m));

    gram_schmidt(std::begin(b), B);
    update_mu(std::begin(b), std::end(b), std::begin(B), std::end(B), mu);

    auto k = 1;
    while (k < n) {
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
                (delta - mu[k][k - 1] * mu[k][k - 1]) *
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
}

void update_mu(std::vector<std::vector<double>>::const_iterator __first1,
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

void gram_schmidt(std::vector<std::vector<double>>::const_iterator __first1,
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

double inner_product(std::vector<double>::const_iterator __first1,
                     std::vector<double>::const_iterator __last1,
                     std::vector<double>::const_iterator __first2)
{
    auto init = 0.0;

    for (; __first1 != __last1; ++__first1, ++__first2)
        init = init + (*__first1 **__first2);

    return init;
}

std::vector<double> scalar_mult(double __x,
                                std::vector<double>::const_iterator __first1,
                                std::vector<double>::const_iterator __last1)
{
    std::vector<double> out(std::distance(__first1, __last1), 0);

    for(auto i = 0; __first1 != __last1; i++, ++__first1)
        out[i] = __x **__first1;

    return out;
}

std::vector<double> vector_sub(std::vector<double>::const_iterator __first1,
                               std::vector<double>::const_iterator __last1,
                               std::vector<double>::const_iterator __first2)
{
    std::vector<double> out(std::distance(__first1, __last1), 0);

    for(int i = 0;  __first1 != __last1; i++, ++__first1, ++__first2)
        out[i] = *__first1 - *__first2;

    return out;
}