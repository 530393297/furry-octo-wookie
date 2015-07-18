#include "lll.h"
#include "polynomial.h"

int main()
{
    std::vector<std::vector<double>> input = {	{1, 1, 1},
        { -1, 0, 2},
        {3, 5, 6}
    };

    std::vector<std::vector<double>>  out = LLL(input, 0.75);
    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << out[j][i] << ",";
        }
        std::cout << std::endl;
    }
    
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1, 100);
    std::vector<std::vector<double>> q(45, std::vector<double>(45));
    for (int i = 0 ; i < 45; ++i)
        for(int j = 0; j < 45; ++j)
            q[i][j] = dist(mt);
    LLL(q, 0.75); 

    std::vector<int> polya = {2, 5, 7};
    std::vector<int> polyb = {4,7, 5};

    Polynomial a(polya);
    Polynomial b(polyb);

    Polynomial c = a.mod(2);

    std::vector<int> polyc = c.get_coefficients();
    for(int i = 0; i < polyc.size(); i++) 
        std::cout << polyc[i] << ",";

}

std::vector<std::vector<double>> LLL(std::vector<std::vector<double>> input, double delta) {
    const int n = input.size();
    const int m = input[0].size();

    std::vector<std::vector<double>> B(n, std::vector<double>(m));
    std::vector<std::vector<double>> mu(n, std::vector<double>(m));
    
    gram_schmidt(std::begin(input), B);
    update_mu(std::begin(input), std::end(input), std::begin(B), std::end(B), mu);

    auto k = 1;
    while (k < n) {
        for(int j = k - 1; j >= 0; j--) {
            if(std::fabs(mu[k][j]) > 0.5) {
                input[k] = vector_sub(input[k], scalar_mult(std::round(mu[k][j]), std::begin(input[j]), std::end(input[j])));
            	gram_schmidt(std::begin(input), B);
                update_mu(std::begin(input), std::end(input), std::begin(B), std::end(B), mu);
            }
        }
        if(inner_product(std::begin(B[k]), std::end(B[k]), std::begin(B[k])) >= (delta - mu[k][k-1] * mu[k][k-1]) * inner_product(B[k -1].begin(), B[k -1].end(), B[k -1].begin())) {
            k++;
        } else {
            input[k].swap(input[k -1]);
            gram_schmidt(std::begin(input), B);
            update_mu(std::begin(input), std::end(input), std::begin(B), std::end(B), mu);
            k = std::max(k - 1, 1);
        }
    }
    return input;
}

void update_mu(std::vector<std::vector<double>>::iterator __first1,
               std::vector<std::vector<double>>::iterator __last1,
               std::vector<std::vector<double>>::iterator __first2,
               std::vector<std::vector<double>>::iterator __last2,
               std::vector<std::vector<double>> &__out) 
{
    auto j = 0;
    
    for(auto i = 0; __first1 != __last1; i++, ++__first1) {
	j = 0;
        for(auto inner = __first2; inner != __last2; j++, ++inner) {
            auto begin = std::begin(*inner);
            auto end = std::end(*inner); 
            __out[i][j] = inner_product(std::begin(*__first1), std::end(*__first1), begin) / inner_product(begin, end, begin);
        }
    }
}

void gram_schmidt(std::vector<std::vector<double>>::iterator __first1, std::vector<std::vector<double>> &__out)
{
    auto first = __first1;    
    int n = __out.size();

    for(auto i = 0; i < n; i++, ++__first1) {
        __out[i] = *__first1;
        for(auto inner = first; inner != __first1; ++inner) { 
	    auto begin = std::begin(*inner);
            auto end = std::end(*inner);
            __out[i] = vector_sub(*__first1, scalar_mult(inner_product(std::begin(*__first1), std::end(*__first1), begin) / inner_product(begin, end, begin), begin, end));
        }
    }
}

double inner_product(std::vector<double>::iterator __first1,
                     std::vector<double>::iterator __last1,
                     std::vector<double>::iterator __first2)
{
    auto init = 0.0;
   
     for (; __first1 != __last1; ++__first1, ++__first2)
        init = init + (*__first1 * *__first2);

    return init;
}

std::vector<double> scalar_mult(double __x, 
                                std::vector<double>::iterator __first1,
                                std::vector<double>::iterator __last1)
{
    std::vector<double> out (std::distance(__first1, __last1), 0);
    
    for(auto i = 0; __first1 != __last1; i++, ++__first1) 
        out[i] = __x * *__first1;
    
    return out;
}

std::vector<double> scalar_div(double __x,
                                std::vector<double>::iterator __first1,
                                std::vector<double>::iterator __last1)
{
    std::vector<double> out (std::distance(__first1, __last1), 0);

    for(auto i = 0; __first1 != __last1; i++, ++__first1)
        out[i] = __x / *__first1;

    return out;
}

std::vector<double> vector_sub(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> ans (b.size(), 0);
    for(int i = 0; i < b.size(); i++) {
        ans[i] = a[i] - b[i];
    }
    return ans;
}

std::vector<double> vector_add(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> ans (b.size(), 0);
    for(int i = 0; i < b.size(); i++) {
        ans[i] = a[i] + b[i];
    }
    return ans;
}

std::vector<double> matrix_mult(std::vector<std::vector<double>> a, std::vector<double> b) 
{
    std::vector<double> c(a[0].size(), 0);
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < a[0].size(); j++) {
            c[j] += a[i][j] * b[i]; 
        }
    }
    return c;
}
