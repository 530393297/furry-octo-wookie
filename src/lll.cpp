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

    std::vector<std::vector<double>> B(n, std::vector<double>(input[0].size()));
    gram_schmidt(input.begin(), input.size(), input[0].size(), B);
    std::vector<std::vector<double>> l(n, std::vector<double>(input[0].size()));
    

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < input[0].size(); j++) {
            l[i][j] = inner_product(input[i].begin(), input[i].end(), B[j].begin()) / inner_product(B[j].begin(), B[j].end(), B[j].begin());
        }
    }


    int k = 1;
    while (k < n) {
        for(int j = k - 1; j >= 0; j--) {
            if(fabs(l[k][j]) > 0.5) {
                input[k] = vector_sub(input[k], scalar_mult(std::round(l[k][j]), input[j]));
                std::vector<double> tt = scalar_mult(std::round(l[k][j]), input[j]);
        
            	gram_schmidt(input.begin(), input.size(), input[0].size(), B);
                for(int ii = 0; ii < n; ii++) {
                    for(int jj = 0; jj < input[0].size(); jj++) {
                        l[ii][jj] = inner_product(input[ii].begin(), input[ii].end(), B[jj].begin()) / inner_product(B[jj].begin(), B[jj].end(), B[jj].begin());
                    }
                }
            }
        }
        if(inner_product(B[k].begin(), B[k].end(), B[k].begin()) >= (delta - l[k][k-1] * l[k][k-1]) * inner_product(B[k -1].begin(), B[k -1].end(), B[k -1].begin())) {
            k++;
        } else {
            input[k].swap(input[k -1]);
            gram_schmidt(input.begin(), input.size(), input[0].size(), B);
            for(int ii = 0; ii < n; ii++) {
                for(int jj = 0; jj < input[0].size(); jj++) {
                    l[ii][jj] = inner_product(input[ii].begin(), input[ii].end(), B[jj].begin()) / inner_product(B[jj].begin(), B[jj].end(), B[jj].begin());
                }
            }
            k = std::max(k - 1, 1);
        }
    }
    return input;
}


void gram_schmidt(std::vector<std::vector<double>>::iterator __begin, int __n, int __m, std::vector<std::vector<double>> &__first2)
{
    auto first = __begin;    

    for(int i = 0; i < __n; i++, ++__begin) {
        __first2[i] = *__begin;
        for(auto inner = first; inner != __begin; ++inner) 
            __first2[i] = vector_sub(*__begin, scalar_mult(inner_product(std::begin(*__begin), std::end(*__begin), std::begin(*inner)) / inner_product(std::begin(*inner), std::end(*inner), std::begin(*inner)), *inner));
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

std::vector<double> scalar_mult(double a, std::vector<double> b)
{
    std::vector<double> ans (b.size(), 0);
    for(int i = 0; i < b.size(); i++) {
        ans[i] = a * b[i];
    }
    return ans;
}

std::vector<double> scalar_div(double a, std::vector<double> b)
{
    std::vector<double> ans (b.size(), 0);
    for(int i = 0; i < b.size(); i++) {
        ans[i] = b[i] / a;
    }
    return ans;
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
