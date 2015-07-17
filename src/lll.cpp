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
    /*std::vector<std::vector<double>> q(40, std::vector<double>(40));
    for (int i = 0 ; i < 40; ++i)
        for(int j = 0; j < 40; ++j)
            q[i][j] = dist(mt);
    LLL(q, 0.75); */

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

    std::vector<std::vector<double>> B = gram_schmidt(input);
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
        
                B = gram_schmidt(input);
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
            B = gram_schmidt(input);
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


std::vector<std::vector<double>> gram_schmidt(std::vector<std::vector<double>> a)
{
    const int n = a.size();
    std::vector<std::vector<double>> r(n, std::vector<double>(a[0].size()));
    std::vector<std::vector<double>> v(n, std::vector<double>(a[0].size(), 0));
    for(int i = 0; i < n; i++) {
        v[i] = a[i];
        for(int j = 0; j < i; j++) {
            r[i][j] = inner_product(a[i].begin(), a[i].end(), v[j].begin()) / inner_product(v[j].begin(), v[j].end(), v[j].begin());
            v[i] = vector_sub(v[i], scalar_mult(r[i][j], v[j]));
        }
    }
    return v;
}

double inner_product(std::vector<double>::iterator first1,
                     std::vector<double>::iterator last1,
                     std::vector<double>::iterator first2)
{
    double init = 0;
    for (; first1 != last1; ++first1, ++first2)
        init = init + (*first1 * *first2);

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
