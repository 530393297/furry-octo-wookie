#include "lll.h"

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
}

/*std::vector<std::vector<double>> LLL(std::vector<std::vector<double>> input, double delta)
{
    const int n = input.size();
    std::vector<std::vector<double>> b(n, std::vector<double>(input[0].size()));
    b = size_reduce(input);
    std::vector<std::vector<double>> B = gram_schmidt(b);

    for(int i = 0; i < n - 1; i++) {
        std::vector<double> project = projection(B, b[i], i);
        std::vector<double> projecttwo = projection(B, b[i + 1], i);
        if(delta * inner_product(project.begin(), project.end(), project.begin()) > inner_product(projecttwo.begin(), projecttwo.end(), projecttwo.begin())) {
            for(int ii = 0; ii < 3; ii++)
                std::cout << b[i][ii] << ",";
            std::cout << std::endl;
            for(int ii = 0; ii < 3; ii++)
                std::cout << b[i+1][ii] << ",";
            std::cout << std::endl;
            b[i].swap(b[i +1]);
            for(int ii = 0; ii < 3; ii++)
                std::cout << b[i][ii] << ",";
            std::cout << std::endl;
            for(int ii = 0; ii < 3; ii++)
                std::cout << b[i+1][ii] << ",";
            std::cout << std::endl;
            return LLL(b, delta);
        }
    }
    return b;
} */

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
        std::cout << k << std::endl;
        for(int j = k - 1; j >= 0; j--) {
            std::cout << j << "," << fabs(l[k][j]) << std::endl;
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
    std::vector<std::vector<double>> q(n, std::vector<double>(a[0].size(), 0));
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

std::vector<double> nearest_plane(std::vector<std::vector<double>> b, std::vector<double> w) {

    const int n = b.size();
    
    std::vector<std::vector<double>> W(n, std::vector<double>(b[0].size()));
    std::vector<std::vector<double>> y(n, std::vector<double>(b[0].size()));

    std::vector<std::vector<double>> B = gram_schmidt(b);

    W[n -1] = w;

    std::vector<double> l(n);
    for(int i = n - 1; i > 0; i--) {
        l[i] = inner_product(W[i].begin(), W[i].end(), B[i].begin()) / inner_product(B[i].begin(), B[i].end(), B[i].begin());
        y[i] = scalar_mult(std::round(l[i]), b[i]);
        W[i - 1] = vector_sub(vector_sub(W[i], scalar_mult(l[i] - std::round(l[i]), B[i])), scalar_mult(std::round(l[i]), b[i]) );
    }

    l[0] = inner_product(W[0].begin(), W[0].end(), B[0].begin()) / inner_product(B[0].begin(), B[0].end(), B[0].begin());
    y[0] = scalar_mult(std::round(l[0]), b[0]);

    std::vector<double> ans(n, 0);
    for(int i = 0; i < n; i++) {
        ans = vector_add(ans, y[i]);
    }
    return ans;
}

std::vector<std::vector<double>> size_reduce(std::vector<std::vector<double>> b) 
{
    std::vector<double> x(b[0].size());
    std::vector<std::vector<double>> B = gram_schmidt(b);

    for(int i = 1; i < b.size(); i++) {
        x = nearest_plane(b, vector_sub(b[i], B[i]));
        b[i] = vector_sub(b[i], matrix_mult(b, x));
    }
    return b;
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

std::vector<double> projection(std::vector<std::vector<double>> B, std::vector<double> x, int i) {
    const int n = x.size(); 
    std::vector<double> p(n, 0);
    for(int j = i; j < n; j++) {
        p = vector_add(scalar_mult(inner_product(x.begin(), x.end(), B[j].begin()) / inner_product(B[j].begin(), B[j].end(), B[j].begin()), B[j]), p);
    }
    return p;
}