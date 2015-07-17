#ifndef __LLL_H__
#define __LLL_H__

#include <vector>
#include <iostream>
#include <numeric>

double dot_prod(std::vector<double> a, std::vector<double> b);
std::vector<std::vector<double>> grant_schmidt(std::vector<std::vector<double>> a);
std::vector<double> scalar_mult(double a, std::vector<double> b);
std::vector<double> vector_sub(std::vector<double> a, std::vector<double> b);
std::vector<double> scalar_div(double a, std::vector<double> b);
void LLL(std::vector<std::vector<double>> input);

#endif