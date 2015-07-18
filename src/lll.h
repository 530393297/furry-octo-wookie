#ifndef __LLL_H__
#define __LLL_H__

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>

double inner_product(std::vector<double>::iterator first1, std::vector<double>::iterator last1, std::vector<double>::iterator first2);
void gram_schmidt(std::vector<std::vector<double>>::iterator a, int, int, std::vector<std::vector<double>> &__first2);
std::vector<double> scalar_mult(double a, std::vector<double> b);
std::vector<double> vector_sub(std::vector<double> a, std::vector<double> b);
std::vector<double> scalar_div(double a, std::vector<double> b);
std::vector<std::vector<double>> LLL(std::vector<std::vector<double>> input, double delta);
std::vector<double> vector_add(std::vector<double> a, std::vector<double> b);
std::vector<double> matrix_mult(std::vector<std::vector<double>> a, std::vector<double> b);

#endif
