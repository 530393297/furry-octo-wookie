#ifndef __LLL_H__
#define __LLL_H__

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>

double inner_product(std::vector<double>::const_iterator first1, 
	                 std::vector<double>::const_iterator last1, 
	                 std::vector<double>::const_iterator first2);
void gram_schmidt(std::vector<std::vector<double>>::const_iterator first1, 
	              std::vector<std::vector<double>> &first2);
void update_mu(std::vector<std::vector<double>>::const_iterator __first1, 
	           std::vector<std::vector<double>>::const_iterator __last1, 
	           std::vector<std::vector<double>>::const_iterator __first2,
               std::vector<std::vector<double>>::const_iterator __last2, 
               std::vector<std::vector<double>> &__out); 

std::vector<double> scalar_mult(double __x,
                                std::vector<double>::const_iterator __first1,
                                std::vector<double>::const_iterator __last1);
std::vector<double> scalar_div(double __x,
                                std::vector<double>::iterator __first1,
                                std::vector<double>::iterator __last1);
std::vector<double> vector_sub(std::vector<double>::const_iterator __first1,
                               std::vector<double>::const_iterator __last1,
                               std::vector<double>::const_iterator __first2);

void LLL(std::vector<std::vector<double>> &input, double delta);


#endif
