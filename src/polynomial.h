#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <vector>
#include <iostream>

class Polynomial {

public:
	Polynomial(int _degree);
	Polynomial(std::vector<int> coefficients);
	Polynomial(const Polynomial &poly);

	Polynomial operator+ (const Polynomial& rhs);
	Polynomial operator* (const Polynomial& rhs);

	Polynomial exp (int exp);
	Polynomial mod (int n); 

	int get_degree() {return degree;};
	std::vector<int> get_coefficients() const;



private:
	int degree;
	std::vector<int> coefficients;

};

#endif