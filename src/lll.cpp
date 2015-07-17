#include "lll.h"

int main()
{
  std::vector<std::vector<double>> input = {	{1,1,1},
                               					{-1,0,2},
                               					{3,5,6}};
  for(int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << input[j][i] << ",";
		}
		std::cout << std::endl;
	}                             					
  LLL(input);
}

void LLL(std::vector<std::vector<double>> input) 
{
	std::vector<std::vector<double>> ans = grant_schmidt(input);
	for(int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << ans[j][i] << ",";
		}
		std::cout << std::endl;
	}
}

std::vector<std::vector<double>> grant_schmidt(std::vector<std::vector<double>> a) 
{
	const int n = a.size();
	std::vector<std::vector<double>> r(n, std::vector<double>(a[0].size()));
	std::vector<std::vector<double>> q(n, std::vector<double>(a[0].size(), 0));
	std::vector<std::vector<double>> v(n, std::vector<double>(a[0].size(), 0));
	for(int i = 0; i < n; i++) 
	{
		v[i] = a[i];
		for(int j = 0; j < i; j++) {
			r[i][j] = dot_prod(a[i], v[j])/dot_prod(v[j], v[j]);
			v[i] = vector_sub(v[i], scalar_mult(r[i][j], v[j]));
		}
		//r[j][j] = dot_prod(v[j], v[j]);
		//q[j] = scalar_div(r[j][j], v[j]);
	}
	return v;
}

// FIX ME
double dot_prod(std::vector<double> a, std::vector<double> b) 
{
	double ans = 0;
	for(int i = 0; i < a.size(); i++) {
		ans += a[i] * b[i];
	}
	return ans;
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