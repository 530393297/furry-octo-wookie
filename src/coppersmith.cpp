#include "coppersmith.h"

CopperSmith::CopperSmith(Polynomial __input, const int __m, const int __n) : _input(__input),
    _m(__m), _n(__n), _d(__input.get_degree()), _polys(__input.get_degree(), std::vector<Polynomial>(__m + 1, Polynomial((__m + 1) * __input.get_degree())))
{

}

std::vector<double> CopperSmith::run_attack() 
{
	generate_polys();
	std::vector<std::vector<long>> ilattice = generate_lattice();
	std::vector<std::vector<double>> dlattice(ilattice.size(), std::vector<double>(ilattice[0].size()));
	
	for(int i = 0; i < ilattice.size(); i++){
		for(int j = 0; j < ilattice[0].size(); j++) {
			dlattice[i][j] = (double)(ilattice[i][j] % _n);
			//std::cout << dlattice[i][j] << ",";
		}
		//
		//std::cout << std::endl;
	}

	LLL l(dlattice);
	std::vector<std::vector<double>> out = l.solve(0.75);
	
	return out[0];
} 

void CopperSmith::generate_polys()
{
	Polynomial x(std::vector<long>({0, 1}));

    for(auto u = 0; u < _d; u++) {
        for(auto v = 0; v <= _m; v++) {
            _polys[u][v] = _input.exp(v) * x.exp(u);
            _polys[u][v] = _polys[u][v].scalar_mult(std::pow(_n, _m - v));
            _polys[u][v].pad(_d * (_m + 1));
            _polys[u][v].mod(_n);
        }
    }
}

std::vector<std::vector<long>> CopperSmith::generate_lattice() 
{
	const int n = _d * (_m + 1);

	const int B = std::pow(std::pow(_n, 1 / _d), _m / (_m + 1));
	std::vector<std::vector<long>> out(n, std::vector<long>(n, 0));

	for(auto u = 0; u < _d; u++) {
        for(auto v = 0; v <= _m; v++) {
			out[u + v * _d] = _polys[u][v].partial_evaluate(B).get_coefficients();
			std::cout << "!!" << out[u + v * _d].size() <<  std::endl;
		}
	}

	return out;
}