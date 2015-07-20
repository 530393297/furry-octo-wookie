#include "coppersmith.h"

#include <iostream>
#include <random>

int main(int, char** argv)
{
    const int n = std::atoi(argv[1]);

    // Wikipedia test
    std::vector<std::vector<double>> input = {
        {1, 1, 1},
        {-1, 0, 2},
        {3, 5, 6},
    };

    auto l_real = LLL(input);
    std::vector<std::vector<double>> ans = l_real.solve(0.75);

    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            std::cout << ans[j][i] << ",";
        std::cout << std::endl;
    } 

    // Large random test
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1, 100);
    std::vector<std::vector<double>> rand(n, std::vector<double>(n));
    for (int i = 0 ; i < n; ++i)
        for(int j = 0; j < n; ++j)
            rand[i][j] = dist(mt);

    auto l_rand = LLL(rand);
    ans = l_rand.solve(0.75);  

    /*std::vector<long> coeff({4, 0, 0,0,0,0,0, 1});
    Polynomial p(coeff);
    CopperSmith cs(p, 10, 33);

    auto ans = cs.run_attack();
    std::vector<long> ians(ans.size());
    for(int i = 0; i < ans.size(); i++) {
        ians[i] = std::round(ans[i]);
        std::cout << ians[i]<< ",";
    }
    std::cout << std::endl;
    Polynomial pp(ians);

    for(int i = 1; i <33; i++)
        std::cout << pp.evaluate(i) % 33  << std::endl;*/
}
