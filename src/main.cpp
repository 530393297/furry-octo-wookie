#include "lll.h"

#include <iostream>
#include <random>

int main(int argc, char* argv[])
{
    const int n = std::atoi(argv[1]);

    // Wikipedia test
    std::vector<std::vector<double>> input = {
        {1, 1, 1},
        { -1, 0, 2},
        {3, 5, 6}
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
}
