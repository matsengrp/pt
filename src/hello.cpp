#include <iostream>
#include <vector>

#include <libpll/pll.h>


int main() {

    std::vector<double> v(4);

    pll_compute_gamma_cats(0.3, 4, &v[0]);

    std::cout << v[1] << std::endl;

}
