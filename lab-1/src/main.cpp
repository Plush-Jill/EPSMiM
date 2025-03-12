#include "poisson_equation_solver.hpp"
#include <nsimd/nsimd-all.h>
#include <iostream>

int main() {
    PoissonEquationSolver solver {"../config.json"};
    solver.solve();
    solver.export_grid_value_as_matrix("../float2.dat");
    solver.check_deltas();


    return 0;
}