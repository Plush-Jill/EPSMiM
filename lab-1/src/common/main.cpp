#ifdef LAB_1
#include "poisson_equation_solver.hpp"
#elif defined(LAB_2)
#include "../lab-2/poisson_equation_solver.hpp"
#elif defined(LAB_3)
#include "../lab-3/poisson_equation_solver.hpp"
#endif


int main(int argc, char *argv[]) {

    // const auto config_path = std::string(argv[1]);


    // PoissonEquationSolver solver {config_path};
    PoissonEquationSolver solver {"/home/plush-jill/All_Random/git/EPSMiM/lab-1/config.json"};

    solver.solve();
    // solver.export_grid_value_as_matrix("/home/plush-jill/All_Random/git/EPSMiM/lab-1/float2.dat");
    // solver.check_deltas();

    return EXIT_SUCCESS;
}