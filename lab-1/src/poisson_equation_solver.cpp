//
// Created by plush-jill on 3/4/25.
//

#include "poisson_equation_solver.hpp"


float PoissonEquationSolver::calc_new_value(const float F_im1_jm1, const float F_im1_j, const float F_im1_jp1,
                                            const float F_i_jm1,                        const float F_i_jp1,
                                            const float F_ip1_jm1, const float F_ip1_j, const float F_ip1_jp1,

                                            const float P_im1_j,
                                            const float P_i_jm1,   const float P_i_j,   const float P_i_jp1,
                                            const float P_ip1_j
) const {

    const float first_line = m_b * (F_i_jm1 + F_i_jp1) + m_b * (F_im1_j + F_ip1_j);
    const float second_line = m_c * (F_im1_jm1 + F_im1_jp1 + F_ip1_jm1 + F_ip1_jp1);
    const float third_line = (2.0f * P_i_j + (1.0f / 4.0f) * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1));
    const float result = m_a * (first_line + second_line + third_line);

    return result;
};

void PoissonEquationSolver::solve() {

    for (int time {}; time < m_Nt; ++time) {
        float delta = INT_MIN;
        std::vector<std::vector<std::pair<float, float>>> prev_grid = m_grid;
        for (int i {1}; i < m_Ny - 1; ++i) {
            for (int j {1}; j < m_Nx - 1; ++j) {
                m_grid[i][j].second = calc_new_value(prev_grid[i-1][j-1].second, prev_grid[i-1][j].second, prev_grid[i-1][j+1].second,
                                              prev_grid[i][j-1].second,                                     prev_grid[i][j+1].second,
                                              prev_grid[i+1][j-1].second,   prev_grid[i+1][j].second, prev_grid[i+1][j+1].second,

                                                                         m_grid[i-1][j].first,
                                              m_grid[i][j-1].first, m_grid[i][j].first, m_grid[i][j+1].first,
                                                                          m_grid[i+1][j].first
                                              );

                delta = std::max(delta, std::abs(prev_grid[i][j].second - m_grid[i][j].second));

            }
        }
        m_deltas[time] = delta;
    }



}

void PoissonEquationSolver::export_grid_value_as_matrix(const std::string &file_path) const {
    std::ofstream output (file_path, std::ios::out);
    if (!output.is_open()) {
        std::cerr << "Failed to open file: " << "float2.dat" << std::endl;
        throw std::exception();
    }


    for (int i {}; i < m_grid.size(); ++i) {
        for (int j {}; j < m_grid[i].size(); ++j) {
            output << m_grid[i][j].second << " ";
        }
        output << "\n";
    }

    output.close();
}


void PoissonEquationSolver::check_deltas() const {
    int k {};
    for (int time {1}; time < m_Nt; ++time) {
        if (!(m_deltas[time] <= m_deltas[time - 1])) {
            std::cerr << std::format("iter {}: delta >= previous delta ({} >= {}).", time, m_deltas[time], m_deltas[time - 1]) << std::endl;
            ++k;
        }
    }
    std::cerr << std::format("delta >= previous delta on {}/{} iterations", k, m_Nt - 1) << std::endl;


}