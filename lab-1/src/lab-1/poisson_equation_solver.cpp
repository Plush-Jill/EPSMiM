//
// Created by plush-jill on 3/4/25.
//

#include "../common/poisson_equation_solver.hpp"
#include <fstream>
#include <iostream>

PoissonEquationSolver::PoissonEquationSolver(const std::string &config_file) :
    m_Xa(0), m_Xb(4.0),
    m_Ya(0), m_Yb(4.0) {
    std::ifstream config(config_file);
    std::stringstream json_buffer;
    json_buffer << config.rdbuf();
    boost::json::object json = boost::json::parse(json_buffer.str()).as_object();
    config.close();

    m_Nx = json["Nx"].as_int64();
    m_Ny = json["Ny"].as_int64();
    m_Nt = json["Nt"].as_int64();
    m_export = json["export"].as_bool();

    m_hx = (m_Xb - m_Xa) / static_cast<float>(m_Nx - 1);
    m_hy = (m_Yb - m_Ya) / static_cast<float>(m_Ny - 1);

    float Xs1 = m_Xa + (m_Xb - m_Xa) / 3;
    float Ys1 = m_Ya + (m_Yb - m_Ya) * 2 / 3;
    float Xs2 = m_Xa + (m_Xb - m_Xa) * 2 / 3;
    float Ys2 = m_Ya + (m_Yb - m_Ya) / 3;
    float R = static_cast<float>(0.1) * std::min(m_Xb - m_Xa, m_Yb - m_Ya);

    m_heat_sources.emplace_back(Xs1, Ys1, R, 0.1);
    m_heat_sources.emplace_back(Xs2, Ys2, R, -0.1);

    m_value_grid = std::make_shared<std::vector<std::vector<float>>> (m_Nx, std::vector<float>(m_Ny, 0.0));
    m_previous_value_grid = std::make_shared<std::vector<std::vector<float>>> (m_Nx, std::vector<float>(m_Ny, 0.0));
    m_heat_grid = std::make_shared<std::vector<std::vector<float>>> (m_Nx, std::vector<float>(m_Ny, 0.0));


    m_deltas = std::vector<float> (m_Nt, 0.0);

    for (auto& heat_source : m_heat_sources) {
        for (int j {}; j < m_Ny; ++j) {
            for (int i {}; i < m_Nx; ++i) {

                float x_i = m_Xa + (static_cast<float>(i) * m_hx);
                float y_j = m_Ya + (static_cast<float>(j) * m_hy);
                if (heat_source.is_has_point(x_i, y_j)) {
                    (*m_heat_grid)[i][j] = heat_source.get_heat();
                }
            }
        }
    }

    const float hx_pow_minus1 = 1.0f / m_hx;
    const float hy_pow_minus1 = 1.0f / m_hy;
    hx_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hx, 2));
    hy_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hy, 2));

    m_hx_hy_2 = hx_pow_minus2 + hy_pow_minus2;
    m_a = 1.0f / (5.0f * m_hx_hy_2);
    m_b = static_cast<float>((1.0f / 2.0f) * (5.0f * hy_pow_minus2 - 1.0f * hx_pow_minus2));
    m_c = m_hx_hy_2 / 4;

    // std::cout << "created solver v1" << std::endl;
}


float PoissonEquationSolver::calc_new_value(const float F_im1_jm1,  const float F_im1_j, const float F_im1_jp1,
                                            const float F_i_jm1,                         const float F_i_jp1,
                                            const float F_ip1_jm1,  const float F_ip1_j, const float F_ip1_jp1,

                                                                    const float P_im1_j,
                                            const float P_i_jm1,    const float P_i_j,   const float P_i_jp1,
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
        // std::cout << std::format("time: {}\n", time);
        float delta = INT_MIN;
        std::swap(m_previous_value_grid, m_value_grid);
        for (int i {1}; i < m_Ny - 1; ++i) {
            for (int j {1}; j < m_Nx - 1; ++j) {
                (*m_value_grid)[i][j] = calc_new_value(
                    (*m_previous_value_grid)[i-1][j-1], (*m_previous_value_grid)[i-1][j], (*m_previous_value_grid)[i-1][j+1],
                      (*m_previous_value_grid)[i][j-1],                                              (*m_previous_value_grid)[i][j+1],
                     (*m_previous_value_grid)[i+1][j-1],  (*m_previous_value_grid)[i+1][j], (*m_previous_value_grid)[i+1][j+1],

                                                                 (*m_heat_grid)[i-1][j],
                      (*m_heat_grid)[i][j-1],                (*m_heat_grid)[i][j],               (*m_heat_grid)[i][j+1],
                                                                  (*m_heat_grid)[i+1][j]
                                              );

                delta = std::max(delta, std::abs((*m_previous_value_grid)[i][j] - (*m_value_grid)[i][j]));

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


    for (int i {}; i < m_value_grid->size(); ++i) {
        for (int j {}; j < (*m_value_grid)[i].size(); ++j) {
            output << (*m_value_grid)[i][j] << " ";
        }
        output << "\n";
    }

    output.close();
}

void PoissonEquationSolver::check_deltas() const {
    int k {};
    for (int time {1}; time < m_Nt; ++time) {
        if (!(m_deltas[time] <= m_deltas[time - 1])) {
            // std::cerr << std::format("iter {}: delta >= previous delta ({} >= {}).", time, m_deltas[time], m_deltas[time - 1]) << std::endl;
            ++k;
        }
    }
    std::cerr << std::format("delta >= previous delta on {}/{} iterations", k, m_Nt - 1) << std::endl;


}