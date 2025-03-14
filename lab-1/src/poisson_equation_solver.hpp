//
// Created by plush-jill on 3/4/25.
//

#ifndef POISSON_EQUATION_SOLVER_HPP
#define POISSON_EQUATION_SOLVER_HPP
#include <iostream>
#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "heat_source_circle.hpp"
#include <format>


class PoissonEquationSolver {
private:
    const float m_Xa;
    const float m_Xb;
    const float m_Ya;
    const float m_Yb;

    long m_Nx;
    long m_Ny;
    long m_Nt;

    float m_hx;
    float m_hy;

    float hx_pow_minus2;
    float hy_pow_minus2;
    float m_hx_hy_2;
    float m_a;
    float m_b;
    float m_c;


    std::vector<HeatSourceCircle> m_heat_sources;
    std::shared_ptr<std::vector<std::vector<std::pair<float, float>>>> m_grid; //first - heat, second - function value
    std::shared_ptr<std::vector<std::vector<std::pair<float, float>>>> m_previous_grid; //first - heat, second - function value

    std::shared_ptr<std::vector<std::vector<float>>> m_value_grid;
    std::shared_ptr<std::vector<std::vector<float>>> m_previous_value_grid;
    std::shared_ptr<std::vector<std::vector<float>>> m_heat_grid;

    std::vector<float> m_deltas;

public:
    explicit PoissonEquationSolver(const std::string& config_file);

    [[nodiscard]] float calc_new_value (float F_im1_jm1, float F_im1_j, float F_im1_jp1,
                                        float F_i_jm1,                  float F_i_jp1,
                                        float F_ip1_jm1, float F_ip1_j, float F_ip1_jp1,

                                                         float P_im1_j,
                                        float P_i_jm1,   float P_i_j,   float P_i_jp1,
                                                         float P_ip1_j
            ) const;

    void solve();

    void solve2();

    void export_grid_value_as_matrix(const std::string& file_path) const;

    void export_grid_value_as_matrix2(const std::string& file_path) const;

    void check_deltas() const;
};



#endif //POISSON_EQUATION_SOLVER_HPP
