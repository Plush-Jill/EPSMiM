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
    std::vector<std::vector<std::pair<float, float>>> m_grid; //first - heat, second - function value
    std::vector<float> m_deltas;

public:
    explicit PoissonEquationSolver(const std::string& config_file) :
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

        m_hx = (m_Xb - m_Xa) / static_cast<float>(m_Nx - 1);
        m_hy = (m_Yb - m_Ya) / static_cast<float>(m_Ny - 1);

        float Xs1 = m_Xa + (m_Xb - m_Xa) / 3;
        float Ys1 = m_Ya + (m_Yb - m_Ya) * 2 / 3;
        float Xs2 = m_Xa + (m_Xb - m_Xa) * 2 / 3;
        float Ys2 = m_Ya + (m_Yb - m_Ya) / 3;
        float R = static_cast<float>(0.1) * std::min(m_Xb - m_Xa, m_Yb - m_Ya);

        m_heat_sources.emplace_back(Xs1, Ys1, R, 0.1);
        m_heat_sources.emplace_back(Xs2, Ys2, R, -0.1);

        m_grid = std::vector<std::vector<std::pair<float, float>>> (m_Nx, std::vector<std::pair<float, float>>(m_Ny, {0.0, 0.0}));
        m_deltas = std::vector<float> (m_Nt, 0.0);

        for (int i {}; i < m_Nx; ++i) {
            for (int j {}; j < m_Ny; ++j) {
                m_grid[i][j].first = 0.0;
                m_grid[i][j].second = 0.0;
            }
        }
        for (auto& heat_source : m_heat_sources) {
            for (int j {}; j < m_Ny; ++j) {
                for (int i {}; i < m_Nx; ++i) {

                    float x_i = m_Xa + (static_cast<float>(i) * m_hx);
                    float y_j = m_Ya + (static_cast<float>(j) * m_hy);
                    if (heat_source.is_has_point(x_i, y_j)) {
                        m_grid[j][i].first = heat_source.get_heat();
                    }
                }
            }
        }



        // m_hx_hy_2 = static_cast<float>(1.0f / std::pow(m_hx, 2) + 1.0f / std::pow(m_hy, 2));
        const float hx_pow_minus1 = 1.0f / m_hx;
        const float hy_pow_minus1 = 1.0f / m_hy;
        hx_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hx, 2));
        hy_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hy, 2));

        m_hx_hy_2 = hx_pow_minus2 + hy_pow_minus2;
        m_a = 1.0f / (5.0f * m_hx_hy_2);
        m_b = static_cast<float>((1.0f / 2.0f) * (5.0f * hy_pow_minus2 - 1.0f * hx_pow_minus2));
        m_c = m_hx_hy_2 / 4;

        std::cout << std::format(
            "m_Xa: {}\n"
            "m_Xb: {}\n"
            "m_Ya: {}\n"
            "m_Yb: {}\n"
            "m_Nx: {}\n"
            "m_Ny: {}\n"
            "m_Nt: {}\n"
            "m_hx: {}\n"
            "m_hy: {}\n"
            "hx_pow_minus2: {}\n"
            "hy_pow_minus2: {}\n"
            "m_hx_hy_2: {}\n"
            "m_a: {}\n"
            "m_b: {}\n"
            "m_c: {}",
            m_Xa,
            m_Xb,
            m_Ya,
            m_Yb,
            m_Nx,
            m_Ny,
            m_Nt,
            m_hx,
            m_hy,
            hx_pow_minus2,
            hy_pow_minus2,
            m_hx_hy_2,
            m_a,
            m_b,
            m_c
            ) << std::endl;

    }

    [[nodiscard]] float calc_new_value (float F_im1_jm1, float F_im1_j, float F_im1_jp1,
                                        float F_i_jm1,                  float F_i_jp1,
                                        float F_ip1_jm1, float F_ip1_j, float F_ip1_jp1,

                                                         float P_im1_j,
                                        float P_i_jm1,   float P_i_j,   float P_i_jp1,
                                                         float P_ip1_j
            ) const;

    void solve();

    void export_grid_value_as_matrix(const std::string& file_path) const;

    void check_deltas() const;
};



#endif //POISSON_EQUATION_SOLVER_HPP
