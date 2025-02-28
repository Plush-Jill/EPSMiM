#include <iostream>
#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "heat_source_circle.hpp"
#include <format>

class PoissonEquationSolver {
// private:
public:
    const float m_Xa;
    const float m_Xb;
    const float m_Ya;
    const float m_Yb;

    long m_Nx;
    long m_Ny;
    long m_Nt;

    float m_hx;
    float m_hy;

    float m_hx_hy_2;
    float m_a;
    float m_b;
    float m_c;


    std::vector<HeatSourceCircle> m_heat_sources;
    std::vector<std::vector<std::pair<float, float>>> m_grid; //first - heat, second - function

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
        m_heat_sources.emplace_back(Xs1, Ys1, R, 1.0f / 10);
        m_heat_sources.emplace_back(Xs2, Ys2, R, - (1.0f / 10));

        m_grid = std::vector<std::vector<std::pair<float, float>>> (m_Nx, std::vector<std::pair<float, float>>(m_Ny, {0.0, 0.0}));


        for (int i {}; i < m_Nx; ++i) {
            for (int j {}; j < m_Ny; ++j) {
                m_grid[i][j].first = 0.0;
                m_grid[i][j].second = 0.0;
            }
        }
        for (auto& heat_source : m_heat_sources) {
            // std::cout << std::format("Heat source with heat {}\n", heat_source.get_heat());
            for (int j {}; j < m_Ny; ++j) {
                for (int i {}; i < m_Nx; ++i) {

                    float x_i = m_Xa + (static_cast<float>(i) * m_hx);
                    float y_j = m_Ya + (static_cast<float>(j) * m_hy);
                    if (heat_source.is_has_point(x_i, y_j)) {
                        // std::cout << std::format("Point ({}, {}) is in heat source with heat {}", i, j, heat_source.get_heat()) << std::endl;
                        m_grid[j][i].first = heat_source.get_heat();
                    }
                }
            }
        }



        m_hx_hy_2 = static_cast<float>(std::pow(m_hx, -2) + std::pow(m_hy, -2));
        // m_a = 0.2f * m_hx_hy_2;
        m_a = static_cast<float>(1 / (5 * (1 / std::pow(m_hx, 2) + 1 / std::pow(m_hy, 2))));
        // m_b = 0.5f * static_cast<float>(5 * std::pow(m_hx, -2) - std::pow(m_hy, -2));
        m_b = static_cast<float>((1 / 2.0) * (5 / std::pow(m_hx, 2) - 1 / std::pow(m_hy, 2)));
        // m_c = 0.25f * m_hx_hy_2;
        m_c = static_cast<float>((1.0 / 4) * (1 / std::pow(m_hx, 2) + 1 / std::pow(m_hy, 2)));

    }

    inline float func (float F_im1_jm1, float F_im1_j, float F_im1_jp1,
                       float F_i_jm1,                  float F_i_jp1,
                       float F_ip1_jm1, float F_ip1_j, float F_ip1_jp1,

                                        float P_im1_j,
                       float P_i_jm1,   float P_i_j,   float P_i_jp1,
                                        float P_ip1_j
            ) {

        const float first_line = m_b * (F_i_jm1 + F_i_jp1) + m_b * (F_im1_j + F_ip1_j);
        const float second_line = m_c * (F_im1_jm1 + F_im1_jp1 + F_ip1_jm1 + F_ip1_jp1);
        const float third_line = static_cast<float>(2 * P_i_j + 0.25 * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1));
        const float result = m_a * (first_line + second_line + third_line);
        // if (F_ip1_jp1 != 0.0 || P_i_jp1 != 0.0) {
        //     std::cout << std::format ("m_hx: {}\nm_hy: {}\nm_a: {}\nm_b: {}\nm_c: {}\nLines:\nfirst: {}\nsecond: {}\nthird: {}\nresult: {}\n",
        //         m_hx, m_hy,
        //         m_a, m_b, m_c,
        //         first_line, second_line, third_line, result);
        // }
        // const float result = static_cast<float> (m_a *
        //                                         (m_b * (F_i_jm1 + F_i_jp1) + m_b * (F_im1_j + F_ip1_j)
        //                                         + m_c * (F_im1_jm1 + F_im1_jp1 + F_ip1_jm1 + F_ip1_jp1)
        //                                         + 2 * P_i_j + 0.25 * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1) ));


        // float third_line = static_cast<float>(2) * P_i_j + 0.25 * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1);
        // if (P_i_j != 0.0 || P_i_jp1 != 0.0 || P_i_jm1 != 0.0 || P_ip1_j != 0.0 || P_im1_j != 0.0) {
        //     std::cout << std::format("P_i_j = {}, P_i_jp1 = {}, P_i_jm1 = {}, P_ip1_j = {}, P_im1_j = {}\n", P_i_j, P_i_jp1, P_i_jm1, P_ip1_j, P_im1_j);
        // }
        // if (third_line != 0) {
        //     std::cout << third_line << std::endl;
        // }

        return result;
    };

    void solve() {

        std::ofstream output ("float2.dat", std::ios::out);
        if (!output.is_open()) {
            std::cerr << "Failed to open file: " << "float2.dat" << std::endl;
            throw std::exception();
        }


        for (int time {}; time < m_Nt; ++time) {
            std::vector prev_grid = m_grid;
            float delta;
            float prev_delta = delta;
            delta = INT_MIN;
            for (int i {1}; i < m_Ny - 1; ++i) {
                for (int j {1}; j < m_Nx - 1; ++j) {
                    // m_grid[i][j].second = func(prev_grid[i-1][j-1].second, prev_grid[i-1][j].second, prev_grid[i-1][j+1].second,
                    //                               prev_grid[i][j-1].second,                                     prev_grid[i][j+1].second,
                    //                               prev_grid[i+1][j-1].second,   prev_grid[i+1][j].second, prev_grid[i+1][j+1].second,
                    //
                    //                                                          prev_grid[i-1][j].first,
                    //                               prev_grid[i][j-1].first, prev_grid[i][j].first, prev_grid[i][j+1].first,
                    //                                                           prev_grid[i+1][j].first
                    //                               );
                    // m_grid[j][i].second = func(prev_grid[j-1][i-1].second, prev_grid[j-1][i].second, prev_grid[j-1][i+1].second,
                    //                               prev_grid[j][i-1].second,                                     prev_grid[j][i+1].second,
                    //                               prev_grid[j+1][i-1].second,   prev_grid[j+1][i].second, prev_grid[j+1][i+1].second,
                    //
                    //                                                          prev_grid[j-1][i].first,
                    //                               prev_grid[j][i-1].first, prev_grid[j][i].first, prev_grid[j][i+1].first,
                    //                                                           prev_grid[j+1][i].first
                    //                               );
                    m_grid[i][j].second = func(prev_grid[i+1][j-1].second, prev_grid[i+1][j].second, prev_grid[i+1][j+1].second,
                                                  prev_grid[i][j-1].second,                                     prev_grid[i][j+1].second,
                                                  prev_grid[i-1][j-1].second,   prev_grid[i-1][j].second, prev_grid[i-1][j+1].second,

                                                                             prev_grid[i+1][j].first,
                                                  prev_grid[i][j-1].first, prev_grid[i][j].first, prev_grid[i][j+1].first,
                                                                              prev_grid[i-1][j].first
                                                  );

                    // std::cout << m_grid[i][j].second << std::endl;
                    delta = std::max(delta, std::abs(m_grid[i][j].second - prev_grid[i][j].second));

                }
            }
            // output << "\n";
            if (delta > prev_delta && time > 0) {
                std::cerr << std::format("Current delta > previous delta ({} > {})\n", delta, prev_delta);
            }
        }


    }

    void print_grid_space_property(const std::string& file_path) {
        std::ofstream output (file_path, std::ios::out);
        if (!output.is_open()) {
            std::cerr << "Failed to open file: " << file_path << std::endl;
            throw std::exception();
        }



        for (int i {}; i < m_grid.size(); ++i) {
            for (int j {}; j < m_grid[i].size(); ++j) {
                output << m_grid[i][j].first << " ";
            }
            output << "\n";
        }

        output.close();
    }

    void print_grid_space_property2(const std::string& file_path) {
        std::ofstream output(file_path, std::ios::out | std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Failed to open file: " << file_path << std::endl;
            throw std::exception();
        }

        for (int i = 0; i < m_grid.size(); ++i) {
            for (int j = 0; j < m_grid[i].size(); ++j) {
                float value = m_grid[i][j].first;
                output.write(reinterpret_cast<const char*>(&value), sizeof(float));
            }
        }

        output.close();
    }

    void export_grid_space_property(const std::string& file_path) {
        std::ofstream output (file_path, std::ios::out);
        if (!output.is_open()) {
            std::cerr << "Failed to open file: " << "float2.dat" << std::endl;
            throw std::exception();
        }



        for (int i {}; i < m_grid.size(); ++i) {
            for (int j {}; j < m_grid[i].size(); ++j) {
                output << std::format("({},{}){}\n", j, i, m_grid[i][j].first);
            }
        }

        output.close();
    }


    void print_grid_value(const std::string& file_path) const {
        std::ofstream output (file_path, std::ios::out);
        if (!output.is_open()) {
            std::cerr << std::format("Failed to open file: {}", file_path) << std::endl;
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
    void export_grid_value(const std::string& file_path) {
        std::ofstream output (file_path, std::ios::out);
        if (!output.is_open()) {
            std::cerr << "Failed to open file: " << "float2.dat" << std::endl;
            throw std::exception();
        }



        // for (int i {}; i < m_grid.size(); ++i) {
        //     for (int j {}; j < m_grid[i].size(); ++j) {
        //         output << std::format("({},{}){}\n", i, j, m_grid[i][j].second);
        //     }
        // }

        for (int i {}; i < m_grid.size(); ++i) {
            for (int j {}; j < m_grid[i].size(); ++j) {
                output << m_grid[i][j].second << " ";
            }
            output << std::endl;
        }

        output.close();
    }
};

int main() {

    PoissonEquationSolver solver {"../config.json"};
    // solver.print_grid_space_property("../float2.dat");
    // solver.export_grid_space_property("../float2.dat");
    solver.solve();
    // solver.print_grid_value("../float2.dat");
    solver.export_grid_value("../float2.dat");


    return 0;
}