#include <iostream>
#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "heat_source_circle.hpp"


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

    float m_hx_hy_2;
    float m_a;
    float m_b;
    float m_c;


    std::vector<HeatSource> m_heat_sources;
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
        m_heat_sources.push_back(HeatSourceCircle {Xs1, Ys1, R, 0.1});
        m_heat_sources.push_back(HeatSourceCircle {Xs2, Ys2, R, -0.1});
        // HeatSourceCircle circle_plus_01 (Xs1, Ys1, R, 0.1);
        // HeatSourceCircle circle_minus_01 (Xs2, Ys2, R, -0.1);

        m_grid = std::vector<std::vector<std::pair<float, float>>> (m_Nx, std::vector<std::pair<float, float>>(m_Ny, {0.0, 0.0}));

        for (int x {}; x < m_Nx; ++x) {
            for (int y {}; y < m_Ny; ++y) {
                m_grid[x][y].first = 0.0;
                for (auto& heat_source : m_heat_sources) {
                    if (heat_source.is_has_point(x, y)) {
                        m_grid[x][y].first = heat_source.get_heat();
                    }
                }
            }
        }
    }



    void solve() {
        m_hx_hy_2 = static_cast<float>(std::pow(m_hx, -2) + std::pow(m_hy, -2));
        m_a = 0.2f * m_hx_hy_2;
        m_b = 0.5f * static_cast<float>(5 * std::pow(m_hx, -2) - std::pow(m_hy, -2));
        m_c = 0.25f * m_hx_hy_2;

        auto function = [this](float F_im1_jm1, float F_im1_j, float F_im1_jp1,
                               float F_i_jm1,                  float F_i_jp1,
                               float F_ip1_jm1, float F_ip1_j, float F_ip1_jp1,

                                                float P_im1_j,
                               float P_i_jm1,   float P_i_j,   float P_i_jp1,
                                                float P_ip1_j
        ) {
            const float result = static_cast<float> (m_a * (m_b * (F_i_jm1 + F_i_jp1) + m_b * (F_im1_j + F_im1_j + F_ip1_j)
                + m_c * (F_im1_jp1 + F_im1_jp1 + F_ip1_jm1 + F_ip1_jp1)
                + 2 * P_i_j + 0.25 * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1) ));

            return result;
        };

        // F_i_jp1, F_i_jm1, F_ip1_jp1, F_ip1_j, F_ip1_jm1, F_im1_jp1, F_im1_j, F_im1_jm1
        for (int time {}; time < m_Nt; ++time) {
            for (int i {1}; i < m_Nx - 1; ++i) {
                for (int j {1}; j < m_Ny - 1; ++j) {
                    m_grid[i][j].first = function(m_grid[i-1][j-1].second, m_grid[i-1][j].second, m_grid[i-1][j+1].second,
                                                  m_grid[i][j-1].second,                                     m_grid[i][j+1].second,
                                                  m_grid[i+1][j-1].second,   m_grid[i+1][j].second, m_grid[i+1][j+1].second,

                                                                             m_grid[i-1][j].second,
                                                  m_grid[i][j-1].second, m_grid[i][j].second, m_grid[i][j+1].second,
                                                                              m_grid[i+1][j].second
                                                  );
                }
            }
        }
    }
};

int main() {

    PoissonEquationSolver solver {"config.json"};



    return 0;
}