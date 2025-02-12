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

    std::vector<HeatSource> m_heat_sources;
    std::vector<std::vector<std::pair<float, float>>> m_grid;

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

};

int main() {

    PoissonEquationSolver solver {"config.json"};



    return 0;
}