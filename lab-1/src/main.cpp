#include <iostream>
#include <boost/json.hpp>
#include <cmath>
#include <fstream>
#include <sstream>



class Circle {
private:
    const float m_center_x;
    const float m_center_y;
    const float m_radius;


public:
    Circle(const float center_x, const float center_y, const float radius) :
    m_center_x(center_x), m_center_y(center_y), m_radius(radius) {}

    [[nodiscard]] bool is_has_point(const double x, const double y) const {
        return (std::pow(x - m_center_x, 2) + std::pow(y - m_center_y, 2)) < std::pow(m_radius, 2);
    }
};
class HeatSourceCircle : public Circle {
private:
    const float m_heat;

public:
    HeatSourceCircle(const float center_x, const float center_y, const float radius, const float heat) :
    Circle(center_x, center_y, radius), m_heat(heat) {}
    [[nodiscard]] float get_heat() const {
        return m_heat;
    }
};

int main() {
    constexpr float Xa {0.0};
    constexpr float Xb {4.0};
    constexpr float Ya {0.0};
    constexpr float Yb {4.0};
    std::ifstream config("config.json");
    std::stringstream json_buffer;
    json_buffer << config.rdbuf();
    boost::json::object json = boost::json::parse(json_buffer.str()).as_object();
    config.close();

    long Nx {json["Nx"].as_int64()};
    long Ny {json["Ny"].as_int64()};
    long Nt {json["Nt"].as_int64()};


    float hx = (Xb - Xa) / (Nx - 1);
    float hy = (Yb - Ya) / (Ny - 1);

    float Xs1 = Xa + (Xb - Xa) / 3;
    float Ys1 = Ya + (Yb - Ya) * 2 / 3;
    float Xs2 = Xa + (Xb - Xa) * 2 / 3;
    float Ys2 = Ya + (Yb - Ya) / 3;
    float R = static_cast<float>(0.1) * std::min(Xb - Xa, Yb - Ya);
    HeatSourceCircle circle_plus_01 (Xs1, Ys1, R, 0.1);
    HeatSourceCircle circle_minus_01 (Xs2, Ys2, R, -0.1);

    std::vector<std::vector<std::pair<float, float>>> grid(Nx, std::vector<std::pair<float, float>>(Ny, {0.0, 0.0}));

    for (int x {}; x < Nx; ++x) {
        for (int y {}; y < Ny; ++y) {
            if (circle_plus_01.is_has_point(x, y)) {
                grid[x][y].first = circle_plus_01.get_heat();
            } else if (circle_minus_01.is_has_point(x, y)) {
                grid[x][y].first = circle_minus_01.get_heat();
            } else {
                grid[x][y].first = 0.0;
            }
        }
    }


    for (int time {0}; time < Nt; ++time) {
        for (int x {}; x < Nx; ++x) {
            for (int y {}; y < Ny; ++y) {

            }
        }
    }


    return 0;
}