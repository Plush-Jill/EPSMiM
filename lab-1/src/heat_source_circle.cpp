//
// Created by plush-jill on 2/12/25.
//

#include "heat_source_circle.hpp"
#include <cmath>


bool HeatSourceCircle::is_has_point(const float x, const float y) const {
    return (std::pow(x - m_center_x, 2) + std::pow(y - m_center_y, 2)) < std::pow(m_radius, 2);
}
