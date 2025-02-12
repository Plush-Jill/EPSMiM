//
// Created by plush-jill on 2/12/25.
//

#ifndef HEAT_SOURCE_CIRCLE_HPP
#define HEAT_SOURCE_CIRCLE_HPP
#include "heat_source.hpp"


class HeatSourceCircle final : public HeatSource {
private:
    const float m_center_x;
    const float m_center_y;
    const float m_radius;
    const float m_heat;

public:
    HeatSourceCircle(const float center_x, const float center_y, const float radius, const float heat) :
    m_center_x(center_x), m_center_y(center_y), m_radius(radius), m_heat(heat) {}

    [[nodiscard]] bool is_has_point(const int x, const int y) const override;

};



#endif //HEAT_SOURCE_CIRCLE_HPP
