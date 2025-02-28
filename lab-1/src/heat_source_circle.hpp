//
// Created by plush-jill on 2/12/25.
//

#ifndef HEAT_SOURCE_CIRCLE_HPP
#define HEAT_SOURCE_CIRCLE_HPP


class HeatSourceCircle final {
private:
    const float m_center_x;
    const float m_center_y;
    const float m_radius;
    const float m_heat;

public:
    HeatSourceCircle(const float center_x, const float center_y, const float radius, const float heat) :
    m_center_x(center_x), m_center_y(center_y), m_radius(radius), m_heat(heat) {
        // float heat = get_heat();
    }
    [[nodiscard]] float get_heat() const { return m_heat; }
    [[nodiscard]] bool is_has_point(float x, float y) const;

};



#endif //HEAT_SOURCE_CIRCLE_HPP
