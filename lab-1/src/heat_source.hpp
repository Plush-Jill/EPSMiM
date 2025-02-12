//
// Created by plush-jill on 2/12/25.
//

#ifndef HEAT_SOURCE_HPP
#define HEAT_SOURCE_HPP



class HeatSource {
    const float m_heat {};
public:
    virtual ~HeatSource() = default;
    [[nodiscard]] virtual bool is_has_point(int x, int y) const = 0;
    [[nodiscard]] virtual float get_heat() const { return m_heat; }
};



#endif //HEAT_SOURCE_HPP
