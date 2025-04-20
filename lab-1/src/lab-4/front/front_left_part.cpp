//
// Created by plush-jill on 4/19/25.
//

#include "front_left_part.hpp"

bool FrontLeftPart::is_index_covered(const long index) const {
    if (index > 0 &&
        index < m_array_size - 1 &&
        (index - m_front_left_edge_position) < m_front_size &&
        (index - m_front_left_edge_position) > 0) {
        return true;
        }
    return false;
}
