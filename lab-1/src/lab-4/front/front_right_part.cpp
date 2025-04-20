//
// Created by plush-jill on 4/20/25.
//

#include "front_right_part.hpp"

bool FrontRightPart::is_index_covered(long index) const {
    if (is_index_in_front(index)) {
        if (!is_front_further_left_border() ||
            index < m_control_time_array_borders.second) {
            return true;
        }
    }
    return false;
}
