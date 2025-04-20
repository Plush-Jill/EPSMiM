//
// Created by plush-jill on 4/20/25.
//

#include "front_middle_part.hpp"

bool FrontMiddlePart::is_index_covered(const long index) const {
    if (is_index_in_front(index)) {
        if (!is_front_further_left_border()) {
            return true;
        }
        if (m_neighbours.second->is_front_further_left_border() &&
            index < m_control_time_array_borders.second &&
            (*m_control_time_array)[index] <= (m_reset_count + 1) * m_front_size) {
            return true;
        }
    }
    return false;
    //     if (is_front_further_left_border()) {
    //         if (m_neighbours.second->is_front_further_left_border() &&
    //             (*m_control_time_array)[index] <= (m_reset_count + 1) * m_front_size) {}
    //     }
    // }
    // if (index - m_front_left_edge_position < m_front_size &&
    //     index - m_front_left_edge_position > 0 &&
    //     (*m_control_time_array)[index] <= (m_reset_count + 1) * m_front_size) {
    //     return true;
    // }
    return false;
}
