//
// Created by plush-jill on 4/20/25.
//

#include "window_right_part.hpp"

bool WindowRightPart::is_index_covered(const long index) const {
    if (is_index_in_window(index) &&
        index < m_control_time_array_borders.second - 1) {
        if (!is_window_further_left_border() ||
            index < m_control_time_array_borders.second) {
            return true;
        }
    }
    return false;
}

// void WindowRightPart::wait_neighbours() {
    // m_barrier->arrive_and_wait();
    // m_ready_to_reset = true;
    // while (
    //     !m_neighbours.first->is_ready_to_reset()
    //     ) {}
    // m_ready_to_reset = false;
// }
