#include "window_middle_part.hpp"


bool WindowMiddlePart::is_index_covered(const long index) const {
    if (is_index_in_array(index) && is_index_in_window(index)) {
        if ((*m_control_time_array)[index] >= (m_reset_count + 1) * m_window_size - 1) {
            return false;
        }
        return true;
    }
    return false;
}

















bool WindowMiddlePart::is_index_near_right_border(const long index) const {
    return (index > m_control_time_array_borders.second - m_window_size);
}

