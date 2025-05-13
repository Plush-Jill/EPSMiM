#include "window_left_part.hpp"


bool WindowLeftPart::is_index_covered(const long index) const {
    if (is_index_in_array(index) && is_index_in_window(index)) {
        if ((*m_control_time_array)[index] >= (m_reset_count + 1) * m_window_size - 1) {
            return false;
        }
        return true;
    }
    return false;
}

