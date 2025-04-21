#include "front_middle_part.hpp"


bool FrontMiddlePart::is_index_covered(const long index) const {
    if (is_index_in_array(index) && is_index_in_front(index)) {
        if ((*m_control_time_array)[index] >= (m_reset_count + 1) * m_front_size - 1) {
            return false;
        }
        return true;
    }
    return false;
}

















bool FrontMiddlePart::is_index_near_right_border(const long index) const {
    return (index > m_control_time_array_borders.second - m_front_size);
}

