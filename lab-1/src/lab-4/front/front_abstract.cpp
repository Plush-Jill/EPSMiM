//
// Created by plush-jill on 4/20/25.
//

#include "front_abstract.hpp"

void FrontAbstract::move_front_to_right() {
    ++m_front_left_edge_position;
}

void FrontAbstract::calc_front_cover_positions() {
    for (long i {m_front_left_edge_position}; i < m_front_left_edge_position + m_front_size; ++i) {
        if (is_index_covered(i) && !is_line_finalized(i)) {
            while (!is_index_ready(i)) {}
            m_functions[(*m_control_time_array)[i] % 2](i);
            ++(*m_control_time_array)[i];
        }
    }
}

bool FrontAbstract::is_index_ready(const long i) const {
    if ((*m_control_time_array)[i] == (*m_control_time_array)[i - 1] &&
        (*m_control_time_array)[i] == (*m_control_time_array)[i + 1]) {
        return true;
        }
    return false;
}

bool FrontAbstract::is_all_ended() const {
    for (int i {1}; i < m_array_size - 1; ++i) {
        if (!is_line_finalized(i)) {
            return false;
        }
    }

    return true;
}

bool FrontAbstract::is_line_finalized(const long index) const {
    if ((*m_control_time_array)[index] < m_total_time) {
        return false;
    }
    return true;
}

bool FrontAbstract::is_front_gone() const {
    return m_front_left_edge_position >= m_control_time_array_borders.second;
}

void FrontAbstract::move_along() {
    // while (!is_all_ended()) {
    while (!is_front_gone()) {
        move_front_to_right();
        calc_front_cover_positions();
    }
}
// }

void FrontAbstract::move_all_times() {
    while (!is_all_ended()) {
        while (!is_front_gone()) {
            move_front_to_right();
            calc_front_cover_positions();
        }
        reset();
    }
}

void FrontAbstract::reset() {
    m_front_left_edge_position = m_control_time_array_borders.first - m_front_size;
    ++m_reset_count;
}

bool FrontAbstract::is_front_further_left_border() const {
    return
        (m_front_left_edge_position > m_control_time_array_borders.first) &&
        (m_front_left_edge_position < m_control_time_array_borders.second - m_front_size);
}

bool FrontAbstract::is_index_in_front(const long index) const {
    return
        index - m_front_left_edge_position < m_front_size &&
        index - m_front_left_edge_position > 0;
}
