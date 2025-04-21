//
// Created by plush-jill on 4/20/25.
//

#include "front_abstract.hpp"

#include <format>
#include <set>
#include <fstream>
void FrontAbstract::move_front_to_right() {
    ++m_front_left_edge_position;
}

void FrontAbstract::calc_front_cover_positions() {
    // std::cout << std::format(
    //     "front ({}, {}) calc positions [{} : {}] in sector {}",
    //     m_control_time_array_borders.first,
    //     m_control_time_array_borders.second,
    //     m_front_left_edge_position,
    //     m_front_left_edge_position + m_front_size - 1,
    //     get_reset_count()
    //     ) << std::endl;
    for (long i {m_front_left_edge_position}; i < m_front_left_edge_position + m_front_size; ++i) {
        if (is_index_covered(i) /*&& (*m_control_time_array)[i] <= m_total_time*/ /*&& !is_line_need_calc_in_current_reset(i)*/) {
            m_functions[(*m_control_time_array)[i] % 2](i);
            ++(*m_control_time_array)[i];
        }
    }
    // sleep(1);
}

bool FrontAbstract::is_index_ready(const long i) const {
    if ((*m_control_time_array)[i] == (*m_control_time_array)[i - 1] &&
        (*m_control_time_array)[i] == (*m_control_time_array)[i + 1]) {
        return true;
        }
    return false;
}

bool FrontAbstract::is_index_in_array(const long i) const {
    return (i > 0 && i < m_array_size - 1);

}

bool FrontAbstract::is_vertical_part_finalized() const {
    const long left {m_control_time_array_borders.first == 0 ? 1 : m_control_time_array_borders.first};
    const long right {m_control_time_array_borders.second == (m_array_size - 1) ? (m_array_size - 1) : m_control_time_array_borders.second};
    for (long i {left}; i < right; ++i) {
        if ((*m_control_time_array)[i] != m_total_time) {
            return false;
        }
    }

    return true;
}

bool FrontAbstract::is_array_finalized() const {
    for (int i {1}; i < m_array_size - 1; ++i) {
        if (!is_line_need_calc_in_current_reset(i)) {
            return false;
        }
    }

    return true;
}

bool FrontAbstract::is_sector_finalized() const {
    // const long left {m_control_time_array_borders.first == 0 ? 1 : m_control_time_array_borders.first};
    // const long right {m_control_time_array_borders.second == (m_array_size - 1) ? (m_array_size - 1) : m_control_time_array_borders.second};
    long left {m_control_time_array_borders.first};
    long right {m_control_time_array_borders.second};

    if (left == 0) {
        ++left;
    }
    if (right == m_array_size) {
        --right;
    }

    // if (m_control_time_array_borders.first != 0) {
        // std::cout << std::format(
        //     "FrontAbstract::is_sector_finalized(): (left, right) = ({}, {})",
        //     left, right
        //     ) << std::endl;
    // }
    for (long i {left}; i < right; ++i) {
        if (is_line_need_calc_in_current_reset(i)) {
            // if (m_control_time_array_borders.first != 0) {
                // std::cout << std::format(
                //     "reset: {}, time[{}] = {}",
                //     m_reset_count,
                //     i,
                //     (*m_control_time_array)[i]
                //     ) << std::endl;
            // }
            return false;
        }
    }

    return true;
}

bool FrontAbstract::is_line_need_calc_in_current_reset(const long index) const {
    if ((*m_control_time_array)[index] <= (m_reset_count + 1) * m_front_size) {
        return false;
    }
    return true;
}

bool FrontAbstract::is_front_gone() const {
    return m_front_left_edge_position >= m_control_time_array_borders.second;
}

bool FrontAbstract::is_ready_to_reset() const {
    return m_ready_to_reset;
}

void FrontAbstract::move_along() {
    // while (!is_array_finalized()) {
    while (!is_front_gone()) {
        move_front_to_right();
        calc_front_cover_positions();
    }
}
// }

void FrontAbstract::move_all_times() {
    for (int i {}; i < m_total_time / m_front_size + 1; ++i) {
        // while (!is_vertical_part_finalized()) {
        while (!is_front_gone()) {
            calc_front_cover_positions();
            move_front_to_right();
        }
        // std::cout << std::format(
        //     "front ({}, {}) finished sector {} and start waiting, sector_finalized = {}",
        //     m_control_time_array_borders.first,
        //     m_control_time_array_borders.second,
        //     get_reset_count(),
        //     is_sector_finalized()) << std::endl;
        // print_time_array_to_file(std::format("log_{}-{}-after-sector-{}.txt", m_control_time_array_borders.first, m_control_time_array_borders.second, get_reset_count()));
        wait_neighbours();
        reset();
    }
}

void FrontAbstract::print_time_array_to_file(const std::string& file_name) const {
    size_t i = 0;

    std::ofstream out(file_name, std::ios::out);
    if (!out) {
        std::cerr << std::format("can't open file {} to write", file_name) << std::endl;
        return;
    }

    while (i < m_control_time_array->size()) {
        constexpr size_t chunk_size = 50;
        size_t end = std::min(i + chunk_size, m_control_time_array->size());
        out << std::format("[{} : {}] ", i, end);
        for (size_t j = i; j < end; ++j) {
            out << (*m_control_time_array)[j];
            if (j + 1 < end) {
                out << ", ";
            }
        }
        out << "\n";
        i = end;
    }

    out.close();
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

bool FrontAbstract::is_index_in_borders(const long i) const {
    return i >= m_control_time_array_borders.first && i < m_control_time_array_borders.second;
}

long FrontAbstract::get_reset_count() const {
    return m_reset_count;
}

void FrontAbstract::set_right_neighbour(const std::shared_ptr<FrontAbstract> &neighbour) {
    m_neighbours.second = neighbour;
}

void FrontAbstract::set_left_neighbour(const std::shared_ptr<FrontAbstract> &neighbour) {
    m_neighbours.first = neighbour;
}


void FrontAbstract::wait_neighbours() {
    m_barrier->arrive_and_wait();
    // print_time_set();
    // m_barrier->arrive_and_wait();
}

void FrontAbstract::print_time_set() const {
    std::set<int> time_set;

    for (int i {}; i < m_array_size - 1; ++i) {
        if (!time_set.contains((*m_control_time_array)[i])) {
            time_set.insert((*m_control_time_array)[i]);
        }
    }
    std::cout << std::format("time_set:") << std::endl;
    for (int i : time_set) {
        std::cout << std::format("{}, ", i);
    }
    std::cout << std::endl;
}
