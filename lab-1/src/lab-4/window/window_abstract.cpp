//
// Created by plush-jill on 4/20/25.
//

#include "window_abstract.hpp"

#include <cstring>
#include <format>
#include <set>
#include <fstream>
void WindowAbstract::move_window_to_right() {
    ++m_left_edge_position;
}

void WindowAbstract::calc_window_cover_positions() {
    // std::cout << std::format(
    //     "window ({}, {}) calc positions [{} : {}] in sector {}",
    //     m_control_time_array_borders.first,
    //     m_control_time_array_borders.second,
    //     m_left_edge_position,
    //     m_left_edge_position + m_window_size - 1,
    //     get_reset_count()
    //     ) << std::endl;
    for (long i {m_left_edge_position + m_window_size - 1}; i >= m_left_edge_position; --i) {
        if (is_index_covered(i) /*&& (*m_control_time_array)[i] <= m_total_time*/ /*&& !is_line_need_calc_in_current_reset(i)*/) {
            m_functions[(*m_control_time_array)[i] % 2](i);
            ++(*m_control_time_array)[i];
        }
    }
    // sleep(1);
}

bool WindowAbstract::is_index_ready(const long i) const {
    if ((*m_control_time_array)[i] == (*m_control_time_array)[i - 1] &&
        (*m_control_time_array)[i] == (*m_control_time_array)[i + 1]) {
        return true;
        }
    return false;
}

bool WindowAbstract::is_index_in_array(const long i) const {
    return (i > 0 && i < m_array_size - 1);

}

bool WindowAbstract::is_vertical_part_finalized() const {
    const long left {m_control_time_array_borders.first == 0 ? 1 : m_control_time_array_borders.first};
    const long right {m_control_time_array_borders.second == (m_array_size - 1) ? (m_array_size - 1) : m_control_time_array_borders.second};
    for (long i {left}; i < right; ++i) {
        if ((*m_control_time_array)[i] != m_total_time) {
            return false;
        }
    }

    return true;
}

bool WindowAbstract::is_array_finalized() const {
    for (int i {1}; i < m_array_size - 1; ++i) {
        if (!is_line_need_calc_in_current_reset(i)) {
            return false;
        }
    }

    return true;
}

bool WindowAbstract::is_sector_finalized() const {
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
        //     "WindowAbstract::is_sector_finalized(): (left, right) = ({}, {})",
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

bool WindowAbstract::is_line_need_calc_in_current_reset(const long index) const {
    if ((*m_control_time_array)[index] <= (m_reset_count + 1) * m_window_size) {
        return false;
    }
    return true;
}

bool WindowAbstract::is_window_gone() const {
    return m_left_edge_position >= m_control_time_array_borders.second;
}

bool WindowAbstract::is_ready_to_reset() const {
    return m_ready_to_reset;
}

void WindowAbstract::move_along() {
    // while (!is_array_finalized()) {
    while (!is_window_gone()) {
        move_window_to_right();
        calc_window_cover_positions();
    }
}
// }

void WindowAbstract::move_all_times() {
    for (int i {}; i < m_total_time / m_window_size + 1; ++i) {
        // while (!is_vertical_part_finalized()) {
        move_along();
        // print_time_array_to_file(std::format("log_{}-{}-after-sector-{}.txt", m_control_time_array_borders.first, m_control_time_array_borders.second, get_reset_count()));
        wait_neighbours();
        reset();
    }
}

void WindowAbstract::start() {
    //TODO: сделать привязку к ядру.

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(m_core_number, &cpuset);

    const pthread_t current_thread = pthread_self();
    if (const int result = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset); result != 0) {
        std::cerr << std::format("Failed to set thread affinity to core {}: {}", m_core_number, strerror(result)) << std::endl;
        throw std::runtime_error("Thread affinity setting failed");
    }


    move_all_times();
}

void WindowAbstract::print_time_array_to_file(const std::string& file_name) const {
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

void WindowAbstract::reset() {
    m_left_edge_position = m_control_time_array_borders.first - m_window_size;
    ++m_reset_count;
}

bool WindowAbstract::is_window_further_left_border() const {
    return
        (m_left_edge_position > m_control_time_array_borders.first) &&
        (m_left_edge_position < m_control_time_array_borders.second - m_window_size);
}

bool WindowAbstract::is_index_in_window(const long index) const {
    return
        index - m_left_edge_position < m_window_size &&
        index - m_left_edge_position > 0;
}

bool WindowAbstract::is_index_in_borders(const long i) const {
    return i >= m_control_time_array_borders.first && i < m_control_time_array_borders.second;
}

long WindowAbstract::get_reset_count() const {
    return m_reset_count;
}

void WindowAbstract::set_right_neighbour(const std::shared_ptr<WindowAbstract> &neighbour) {
    m_neighbours.second = neighbour;
}

void WindowAbstract::set_left_neighbour(const std::shared_ptr<WindowAbstract> &neighbour) {
    m_neighbours.first = neighbour;
}


void WindowAbstract::wait_neighbours() {
    m_barrier->arrive_and_wait();
    // print_time_set();
    // m_barrier->arrive_and_wait();
}

void WindowAbstract::print_time_set() const {
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
