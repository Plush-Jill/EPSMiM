//
// Created by plush-jill on 4/19/25.
//

#ifndef MOVING_FRONT_ALONG_ARRAY_HPP
#define MOVING_FRONT_ALONG_ARRAY_HPP
#include <format>
#include <functional>
#include <memory>
#include <bits/ranges_algo.h>
#include "../common/aligned_allocator.hpp"



class FrontMovingAlongArray {
private:
    const long m_array_size;
    const long m_total_time;
    const long m_front_size;
    std::vector<bool> m_front;
    std::vector<int> m_control_time_array;


    int m_front_left_edge_position;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_value_grid_a;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_value_grid_b;
    std::unordered_map<int, std::function<void(
        int
        // std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
        // std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
        )>> m_functions;

    void move_front_to_right() {
        ++m_front_left_edge_position;
    }

    void calc_front_cover_positions() {
        std::cout << "calc..." << std::endl;
        for (int i {m_front_left_edge_position}; i < m_front_left_edge_position + m_front_size; ++i) {
            if (is_index_covered(i) && !is_ended(i)) {
                std::cout << std::format("calc function call will be on i = {}", i) << std::endl;
                m_functions[i % 2](i);
                ++m_control_time_array[i];
            }
            std::cout << std::format("calc function call skipped on i = {}", i) << std::endl;
        }
    }

    [[nodiscard]] bool is_all_ended () const {
        return std::ranges::all_of(m_control_time_array, [this](const int i) { return is_ended(i); });
    }
    [[nodiscard]] bool is_ended (const int index) const {
        if (m_control_time_array[index] < m_total_time) {
            return false;
        }
        return true;
    }

    [[nodiscard]] bool is_index_covered(const int index) const {
        if (index >= 0 && index <= m_array_size && (index - m_front_left_edge_position) < m_front_size && (index - m_front_left_edge_position) > 0) {
            return true;
        }
        return false;
    }
    [[nodiscard]] bool is_front_gone() const {
        return m_front_left_edge_position >= m_array_size;
    }

public:
    FrontMovingAlongArray(
        const long array_size,
        const long total_time,
        const long front_size,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_a,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_b,
        const std::function<void(
            int,
            std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
            std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
        )>& function
        ) :
    m_array_size(array_size),
    m_total_time(total_time),
    m_front_size(front_size),
    m_front_left_edge_position (-front_size) {
        m_front = std::vector<bool>(m_front_size, false);
        m_control_time_array = std::vector<int>(m_array_size, 0);

        m_value_grid_a = value_grid_a;
        m_value_grid_b = value_grid_b;

        auto place_held_function_1 = std::bind(function, std::placeholders::_1, m_value_grid_a, m_value_grid_b);
        auto place_held_function_2 = std::bind(function, std::placeholders::_1, m_value_grid_b, m_value_grid_a);


        m_functions[0] = place_held_function_1;
        m_functions[1] = place_held_function_2;

        std::cout << "front created" << std::endl;
    }

    void move_along() {
        // while (!is_all_ended()) {
            while (!is_front_gone()) {
                move_front_to_right();
                calc_front_cover_positions();
            }
        }
    // }

    void move_all_times () {
        std::cout << "start moving front" << std::endl;
        while (!is_all_ended()) {
            while (!is_front_gone()) {
                move_front_to_right();
                calc_front_cover_positions();
            }
            reset();
        }
    }

    void reset() {
        m_front_left_edge_position = -m_front_size;
    }

};



#endif //MOVING_FRONT_ALONG_ARRAY_HPP
