//
// Created by plush-jill on 4/20/25.
//

#ifndef FRONT_MOVING_ALONG_ARRAY_ABSTRACT_HPP
#define FRONT_MOVING_ALONG_ARRAY_ABSTRACT_HPP

#include <functional>
#include <memory>
#include <barrier>
#include "../../common/aligned_allocator.hpp"



class WindowAbstract {
protected:
    const long m_array_size;
    const long m_total_time;
    const long m_window_size;
    std::vector<bool> m_window;
    std::shared_ptr<std::vector<int>> m_control_time_array;
    std::pair<int, int> m_control_time_array_borders;
    long m_left_edge_position;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_value_grid_a;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_value_grid_b;
    std::unordered_map<long, std::function<void(long)>> m_functions;
    std::pair<std::shared_ptr<WindowAbstract>, std::shared_ptr<WindowAbstract>> m_neighbours;
    long m_reset_count;
    bool m_ready_to_reset;
    std::shared_ptr<std::barrier<>> m_barrier;
    int m_core_number;

    virtual void move_window_to_right();
    virtual void calc_window_cover_positions();
    [[nodiscard]] bool is_index_ready(long i) const;
    [[nodiscard]] bool is_index_in_array(long i) const;
    [[nodiscard]] bool is_index_in_borders(long i) const;
    [[nodiscard]] bool is_vertical_part_finalized() const;
    [[nodiscard]] virtual bool is_array_finalized () const;
    [[nodiscard]] virtual bool is_sector_finalized () const;
    [[nodiscard]] virtual bool is_line_need_calc_in_current_reset (long index) const;
    [[nodiscard]] virtual bool is_index_covered(long index) const = 0;
    [[nodiscard]] virtual bool is_window_gone() const;
    virtual void wait_neighbours();
    void print_time_set() const;
    void print_time_array_to_file(const std::string& file_name) const;

public:
    virtual ~WindowAbstract() = default;

    WindowAbstract(
        const long array_size,
        const long total_time,
        const long window_size,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_a,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_b,
        const std::shared_ptr<std::vector<int>>& control_time_array,
        const std::pair<int, int>& control_time_array_borders,
        const std::shared_ptr<std::barrier<>>& barrier,
        const int core_number,
        const std::function<void(
            int,
            std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
            std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
        )>& function
        ) :
    m_array_size(array_size),
    m_total_time(total_time),
    m_window_size(window_size),
    m_control_time_array_borders(control_time_array_borders),
    m_left_edge_position (m_control_time_array_borders.first - m_window_size),
    m_reset_count(0),
    m_barrier(barrier),
    m_core_number(core_number) {
        m_window = std::vector<bool>(m_window_size, false);
        m_control_time_array = control_time_array;

        m_value_grid_a = value_grid_a;
        m_value_grid_b = value_grid_b;
        m_neighbours = std::pair<std::shared_ptr<WindowAbstract>, std::shared_ptr<WindowAbstract>> (nullptr, nullptr);
        auto place_held_function_1 = std::bind(function, std::placeholders::_1, m_value_grid_a, m_value_grid_b);
        auto place_held_function_2 = std::bind(function, std::placeholders::_1, m_value_grid_b, m_value_grid_a);


        m_functions[0] = place_held_function_2;
        m_functions[1] = place_held_function_1;
        m_ready_to_reset = false;

    }

    virtual void move_along();
    virtual void move_all_times ();
    void start();
    void reset();
    [[nodiscard]] bool is_window_further_left_border() const;
    [[nodiscard]] bool is_index_in_window(long index) const;
    [[nodiscard]] virtual long get_reset_count() const;
    void set_right_neighbour(const std::shared_ptr<WindowAbstract> &neighbour);
    void set_left_neighbour(const std::shared_ptr<WindowAbstract> &neighbour);
    [[nodiscard]] virtual bool is_ready_to_reset() const;
};


#endif //FRONT_MOVING_ALONG_ARRAY_ABSTRACT_HPP
