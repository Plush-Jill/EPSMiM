//
// Created by plush-jill on 4/20/25.
//

#ifndef FRONT_RIGHT_PART_HPP
#define FRONT_RIGHT_PART_HPP
#include <format>

#include "../window_abstract.hpp"


class WindowRightPart final : public WindowAbstract {
private:
    [[nodiscard]] bool is_index_covered(long index) const override;
    // void wait_neighbours() override;


public:
    WindowRightPart(
    const long array_size,
    const long total_time,
    const long front_size,
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
    WindowAbstract(
        array_size,
        total_time,
        front_size,
        value_grid_a,
        value_grid_b,
        control_time_array,
        control_time_array_borders,
        barrier,
        core_number,
        function
        ) {
        // ++m_left_edge_position;
        // std::cout << std::format("created right window at pos {}", m_left_edge_position) << std::endl;

    }
};



#endif //FRONT_RIGHT_PART_HPP
