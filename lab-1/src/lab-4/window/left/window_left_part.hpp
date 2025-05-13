//
// Created by plush-jill on 4/19/25.
//

#ifndef MOVING_FRONT_ALONG_ARRAY_HPP
#define MOVING_FRONT_ALONG_ARRAY_HPP
#include <format>

#include "../window_abstract.hpp"



class WindowLeftPart final : public WindowAbstract {
private:
    [[nodiscard]] bool is_index_covered(long index) const override;
    // void wait_neighbours() override;

public:
    WindowLeftPart(
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
        // std::cout << std::format("created left window at pos {}", m_left_edge_position) << std::endl;

    }

};



#endif //MOVING_FRONT_ALONG_ARRAY_HPP
