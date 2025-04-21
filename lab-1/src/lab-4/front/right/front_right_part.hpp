//
// Created by plush-jill on 4/20/25.
//

#ifndef FRONT_RIGHT_PART_HPP
#define FRONT_RIGHT_PART_HPP
#include <format>

#include "../front_abstract.hpp"


class FrontRightPart final : public FrontAbstract {
private:
    [[nodiscard]] bool is_index_covered(long index) const override;
    // void wait_neighbours() override;


public:
    FrontRightPart(
    const long array_size,
    const long total_time,
    const long front_size,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_a,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_b,
    const std::shared_ptr<std::vector<int>>& control_time_array,
    const std::pair<int, int>& control_time_array_borders,
    const std::shared_ptr<std::barrier<>>& barrier,
    const std::function<void(
        int,
        std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
        std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
    )>& function
    ) :
    FrontAbstract(
        array_size,
        total_time,
        front_size,
        value_grid_a,
        value_grid_b,
        control_time_array,
        control_time_array_borders,
        barrier,
        function
        ) {
        // ++m_front_left_edge_position;
        // std::cout << std::format("created right front at pos {}", m_front_left_edge_position) << std::endl;

    }
};



#endif //FRONT_RIGHT_PART_HPP
