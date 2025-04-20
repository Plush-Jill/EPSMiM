//
// Created by plush-jill on 4/19/25.
//

#ifndef MOVING_FRONT_ALONG_ARRAY_HPP
#define MOVING_FRONT_ALONG_ARRAY_HPP
#include "front_abstract.hpp"



class FrontLeftPart final : protected FrontAbstract {
private:
    [[nodiscard]] bool is_index_covered(long index) const override;


public:
    FrontLeftPart(
        const long array_size,
        const long total_time,
        const long front_size,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_a,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_b,
        const std::shared_ptr<std::vector<int>>& control_time_array,
        const std::pair<int, int>& control_time_array_borders,
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
        function) {

    }

};



#endif //MOVING_FRONT_ALONG_ARRAY_HPP
