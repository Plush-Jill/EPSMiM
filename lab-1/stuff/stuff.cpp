#include <climits>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <format>
#include <functional>
#include <iostream>
#include <memory>
#include <bits/ranges_algo.h>

#include "../src/common/aligned_allocator.hpp"


using Grid = std::vector<std::vector<float>>;

void solve_with_wide_front(
) {
    int Nt = 100;
    int Nx = 40;
    int Ny = 40;
    int front_height = 10;

    // std::vector<Grid*> buffers = {&A, &B};
    std::string A = "A";
    std::string B = "B";
    std::vector<std::string*> buffers = {&A, &B};

    // Сдвиг времени - отложенный по фронту
    for (int time = 0; time < Nt + front_height - 1; ++time) {
        float delta = INT_MIN;

        for (int local_step = 0; local_step < front_height; ++local_step) {
            int current_time = time - local_step;
            int row = local_step;

            if (current_time < 0 || current_time >= Nt) continue;
            if (row <= 0 || row >= Ny - 1) continue;

            // Определим, какие строки мы можем считать на этом уровне фронта
            int base_row = time - front_height + 1;
            int i = base_row + local_step;

            if (i <= 0 || i >= Ny - 1) continue;

            // Определяем текущий буфер для этого временного шага
            std::string* from = buffers[(current_time + 1) % 2];
            std::string* to   = buffers[current_time % 2];
            std::cout << std::format("calc string {} for time {}, from matrix {} to matrix {}\n", i, time, *from, *to);

            // horizontal_step(delta, i, *from, *to);
        }
    }
}

class MovingFrontUnderArray {
private:
    const int m_array_size;
    const int m_total_time;
    const int m_front_size;
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

    MovingFrontUnderArray(
        const int array_size,
        const int total_time,
        const int front_size,
        const std::function<void(
                int,
                std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
                std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
                )>& function_1,
        const std::function<void(
                int,
                std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>,
                std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>
                )>& function_2,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_a,
        const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid_b
        ) :
    m_array_size(array_size),
    m_total_time(total_time),
    m_front_size(front_size),
    m_front_left_edge_position (-front_size) {
        m_front = std::vector<bool>(m_front_size, false);
        m_control_time_array = std::vector<int>(m_array_size, 0);

        m_value_grid_a = value_grid_a;
        m_value_grid_b = value_grid_b;

        auto place_held_function_1 = std::bind(function_1, std::placeholders::_1, m_value_grid_a, m_value_grid_b);
        auto place_held_function_2 = std::bind(function_2, std::placeholders::_1, m_value_grid_b, m_value_grid_a);


        m_functions[0] = place_held_function_1;
        m_functions[1] = place_held_function_2;
    }

    void move_along() {
        while (!is_all_ended()) {
            while (!is_front_gone()) {
                move_front_to_right();
                calc_front_cover_positions();
            }
        }
    }

    void reset() {
        m_front_left_edge_position = -m_front_size;
    }

    void move_front_to_right() {
        ++m_front_left_edge_position;
    }

    void calc_front_cover_positions() {
        for (int i {m_front_left_edge_position}; i < m_front_left_edge_position + m_front_size; ++i) {
            if (is_index_covered(i) && !is_ended(i)) {
                m_functions[i % 2](i);
                ++m_control_time_array[i];
            }
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


};

void front_trying() {
    int size_y = 100;
    int total_time = 100;
    int front_size = 20;
    std::vector<bool> front(front_size, false);
    std::vector<int> array (size_y, 0);



}
class A {

public:
    virtual ~A() = default;

    virtual void foo() = 0;
    virtual bool is_virtual() {
        return true;
    }
    [[nodiscard]] const std::type_info& get_type () const {
        return typeid(*this);
    }
};

class B final : public A {
public:
    B() : A() {};
    ~B() override = default;
    void foo() override {
        std::cout << "B::foo" << std::endl;
    }
};

int main() {
    const B b {};
    std::cout << b.get_type().name() << std::endl;
    std::cout << std::format("is same type: {}", b.get_type() == typeid(B));
    // solve_with_wide_front();
    return 0;

    constexpr int Nt = 100;
    constexpr int Nx = 40;
    constexpr int Ny = 40;
    constexpr int front_height = 10;
    std::string matrix_a = "A";
    std::string matrix_b = "B";
    std::vector<bool> mask(front_height, false);

    for (int time {}; time < Nt; ++time) {
        float delta = INT_MIN;
        std::swap(matrix_a, matrix_b);

        for (int i = 1 - front_height + 1; i < Ny - 1; ++i) {
            // Обновляем маску
            for (int j = 0; j < front_height; ++j) {
                int row = i + j;
                if (row <= 0 || row >= Ny - 1) {
                    mask[j] = false;
                } else if (i < front_height) {
                    mask[j] = (j <= i);  // наращиваем фронт
                } else if (i >= Ny - 1 - front_height) {
                    mask[j] = (j < Ny - 1 - i);  // сужаем фронт
                } else {
                    mask[j] = true;  // постоянный фронт
                }
            }

            // Применяем шаг по всем активным строкам
            for (int j = 0; j < front_height; ++j) {
                if (mask[j]) {
                    std::cout << std::format("calc string {} for step {}, from matrix {} to matrix {}\n", i + j, time, matrix_a, matrix_b);
                    // horizontal_step(delta, i + j, m_previous_value_grid, m_value_grid);
                }
            }
        }

        // m_deltas[time] = delta;
    }

    return 0;
}
