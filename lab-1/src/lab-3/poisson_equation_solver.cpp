//
// Created by plush-jill on 3/4/25.
//

#include "poisson_equation_solver.hpp"
#include <fstream>
#include <iostream>
#include <immintrin.h>
#include "front_moving_along_array.hpp"

PoissonEquationSolver::PoissonEquationSolver(const std::string &config_file) :
    m_Xa(0), m_Xb(4.0),
    m_Ya(0), m_Yb(4.0) {
    std::ifstream config(config_file);
    std::stringstream json_buffer;
    json_buffer << config.rdbuf();
    boost::json::object json;
    json = boost::json::parse(json_buffer.str()).as_object();
    config.close();

    m_Nx = json["Nx"].as_int64();
    m_Ny = json["Ny"].as_int64();
    m_Nt = json["Nt"].as_int64();

    m_hx = (m_Xb - m_Xa) / static_cast<float>(m_Nx - 1);
    m_hy = (m_Yb - m_Ya) / static_cast<float>(m_Ny - 1);

    float Xs1 = m_Xa + (m_Xb - m_Xa) / 3;
    float Ys1 = m_Ya + (m_Yb - m_Ya) * 2 / 3;
    float Xs2 = m_Xa + (m_Xb - m_Xa) * 2 / 3;
    float Ys2 = m_Ya + (m_Yb - m_Ya) / 3;
    float R = static_cast<float>(0.1) * std::min(m_Xb - m_Xa, m_Yb - m_Ya);

    m_heat_sources.emplace_back(Xs1, Ys1, R, 0.1);
    m_heat_sources.emplace_back(Xs2, Ys2, R, -0.1);

    m_value_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));
    m_previous_value_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));
    m_heat_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));


    m_deltas = std::vector<float> (m_Nt, 0.0);

    for (auto& heat_source : m_heat_sources) {
        for (int j {}; j < m_Ny; ++j) {
            for (int i {}; i < m_Nx; ++i) {

                float x_i = m_Xa + (static_cast<float>(i) * m_hx);
                float y_j = m_Ya + (static_cast<float>(j) * m_hy);
                if (heat_source.is_has_point(x_i, y_j)) {
                    (*m_heat_grid)[i][j] = heat_source.get_heat();
                }
            }
        }
    }

    const float hx_pow_minus1 = 1.0f / m_hx;
    const float hy_pow_minus1 = 1.0f / m_hy;
    hx_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hx, 2));
    hy_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hy, 2));

    m_hx_hy_2 = hx_pow_minus2 + hy_pow_minus2;
    m_a = 1.0f / (5.0f * m_hx_hy_2);
    m_b = static_cast<float>((1.0f / 2.0f) * (5.0f * hy_pow_minus2 - 1.0f * hx_pow_minus2));
    m_c = m_hx_hy_2 / 4;


    m_a_m512 = _mm512_set1_ps(m_a);
    m_b_m512 = _mm512_set1_ps(m_b);
    m_c_m512 = _mm512_set1_ps(m_c);


    std::cout << std::format("Address of begin of m_value_grid: {}", static_cast<const void*>(&(*m_value_grid)[0][0])) << std::endl;
    std::cout << std::format("Address of begin of m_previous_value_grid: {}", static_cast<const void*>(&(*m_previous_value_grid)[0][0])) << std::endl;
    std::cout << "created solver v3" << std::endl;
}
PoissonEquationSolver::PoissonEquationSolver(const int Nx, const int Ny, const int Nt) :
    m_Xa(0), m_Xb(4.0),
    m_Ya(0), m_Yb(4.0) {

    m_Nx = Nx;
    m_Ny = Ny;
    m_Nt = Nt;

    m_hx = (m_Xb - m_Xa) / static_cast<float>(m_Nx - 1);
    m_hy = (m_Yb - m_Ya) / static_cast<float>(m_Ny - 1);

    float Xs1 = m_Xa + (m_Xb - m_Xa) / 3;
    float Ys1 = m_Ya + (m_Yb - m_Ya) * 2 / 3;
    float Xs2 = m_Xa + (m_Xb - m_Xa) * 2 / 3;
    float Ys2 = m_Ya + (m_Yb - m_Ya) / 3;
    float R = static_cast<float>(0.1) * std::min(m_Xb - m_Xa, m_Yb - m_Ya);

    m_heat_sources.emplace_back(Xs1, Ys1, R, 0.1);
    m_heat_sources.emplace_back(Xs2, Ys2, R, -0.1);

    m_value_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));
    m_previous_value_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));
    m_heat_grid = std::make_shared<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> (m_Nx, std::vector<float, AlignedAllocator<float, 64>>(m_Ny, 0.0));


    m_deltas = std::vector<float> (m_Nt, 0.0);

    for (auto& heat_source : m_heat_sources) {
        for (int j {}; j < m_Ny; ++j) {
            for (int i {}; i < m_Nx; ++i) {

                float x_i = m_Xa + (static_cast<float>(i) * m_hx);
                float y_j = m_Ya + (static_cast<float>(j) * m_hy);
                if (heat_source.is_has_point(x_i, y_j)) {
                    (*m_heat_grid)[i][j] = heat_source.get_heat();
                }
            }
        }
    }

    const float hx_pow_minus1 = 1.0f / m_hx;
    const float hy_pow_minus1 = 1.0f / m_hy;
    hx_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hx, 2));
    hy_pow_minus2 = static_cast<float>(std::pow(1.0f / m_hy, 2));

    m_hx_hy_2 = hx_pow_minus2 + hy_pow_minus2;
    m_a = 1.0f / (5.0f * m_hx_hy_2);
    m_b = static_cast<float>((1.0f / 2.0f) * (5.0f * hy_pow_minus2 - 1.0f * hx_pow_minus2));
    m_c = m_hx_hy_2 / 4;


    m_a_m512 = _mm512_set1_ps(m_a);
    m_b_m512 = _mm512_set1_ps(m_b);
    m_c_m512 = _mm512_set1_ps(m_c);

    std::cout << "created solver v3" << std::endl;
}
float PoissonEquationSolver::calc_new_value(const float F_im1_jm1,  const float F_im1_j, const float F_im1_jp1,
                                            const float F_i_jm1,    /*target*/           const float F_i_jp1,
                                            const float F_ip1_jm1,  const float F_ip1_j, const float F_ip1_jp1,

                                                                    const float P_im1_j,
                                            const float P_i_jm1,    const float P_i_j,   const float P_i_jp1,
                                                                    const float P_ip1_j
                                                                    ) const {

    const auto first_line = m_b * (F_i_jm1 + F_i_jp1) + m_b * (F_im1_j + F_ip1_j);
    const auto second_line = m_c * (F_im1_jm1 + F_im1_jp1 + F_ip1_jm1 + F_ip1_jp1);
    const auto third_line = 2.0f * P_i_j + (1.0f / 4.0f) * (P_im1_j + P_ip1_j + P_i_jm1 + P_i_jp1);
    const auto result = m_a * (first_line + second_line + third_line);

    return result;
}

__m512 PoissonEquationSolver::calc_new_value (const __m512 F_im1_jm1, const __m512 F_im1_j, const __m512 F_im1_jp1,
                                              const __m512 F_i_jm1,   /*target*/            const __m512 F_i_jp1,
                                              const __m512 F_ip1_jm1, const __m512 F_ip1_j, const __m512 F_ip1_jp1,

                                                                      const __m512 P_im1_j,
                                              const __m512 P_i_jm1,   const __m512 P_i_j,   const __m512 P_i_jp1,
                                                                      const __m512 P_ip1_j
                                                                      ) const {

    const __m512 first_line  = _mm512_fmadd_ps(m_b_m512, _mm512_add_ps(F_i_jm1, F_i_jp1), _mm512_mul_ps(m_b_m512, _mm512_add_ps(F_im1_j, F_ip1_j)));
    const __m512 second_line = _mm512_mul_ps(m_c_m512, _mm512_add_ps(_mm512_add_ps(F_im1_jm1, F_im1_jp1), _mm512_add_ps(F_ip1_jm1, F_ip1_jp1)));
    const __m512 third_line = _mm512_add_ps(_mm512_mul_ps(m_mul_2, P_i_j), _mm512_mul_ps(m_mul_025, _mm512_add_ps(P_im1_j, _mm512_add_ps(P_ip1_j, _mm512_mul_ps(P_i_jm1, P_i_jp1)))));
    const __m512 result = _mm512_mul_ps(m_a_m512, _mm512_add_ps(first_line, _mm512_add_ps(second_line, third_line)));

    return result;
};

__m512 PoissonEquationSolver::calc_new_value (
    __m512* p_F_im1_jm1, __m512* p_F_im1_j, __m512* p_F_im1_jp1,
    __m512* p_F_i_jm1,   /*target*/         __m512* p_F_i_jp1,
    __m512* p_F_ip1_jm1, __m512* p_F_ip1_j, __m512* p_F_ip1_jp1,

                         __m512* p_P_im1_j,
    __m512* p_P_i_jm1,   __m512* p_P_i_j,   __m512* p_P_i_jp1,
                         __m512* p_P_ip1_j
                         ) const {

    const __m512 first_line  = _mm512_fmadd_ps(m_b_m512, _mm512_add_ps(*p_F_i_jm1, *p_F_i_jp1), _mm512_mul_ps(m_b_m512, _mm512_add_ps(*p_F_im1_j, *p_F_ip1_j)));
    const __m512 second_line = _mm512_mul_ps(m_c_m512, _mm512_add_ps(_mm512_add_ps(*p_F_im1_jm1, *p_F_im1_jp1), _mm512_add_ps(*p_F_ip1_jm1, *p_F_ip1_jp1)));
    const __m512 third_line = _mm512_add_ps(_mm512_mul_ps(m_mul_2, *p_P_i_j), _mm512_mul_ps(m_mul_025, _mm512_add_ps(*p_P_im1_j, _mm512_add_ps(*p_P_ip1_j, _mm512_mul_ps(*p_P_i_jm1, *p_P_i_jp1)))));
    const __m512 result = _mm512_mul_ps(m_a_m512, _mm512_add_ps(first_line, _mm512_add_ps(second_line, third_line)));

    return result;
};

void PoissonEquationSolver::make_one_calc_vectorized_512(float &delta, int i, int j,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const {


    const __m512 F_im1_jm1      = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j-1));
    const __m512 F_im1_j        = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j));
    const __m512 F_im1_jp1      = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j+1));
    const __m512 F_i_jm1        = _mm512_loadu_ps((*previous_value_grid)[i].data() + (j-1));
    const __m512 F_i_jp1        = _mm512_loadu_ps((*previous_value_grid)[i].data() + (j+1));
    const __m512 F_ip1_jm1      = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j-1));
    const __m512 F_ip1_j        = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j));
    const __m512 F_ip1_jp1      = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j+1));
    const __m512 P_im1_j        = _mm512_loadu_ps((*m_heat_grid)[i-1].data() + (j));
    const __m512 P_i_jm1        = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j-1));
    const __m512 P_i_j          = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j));
    const __m512 P_i_jp1        = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j+1));
    const __m512 P_ip1_j        = _mm512_loadu_ps((*m_heat_grid)[i+1].data() + (j));

    __m512 result = calc_new_value(
        F_im1_jm1, F_im1_j,   F_im1_jp1,
        F_i_jm1,              F_i_jp1,
        F_ip1_jm1, F_ip1_j,   F_ip1_jp1,

        P_im1_j,
        P_i_jm1,   P_i_j,     P_i_jp1,
        P_ip1_j
    );


    _mm512_storeu_ps((*previous_value_grid)[i].data() + j, result);


    for (int k = j; k < j + 16; ++k) {
        delta = std::max(delta, std::abs((*previous_value_grid)[i][k] - (*value_grid)[i][k]));
    }
}

void PoissonEquationSolver::make_one_calc_vectorized_512(const int i, const int j,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const {
    std::cout << std::format("starting make_one_calc_vectorized_512()...") << std::endl;
    if (previous_value_grid == nullptr || value_grid == nullptr || value_grid->empty()) {
        std::cerr << "NIGGERS" << std::endl;
    }
    std::cout << std::format("Address of begin of value_grid: {}", static_cast<const void*>(&(*value_grid)[0][0])) << std::endl;
    std::cout << std::format("Address of begin of previous_value_grid: {}", static_cast<const void*>(&(*previous_value_grid)[0][0])) << std::endl;
    std::cout << std::format("i = {}, j = {}", i, j) << std::endl;
    std::cout << std::format("trying to load F vectors...") << std::endl;
    //TODO: тут таргет это нулевая строка, следовательно чтение с i = -1, нужно в FrontMovingAlongArray учесть это.
    const __m512 F_im1_jm1      = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j-1));
    const __m512 F_im1_j        = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j));
    const __m512 F_im1_jp1      = _mm512_loadu_ps((*previous_value_grid)[i-1].data() + (j+1));
    const __m512 F_i_jm1        = _mm512_loadu_ps((*previous_value_grid)[i].data() + (j-1));
    const __m512 F_i_jp1        = _mm512_loadu_ps((*previous_value_grid)[i].data() + (j+1));
    const __m512 F_ip1_jm1      = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j-1));
    const __m512 F_ip1_j        = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j));
    const __m512 F_ip1_jp1      = _mm512_loadu_ps((*previous_value_grid)[i+1].data() + (j+1));
    std::cout << std::format("F vectors are loaded successfully") << std::endl;
    const __m512 P_im1_j        = _mm512_loadu_ps((*m_heat_grid)[i-1].data() + (j));
    const __m512 P_i_jm1        = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j-1));
    const __m512 P_i_j          = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j));
    const __m512 P_i_jp1        = _mm512_loadu_ps((*m_heat_grid)[i].data() + (j+1));
    const __m512 P_ip1_j        = _mm512_loadu_ps((*m_heat_grid)[i+1].data() + (j));

    std::cout << std::format("calling calc_new_value()...") << std::endl;
    const __m512 result = calc_new_value(
        F_im1_jm1, F_im1_j,   F_im1_jp1,
        F_i_jm1,              F_i_jp1,
        F_ip1_jm1, F_ip1_j,   F_ip1_jp1,

        P_im1_j,
        P_i_jm1,   P_i_j,     P_i_jp1,
        P_ip1_j
    );


    _mm512_storeu_ps((*value_grid)[i].data() + j, result);
}

void PoissonEquationSolver::horizontal_step(
    float& delta,
    const int i,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const {
    for (int j = 1; j < m_Nx - 1; j += m_alignment_float) {
        make_one_calc_vectorized_512(delta, i, j, previous_value_grid, value_grid);
    }
}

void PoissonEquationSolver::horizontal_step(
    const int i,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const {
    std::cout << std::format("calc {}'th line in horizontal_step()", i) << std::endl;

    for (int j = 1; j < m_Nx - 1; j += m_alignment_float) {
        make_one_calc_vectorized_512(i, j, previous_value_grid, value_grid);
    }
}

// void PoissonEquationSolver::solve() {
//
//     for (int time {}; time < m_Nt; ++time) {
//         float delta = INT_MIN;
//         std::swap(m_previous_value_grid, m_value_grid);
//         // #pragma omp simd for reduction(max: delta)
//         // #pragma omp parallel for reduction(+:delta)
//         // #pragma omp parallel for
//         for (int i = 1; i < m_Ny - 1; ++i) {
//
//             horizontal_step(delta, i, m_previous_value_grid, m_value_grid);
//         }
//         m_deltas[time] = delta;
//     }
// }

void PoissonEquationSolver::solve() const {
    // float delta {};
    auto front = FrontMovingAlongArray(
        m_Ny,
        m_Nt,
        100,
        m_previous_value_grid,
        m_value_grid,
        [this](
        const int index,
            const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& a,
            const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& b) {
            this->horizontal_step(index, a, b);
        }
    );

    front.move_all_times();
}

/*
void PoissonEquationSolver::solve() {
    constexpr int front_height = 100;
    std::vector<bool> mask(front_height, false);

    for (int time {}; time < m_Nt; ++time) {
        float delta = INT_MIN;
        std::swap(m_previous_value_grid, m_value_grid);

        for (int i = 1 - front_height + 1; i < m_Ny - 1; ++i) {
            // Обновляем маску
            for (int j = 0; j < front_height; ++j) {
                int row = i + j;
                if (row <= 0 || row >= m_Ny - 1) {
                    mask[j] = false;
                } else if (i < front_height) {
                    mask[j] = (j <= i);  // наращиваем фронт
                } else if (i >= m_Ny - 1 - front_height) {
                    mask[j] = (j < m_Ny - 1 - i);  // сужаем фронт
                } else {
                    mask[j] = true;  // постоянный фронт
                }
            }

            // Применяем шаг по всем активным строкам
            for (int j = 0; j < front_height; ++j) {
                if (mask[j]) {
                    horizontal_step(delta, i + j, m_previous_value_grid, m_value_grid);
                }
            }
        }

        m_deltas[time] = delta;
    }
}
*/
void PoissonEquationSolver::export_grid_value_as_matrix(const std::string &file_path) const {
    std::ofstream output (file_path, std::ios::out);
    if (!output.is_open()) {
        std::cerr << "Failed to open file: " << "float2.dat" << std::endl;
        throw std::exception();
    }


    for (int i {}; i < m_value_grid->size(); ++i) {
        for (int j {}; j < (*m_value_grid)[i].size(); ++j) {
            output << (*m_value_grid)[i][j] << " ";
        }
        output << "\n";
    }

    output.close();
}
void PoissonEquationSolver::check_deltas() const {
    int k {};
    for (int time {1}; time < m_Nt; ++time) {
        if (!(m_deltas[time] <= m_deltas[time - 1])) {
            std::cerr << std::format("iter {}: delta >= previous delta ({} >= {}).", time, m_deltas[time], m_deltas[time - 1]) << std::endl;
            ++k;
        }
    }
    std::cerr << std::format("delta >= previous delta on {}/{} iterations", k, m_Nt - 1) << std::endl;


}