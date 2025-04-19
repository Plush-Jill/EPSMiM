//
// Created by plush-jill on 3/4/25.
//

#include "poisson_equation_solver.hpp"
#include <fstream>
#include <iostream>
#include <immintrin.h>

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

    // std::cout << "created solver v2" << std::endl;
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

    // std::cout << "created solver v2" << std::endl;
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
    // __m512 result = {};
#ifdef ONLY_ALIGNED_LOAD
                std::swap(p_F_im1_next, p_F_im1_current);
                std::swap(p_F_i_next, p_F_i_current);
                std::swap(p_F_ip1_next, p_F_ip1_current);
                std::swap(p_P_im1_next, p_P_im1_current);
                std::swap(p_P_i_next, p_P_i_current);
                std::swap(p_P_ip1_next, p_P_ip1_current);

                F_im1_next     = _mm512_load_ps((*m_previous_value_grid)[i-1].data() + (j-1) + m_alignment_float);     /*p_F_im1_next = &F_im1_next;*/  p_F_im1_jm1 = p_F_im1_current;
                F_i_next       = _mm512_load_ps((*m_previous_value_grid)[i].data() + (j-1) + m_alignment_float);       /*p_F_i_next = &F_i_next;*/      p_F_i_jm1 = p_F_i_current;
                F_ip1_next     = _mm512_load_ps((*m_previous_value_grid)[i+1].data() + (j-1) + m_alignment_float);     /*p_F_ip1_next = &F_ip1_next;*/  p_F_ip1_jm1 = p_F_ip1_current;
                P_im1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);               /*p_P_im1_next = &P_im1_next;*/
                P_i_next       = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);               /*p_P_i_next = &P_i_next;*/      p_P_i_jm1 = p_P_i_current;
                P_ip1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);               /*p_P_ip1_next = &P_ip1_next;*/

                F_im1_next     = _mm512_load_ps((*m_previous_value_grid)[i-1].data() + (j-1) + m_alignment_float);                           p_F_im1_next = &F_im1_next;
                F_i_next       = _mm512_load_ps((*m_previous_value_grid)[i].data() + (j-1) + m_alignment_float);                             p_F_i_next = &F_i_next;
                F_ip1_next     = _mm512_load_ps((*m_previous_value_grid)[i+1].data() + (j-1) + m_alignment_float);                           p_F_ip1_next = &F_ip1_next;
                P_im1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);                                     p_P_im1_next = &P_im1_next;
                P_i_next       = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);                                     p_P_i_next = &P_i_next;
                P_ip1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data() + (j-1) + m_alignment_float);                                     p_P_ip1_next = &P_ip1_next;

                const __m512 F_im1_j      = _mm512_mask_blend_ps(m_current_p1_mask, *p_F_im1_current, F_im1_next);   p_F_im1_j = &F_im1_j;
                const __m512 F_im1_jp1    = _mm512_mask_blend_ps(m_current_p2_mask, *p_F_im1_current, F_im1_next);   p_F_im1_jp1 = &F_im1_jp1;
                const __m512 F_i_jp1      = _mm512_mask_blend_ps(m_current_p2_mask, *p_F_i_current, F_i_next);       p_F_i_jp1 = &F_i_jp1;
                const __m512 F_ip1_j      = _mm512_mask_blend_ps(m_current_p1_mask, *p_F_ip1_current, F_ip1_next);   p_F_ip1_j = &F_ip1_j;
                const __m512 F_ip1_jp1    = _mm512_mask_blend_ps(m_current_p2_mask, *p_F_ip1_current, F_ip1_next);   p_F_ip1_jp1 = &F_ip1_jp1;
                const __m512 P_im1_j      = _mm512_mask_blend_ps(m_current_p1_mask, *p_P_im1_current, P_im1_next);   p_P_im1_j = &P_im1_j;
                const __m512 P_i_j        = _mm512_mask_blend_ps(m_current_p1_mask, *p_P_i_current, P_i_next);       p_P_i_j = &P_i_j;
                const __m512 P_i_jp1      = _mm512_mask_blend_ps(m_current_p2_mask, *p_P_i_current, P_i_next);       p_P_i_jp1 = &P_i_jp1;
                const __m512 P_ip1_j      = _mm512_mask_blend_ps(m_current_p1_mask, *p_P_ip1_current, P_ip1_next);   p_P_ip1_j = &P_ip1_j;

                __m512 result = calc_new_value(
                    p_F_im1_jm1, p_F_im1_j, p_F_im1_jp1, // <- <-- p_F_im1_j_next
                    p_F_i_jm1,              p_F_i_jp1,   // <- <-- p_F_i_j_next
                    p_F_ip1_jm1, p_F_ip1_j, p_F_ip1_jp1, // <- <-- p_F_ip1_j_next

                                 p_P_im1_j,              // <- <-- p_P_im1_j_next
                    p_P_i_jm1,   p_P_i_j,   p_P_i_jp1,   // <- <-- p_P_i_j_next
                                 p_P_ip1_j               // <- <-- p_P_ip1_j_next
                );
#endif

// #ifdef  UNALIGNED_LOAD

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

// #endif

    _mm512_storeu_ps((*previous_value_grid)[i].data() + j, result);


    for (int k = j; k < j + 16; ++k) {
        delta = std::max(delta, std::abs((*previous_value_grid)[i][k] - (*value_grid)[i][k]));
    }
}

void PoissonEquationSolver::horizontal_step(float delta, const int i,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const {
    for (int j = 1; j < m_Nx - 1; j += m_alignment_float) {
        make_one_calc_vectorized_512(delta, i, j, previous_value_grid, value_grid);
    }
}

void PoissonEquationSolver::solve() {

// #define ONLY_ALIGNED_LOAD
#define UNALIGNED_LOAD

#ifdef ONLY_ALIGNED_LOAD
    const __m512* p_F_im1_jm1 = nullptr;
    const __m512* p_F_im1_j = nullptr;
    const __m512* p_F_im1_jp1 = nullptr;
    const __m512* p_F_i_jm1 = nullptr;
    const __m512* p_F_i_jp1 = nullptr;
    const __m512* p_F_ip1_jm1 = nullptr;
    const __m512* p_F_ip1_j = nullptr;
    const __m512* p_F_ip1_jp1 = nullptr;
    const __m512* p_P_im1_j = nullptr;
    const __m512* p_P_i_jm1 = nullptr;
    const __m512* p_P_i_j = nullptr;
    const __m512* p_P_i_jp1 = nullptr;
    const __m512* p_P_ip1_j = nullptr;

    const __m512* p_F_im1_current = nullptr;
    const __m512* p_F_i_current = nullptr;
    const __m512* p_F_ip1_current = nullptr;
    const __m512* p_P_im1_current = nullptr;
    const __m512* p_P_i_current = nullptr;
    const __m512* p_P_ip1_current = nullptr;

    const __m512* p_F_im1_next = nullptr;
    const __m512* p_F_i_next = nullptr;
    const __m512* p_F_ip1_next = nullptr;
    const __m512* p_P_im1_next = nullptr;
    const __m512* p_P_i_next = nullptr;
    const __m512* p_P_ip1_next = nullptr;
#endif

    for (int time {}; time < m_Nt; ++time) {
        float delta = INT_MIN;
        std::swap(m_previous_value_grid, m_value_grid);
        // #pragma omp simd for reduction(max: delta)
        // #pragma omp parallel for reduction(+:delta)
        for (int i = 1; i < m_Ny - 1; ++i) {
            // #pragma omp simd
            #ifdef ONLY_ALIGNED_LOAD
            __m512 F_im1_next     = _mm512_load_ps((*m_previous_value_grid)[i-1].data());                               p_F_im1_next = &F_im1_next;
            __m512 F_i_next       = _mm512_load_ps((*m_previous_value_grid)[i].data());                                 p_F_i_next = &F_i_next;
            __m512 F_ip1_next     = _mm512_load_ps((*m_previous_value_grid)[i+1].data());                               p_F_ip1_next = &F_ip1_next;
            __m512 P_im1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data());                                         p_P_im1_next = &P_im1_next;
            __m512 P_i_next       = _mm512_load_ps((*m_heat_grid)[i+1].data());                                         p_P_i_next = &P_i_next;
            __m512 P_ip1_next     = _mm512_load_ps((*m_heat_grid)[i+1].data());                                         p_P_ip1_next = &P_ip1_next;
            #endif

            horizontal_step(delta, i, m_previous_value_grid, m_value_grid);
            // for (int k {}; k < m_Nx - 1; ++k) {
            //     if ((*m_value_grid)[i][k] != 0) {
            //         std::cerr << "Non zero value" << std::endl;
            //     }
            // }
        }
        m_deltas[time] = delta;
    }
}
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