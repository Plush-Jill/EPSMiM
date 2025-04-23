#ifndef POISSON_EQUATION_SOLVER_HPP
#define POISSON_EQUATION_SOLVER_HPP
#include <boost/json.hpp>
#include <filesystem>
#include "../common/heat_source_circle.hpp"
#include <format>
#include "../common/aligned_allocator.hpp"
#include <immintrin.h>



class PoissonEquationSolver {
private:
    const float m_Xa;
    const float m_Xb;
    const float m_Ya;
    const float m_Yb;

    long m_Nx;
    long m_Ny;
    long m_Nt;

    float m_hx;
    float m_hy;

    float hx_pow_minus2;
    float hy_pow_minus2;
    float m_hx_hy_2;
    float m_a;
    float m_b;
    float m_c;

    __m512 m_a_m512;
    __m512 m_b_m512;
    __m512 m_c_m512;
    const __m512 m_mul_2  = _mm512_set1_ps(2.0f);
    const __m512 m_mul_025 = _mm512_set1_ps(0.25f);

    const int m_alignment_byte = 64;
    const int m_alignment_float = 16;

    const unsigned short m_current_p1_mask = 0b1111'1111'1111'1110;
    const unsigned short m_current_p2_mask = 0b1111'1111'1111'1100;

    std::vector<HeatSourceCircle> m_heat_sources;
    std::shared_ptr<std::vector<std::vector<std::pair<float, float>>>> m_grid; //first - heat, second - function value
    std::shared_ptr<std::vector<std::vector<std::pair<float, float>>>> m_previous_grid; //first - heat, second - function value

    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_value_grid;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_previous_value_grid;
    std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>> m_heat_grid;

    std::vector<float> m_deltas;
    bool m_export;

public:
    [[nodiscard]] bool is_should_export() const;
    explicit PoissonEquationSolver(const std::string& config_file);

    [[nodiscard]] float calc_new_value (
        float F_im1_jm1, float F_im1_j, float F_im1_jp1,
        float F_i_jm1,                  float F_i_jp1,
        float F_ip1_jm1, float F_ip1_j, float F_ip1_jp1,

                         float P_im1_j,
        float P_i_jm1,   float P_i_j,   float P_i_jp1,
                         float P_ip1_j
            ) const;

    [[nodiscard]] __m512 calc_new_value (
        __m512 F_im1_jm1, __m512 F_im1_j,  __m512 F_im1_jp1,
        __m512 F_i_jm1,                    __m512 F_i_jp1,
        __m512 F_ip1_jm1, __m512 F_ip1_j, __m512 F_ip1_jp1,

                          __m512 P_im1_j,
        __m512 P_i_jm1,   __m512 P_i_j,   __m512 P_i_jp1,
                          __m512 P_ip1_j

        ) const;

    [[nodiscard]] __m512 calc_new_value (
        __m512* p_F_im1_jm1, __m512* p_F_im1_j, __m512* p_F_im1_jp1,
        __m512* p_F_i_jm1,   /*target*/         __m512* p_F_i_jp1,
        __m512* p_F_ip1_jm1, __m512* p_F_ip1_j, __m512* p_F_ip1_jp1,

                             __m512* p_P_im1_j,
        __m512* p_P_i_jm1,   __m512* p_P_i_j,   __m512* p_P_i_jp1,
                             __m512* p_P_ip1_j
            ) const;

    void make_one_calc_vectorized_512(float &delta, int i, int j,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid,
    const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const;

    void horizontal_step(float delta, int i, const std::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& previous_value_grid, const std
                         ::shared_ptr<std::vector<std::vector<float, AlignedAllocator<float, 64>>>>& value_grid) const;

    void solve();

    void export_grid_value_as_matrix(const std::string& file_path) const;

    void check_deltas() const;
};


#endif //POISSON_EQUATION_SOLVER_HPP
