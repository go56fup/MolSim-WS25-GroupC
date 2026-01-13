//
// Created by Gabriel Ribeiro Fernandes on 10.01.26.
//
#include <valarray>
#include "vec_3d.hpp"
#include "soa.h"

class particle
/// @cond DO_NOT_DOCUMENT
/// @endcond
{
    using vec = vec_3d<double>;
private:
    using vec_args = std::tuple<double &, double &, double &>;

public:
    /// Position of the particle
    vec x{0, 0, 0};

    /// Velocity of the particle
    vec v{0, 0, 0};

    /// Force effective on this particle
    vec f{0, 0, 0};

    /// Force which was effective on this particle
    vec old_f{0, 0, 0};

    /// Mass of this particle
    double m{};

    double sigma{};
    double epsilon{};

    constexpr particle(
            vec x_arg, vec v_arg, double m_arg, double sigma_, double epsilon_
    ) noexcept
            : x(x_arg)
            , v(v_arg)
            , m(m_arg)
            , sigma(sigma_)
            , epsilon(epsilon_){}
};

class particle_container {
    using cell = std::vector<particle>;
    using value_type = cell;
    using size_type = std::uint32_t;
    using difference_type = std::int32_t;
    using index = vec_3d<size_type>;

    using cell_soa = std::vector<size_type>;

public:
    std::vector<cell> grid;
    std::vector<cell_soa> grid_soa;

    ParticleSystem particle_system{};

    vec domain_;
    index grid_size_;
    double cutoff_radius_;
    std::size_t size_ = 0;

    size_type div_round_up(double dividend, double divisor) noexcept {
        return static_cast<particle_container::size_type>(std::ceil(dividend / divisor));
    }

    constexpr particle_container(vec &&domain_arg, double cutoff_radius_arg)
            : domain_{std::forward<vec>(domain_arg)}, grid_size_{div_round_up(domain_.x, cutoff_radius_arg)
            ,div_round_up(domain_.y, cutoff_radius_arg), div_round_up(domain_.z, cutoff_radius_arg)}
            , cutoff_radius_{cutoff_radius_arg} {
        // Prevent implicit widening later, if x * y * z overflows for uint32, this takes care of it
        const auto grid_x = static_cast<std::size_t>(grid_size_.x);
        const auto grid_y = static_cast<std::size_t>(grid_size_.y);
        const auto grid_z = static_cast<std::size_t>(grid_size_.z);
        grid.resize(grid_x * grid_y * grid_z);
        grid_soa.resize(grid_x * grid_y * grid_z);
    }

    constexpr size_type linear_index(
            size_type x, size_type y,
            size_type z
    ) const noexcept {
        const auto result = (x * grid_size_.y * grid_size_.z) + (y * grid_size_.z) + z;
        return result;
    }

    constexpr size_type pos_to_linear_index(const vec& pos
    ) const noexcept {
        const auto x = static_cast<size_type>(pos.x / cutoff_radius_);
        const auto y = static_cast<size_type>(pos.y / cutoff_radius_);
        const auto z = static_cast<size_type>(pos.z / cutoff_radius_);
        return linear_index(x, y, z);
    }

    template <typename... Args>
    constexpr void emplace(vec&& position, Args&&... args) {
        grid[pos_to_linear_index(position)].emplace_back(
                std::forward<vec>(position), std::forward<Args>(args)...
        );
        ++size_;
    }

    template <typename... Args>
    constexpr void emplace_soa(size_type idx) {
        grid_soa[pos_to_linear_index(vec{particle_system.x[idx],
                                         particle_system.y[idx], particle_system.z[idx]})].emplace_back(idx);
        ++size_;
    }

    template <typename... Args>
    void emplace_random(const vec& min, const vec& max, Args&&... args) {
        // Static engine ensures we don't re-seed on every function call
        static std::random_device rd;
        static std::mt19937 gen(rd());

        // Use uniform_real_distribution<double> for high precision
        std::uniform_real_distribution<double> disX(min.x, max.x);
        std::uniform_real_distribution<double> disY(min.y, max.y);
        std::uniform_real_distribution<double> disZ(min.z, max.z);

        // Generate the random 3D position
        vec random_pos{
                disX(gen),
                disY(gen),
                disZ(gen)
        };

        // Forward the position and all extra arguments to your emplace function
        emplace(std::move(random_pos), std::forward<Args>(args)...);
    }

    template <typename... Args>
    void emplace_random_soa(size_t n, const double min, const double max,
                            vec velocity, double m, double sigma, double epsilon) {
        // 1. Resize vectors to the count
        particle_system.x.resize(n);
        particle_system.y.resize(n);
        particle_system.z.resize(n);

        particle_system.vx.resize(n);
        particle_system.vy.resize(n);
        particle_system.vz.resize(n);

        particle_system.fx.resize(n);
        particle_system.fy.resize(n);
        particle_system.fz.resize(n);

        particle_system.old_fx.resize(n);
        particle_system.old_fy.resize(n);
        particle_system.old_fz.resize(n);

        particle_system.m.resize(n);
        particle_system.sigma.resize(n);
        particle_system.epsilon.resize(n);

        // 2. Setup Random Engine (Mersenne Twister)
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(min, max);

        // 3. Fill the SoA
        for (size_type i = 0; i < n; ++i) {
            particle_system.x[i] = dis(gen);
            particle_system.y[i] = dis(gen);
            particle_system.z[i] = dis(gen);

            particle_system.vx[i] = velocity.x;
            particle_system.vy[i] = velocity.y;
            particle_system.vz[i] = velocity.z;

            particle_system.m[i] = m;
            particle_system.sigma[i] = sigma;
            particle_system.epsilon[i] = epsilon;
            emplace_soa(i);
        }

    }
};

using size_type = std::uint32_t;


double calc_norm(const size_type p1, const size_type p2, const ParticleSystem& sys) noexcept {
    double dx = sys.x[p1] - sys.x[p2];
    double dy = sys.x[p1] - sys.x[p2];
    double dz = sys.x[p1] - sys.x[p2];

    return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));
}

vec lennard_jones_force_soa(const size_type p1, const size_type p2, const particle_container& particleContainer) noexcept {

    ParticleSystem sys = particleContainer.particle_system;

    const double norm = calc_norm(p1, p2, sys);
    const double sigma = (sys.sigma[p1] + sys.sigma[p2]) / 2;
    const double eps = std::sqrt(sys.epsilon[p1] * sys.epsilon[p2]);
    const double scaling_factor =
            24 * eps / std::pow(norm, 2) * (std::pow(sigma / norm, 6) - 2 * (std::pow(sigma / norm, 12)));
    const auto result = scaling_factor * sys.particle_coord_diff(p1, p2);
    return result;
}

vec lennard_jones_force(const particle& p1, const particle& p2) noexcept {
    const auto& xi = p1.x;
    const auto& xj = p2.x;
    const double norm = (xi - xj).euclidian_norm();
    const double sigma = (p1.sigma + p2.sigma) / 2;
    const double eps = std::sqrt(p1.epsilon * p2.epsilon);
    const double scaling_factor =
            24 * eps / std::pow(norm, 2) * (std::pow(sigma / norm, 6) - 2 * (std::pow(sigma / norm, 12)));
    const auto result = scaling_factor * (xj - xi);
    return result;
}
