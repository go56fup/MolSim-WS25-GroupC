//
// Created by Gabriel Ribeiro Fernandes on 10.01.26.
//

#include <vector>
#include <__random/random_device.h>
#include <random>

struct ParticleSystem {
    using vector = std::vector<double>;
    using size_type = std::uint32_t;

    vector x, y, z;
    vector vx, vy, vz;
    vector fx, fy, fz;
    vector old_fx, old_fy, old_fz;
    vector m;
    vector sigma;
    vector epsilon;

    constexpr ParticleSystem(size_t  particles_num) {
        x.assign(particles_num, 0.0);
        y.assign(particles_num, 0.0);
        z.assign(particles_num, 0.0);

        vx.assign(particles_num, 0.0);
        vy.assign(particles_num, 0.0);
        vz.assign(particles_num, 0.0);

        fx.assign(particles_num, 0.0);
        fy.assign(particles_num, 0.0);
        fz.assign(particles_num, 0.0);

        old_fx.assign(particles_num, 0.0);
        old_fy.assign(particles_num, 0.0);
        old_fz.assign(particles_num, 0.0);

        m.assign(particles_num, 0.0);
        sigma.assign(particles_num, 0.0);
        epsilon.assign(particles_num, 0.0);
    }
    constexpr ParticleSystem( ) {
            x = {};
            y = {};
            z = {};

            vx = {};
            vy = {};
            vz = {};

            fx = {};
            fy = {};
            fz = {};

            old_fx = {};
            old_fy = {};
            old_fz = {};

            m = {};
            sigma = {};
            epsilon = {};
    }

public:

    constexpr vec particle_coord_diff(size_type p1, size_type p2) {

            return vec{x[p1] - x[p2], y[p1] - y[p2], z[p1] - z[p2]};
    }
};
