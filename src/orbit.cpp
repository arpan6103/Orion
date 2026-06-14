#include "../include/orbit.h"

namespace orion{
    StateVector keplerian_to_state(const KeplerianElements& k,double mu){
        double p=k.a*(1.0-k.e*k.e);
        double r_mag=p/(1.0+k.e*std::cos(k.nu));

        Vec3 r_pqw={r_mag*std::cos(k.nu), r_mag*std::sin(k.nu), 0.0};

        Vec3 v_pqw={-std::sqrt(mu/p)*std::sin(k.nu), std::sqrt(mu/p)*(k.e+std::cos(k.nu)), 0.0};

        double cO=std::cos(k.raan) , sO=std::sin(k.raan);
        double ci=std::cos(k.i) , si=std::sin(k.i);
        double co=std::cos(k.omega) , so=std::sin(k.omega);

        Vec3 Qx={cO*co-sO*so*si, -cO*so-sO*co*ci, sO*si};
        Vec3 Qy={sO*co+cO*so*ci, -sO*so+cO*co*ci, -cO*si};
        Vec3 Qz={so*si, co*si, ci};

        auto rotate=[&](const Vec3& v)->Vec3{
            return {Qx.dot(v), Qy.dot(v), Qz.dot(v)};
        };

        return {rotate(r_pqw),rotate(v_pqw)};
    }
}