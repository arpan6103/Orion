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

        Vec3 Qx={cO*co-sO*so*ci, -cO*so-sO*co*ci, sO*si};
        Vec3 Qy={sO*co+cO*so*ci, -sO*so+cO*co*ci, -cO*si};
        Vec3 Qz={so*si, co*si, ci};

        auto rotate=[&](const Vec3& v)->Vec3{
            return {Qx.dot(v), Qy.dot(v), Qz.dot(v)};
        };

        return {rotate(r_pqw),rotate(v_pqw)};
    }

    KeplerianElements state_to_keplerian(const StateVector& sv,double mu){
        Vec3 r=sv.r;
        Vec3 v=sv.v;
        double r_mag=r.norm();
        double v_mag=v.norm();

        Vec3 h=r.cross(v);
        double h_mag=h.norm();

        Vec3 k={0.0,0.0,1.0};
        Vec3 n=k.cross(h);
        double n_mag=n.norm();

        Vec3 e_vec=(r*(v_mag*v_mag-mu/r_mag)-v*r.dot(v))*(1.0/mu);
        double e=e_vec.norm();

        double energy=(v_mag*v_mag/2.0)-(mu/r_mag);
        double a=-mu/(2.0*energy);

        double i=std::acos(h.z/h_mag);

        double raan=std::acos(n.x/n_mag);
        if(n.y<0) raan=2*M_PI-raan;

        double omega=std::acos(n.dot(e_vec)/(n_mag*e));
        if(e_vec.z<0) omega=2*M_PI-omega;

        double nu=std::acos(e_vec.dot(r)/(e*r_mag));
        if(r.dot(v)<0) nu=2*M_PI-nu;

        return {a,e,i,raan,omega,nu};
    }
}