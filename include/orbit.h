#pragma once
#include<cmath>

namespace orion{
    const double MU_SUN=1.32712440018e11;

    struct Vec3{
        double x,y,z;

        Vec3 operator +(const Vec3& o) const {return {x+o.x, y+o.y, z+o.z};}
        Vec3 operator -(const Vec3& o) const {return {x-o.x, y-o.y, z-o.z};}
        Vec3 operator *(double s) const {return {x*s, y*s, z*s};}

        double dot(const Vec3& o) const {return x*o.x+y*o.y+z*o.z;}
        Vec3 cross(const Vec3& o)const {
            return {y*o.z-z*o.y,
                    z*o.x-x*o.z,
                    x*o.y-y*o.x
            };
        }
        double norm() const {return std::sqrt(x*x+y*y+z*z);}
        Vec3 unit() const {
            double n=norm();
            return {x/n,y/n,z/n};
        }
    };

    struct StateVector{
        Vec3 r;  //pos
        Vec3 v;  // vel
    };

    struct KeplerianElements{
        double a;  
        double e;  
        double i;
        double raan;
        double omega;
        double nu;
    };

    StateVector keplerian_to_state(const KeplerianElements& k,double mu=MU_SUN);

    KeplerianElements state_to_keplerian(const StateVector& sv,double mu=MU_SUN);

    double solve_kepler(double M,double e, double tol=1e-10, int max_iter=100);
    double mean_anomaly_to_nu(double M, double e);
    double compute_mean_anomaly(double t,double t0, double a, double mu=MU_SUN);
    StateVector propagate(const KeplerianElements& k, double t, double t0=0.0, double mu=MU_SUN);
}