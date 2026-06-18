#include "../include/orbit.h"
#include<iostream>
#include <algorithm>

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

    double solve_kepler(double M, double e, double tol, int max_iter){
        M=std::fmod(M,2*M_PI);
        if(M<0){
            M+=2*M_PI;
        }
        double E=M;
        for(int i=0;i<max_iter;i++){
            double f=E-e*std::sin(E)-M;
            double fp=1.0-e*std::cos(E);
            double dE=-f/fp;
            E+=dE;
            if(std::abs(dE)<tol) break;
        }
        return E;
    }

    double mean_anomaly_to_nu(double M, double e){
        double E=solve_kepler(M,e);
        double nu = 2.0 * std::atan2(
            std::sqrt(1+e) * std::sin(E/2),
            std::sqrt(1-e) * std::cos(E/2)
        );
        if(nu<0){
            nu+=2*M_PI;
        }
        return nu;
    }

    double compute_mean_anomaly(double t,double t0,double a,double mu){
        double n=std::sqrt(mu/(a*a*a));
        return n*(t-t0);
    }

    StateVector propagate(const KeplerianElements& k, double t, double t0, double mu){
        double M=compute_mean_anomaly(t,t0,k.a,mu);
        double nu=mean_anomaly_to_nu(M,k.e);
        KeplerianElements k_at_t=k;
        k_at_t.nu=nu;

        return keplerian_to_state(k_at_t,mu);
    }

    LambertResult lambert(const Vec3& r1,const Vec3& r2, double dt, double mu, bool prograde){
        LambertResult result;
        result.converged=false;

        //step1 geometry consts
        double r1_mag=r1.norm();
        double r2_mag=r2.norm();

        Vec3 r12=r2-r1;
        double c=r12.norm();            // chord
        double s=(r1_mag+r2_mag+c)/2.0; // semi perimeter

        Vec3 r1xr2=r1.cross(r2);
        double sin_dnu = std::sqrt(1.0 - std::pow(r1.dot(r2)/(r1_mag*r2_mag), 2));
        if (prograde ? (r1xr2.z < 0) : (r1xr2.z >= 0)) sin_dnu = -sin_dnu;

        double lam2=1.0-c/s;
        double lambda=std::sqrt(lam2);
        if (sin_dnu<0) lambda=-lambda;

        double rho=(r1_mag-r2_mag)/c;
        double sigma=std::sqrt(1.0-rho*rho);
        double gamma=std::sqrt(mu*s/2.0);

        double T_target=dt*std::sqrt(2.0*mu/(s*s*s));

        //step2 initial guess for x
        double beta_0=2.0*std::asin(std::sqrt(lam2));
        if(lambda<0) beta_0=-beta_0;

        double T0=(std::acos(lambda)+lambda*std::sqrt(1.0-lam2));
        double x=(T_target>T0) ? 0.0 : 0.5;

        //step3 Newton-Raphson
        for(int iter=0;iter<50;iter++){
            x=std::max(-0.99,std::min(0.99,x));

            double y_ = std::sqrt(1.0 - lam2*(1.0 - x*x));

            // psi must be in [0, pi] for elliptic transfers
            double psi;
            if (x < 1.0) {
                // elliptic
                psi = std::acos(x*y_ + lambda*(1.0 - x*x));
            } else {
                psi = 0.0; // hyperbolic not needed for Earth-Mars
            }

            double T_x = (psi/std::sqrt(1.0 - x*x) - x + lambda*y_) / (1.0 - x*x);
            double dT  = (3.0*x*T_x - 2.0*(1.0 - lambda*x/y_*(1.0-x*x) - lambda*y_)) / (1.0 - x*x);

            double err=T_x-T_target;
            if(std::abs(err)<1e-10){
                result.converged=true;
                break;
            }
            if(std::abs(dT)<1e-15) break;
            x-=err/dT;
        }

        //step4 recover velocities from x
        double y = std::sqrt(1.0 - lam2*(1.0 - x*x));

        double gamma_v = std::sqrt(mu * s / 2.0);
        double rho_v   = (r1_mag - r2_mag) / c;
        double sigma_v = std::sqrt(1.0 - rho_v*rho_v);

        // Radial and tangential velocity components
        double vr1 =  gamma_v * ((lambda*y - x) - rho_v*(lambda*y + x)) / r1_mag;
        double vr2 = -gamma_v * ((lambda*y - x) + rho_v*(lambda*y + x)) / r2_mag;
        double vt1 =  gamma_v * sigma_v * (y + lambda*x) / r1_mag;
        double vt2 =  gamma_v * sigma_v * (y + lambda*x) / r2_mag;

        // Unit vectors: radial, and tangential built from unit normal
        Vec3 ir1 = r1 * (1.0/r1_mag);
        Vec3 ir2 = r2 * (1.0/r2_mag);

        Vec3 ih = r1.cross(r2);
        ih = ih * (1.0/ih.norm());
        if (lambda < 0) ih = ih * -1.0;   // flip normal for retrograde

        Vec3 it1 = ih.cross(ir1);
        Vec3 it2 = ih.cross(ir2);

        result.v1 = ir1*vr1 + it1*vt1;
        result.v2 = ir2*vr2 + it2*vt2;

        return result;
    }
}