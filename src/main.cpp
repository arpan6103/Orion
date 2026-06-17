#include <iostream>
#include "../include/orbit.h"

int main() {
    // Earth's orbital elements (nu=0 means we start at periapsis, t0=0)
orion::KeplerianElements earth {
    .a=1.496e8, .e=0.0167, .i=0.0,
    .raan=0.0, .omega=1.7966, .nu=0.0
};
orion::KeplerianElements mars {
    .a=2.279e8, .e=0.0934, .i=0.03229,
    .raan=0.8650, .omega=5.0004, .nu=0.0
};

double dt = 259.0 * 24.0 * 3600.0;

// Mars orbital period
double Tm = 2*M_PI*std::sqrt(std::pow(mars.a,3)/orion::MU_SUN);

// Departure: Earth at t=0
auto sv_earth = orion::propagate(earth, 0.0);

// Arrival: Mars at t=dt, but Mars began 0.45 of an orbit earlier
// so the transfer sweeps a realistic ~161 degrees
auto sv_mars = orion::propagate(mars, dt, -0.45*Tm);

orion::Vec3 r1 = sv_earth.r;
orion::Vec3 r2 = sv_mars.r;

auto res = orion::lambert(r1, r2, dt, orion::MU_SUN, true);

double dv1 = (res.v1 - sv_earth.v).norm();
double dv2 = (sv_mars.v - res.v2).norm();

std::cout << "\nLambert test (realistic geometry):\n";
std::cout << "  converged: " << (res.converged ? "yes" : "no") << "\n";
std::cout << "  total dv = " << dv1+dv2 << " km/s (expect ~9)\n";

    return 0;
}