#include <iostream>
#include "../include/orbit.h"

int main() {
    // Earth's orbital elements (nu=0 means we start at periapsis, t0=0)
orion::KeplerianElements earth {
    .a     = 1.496e8,
    .e     = 0.0167,
    .i     = 0.0,
    .raan  = 0.0,
    .omega = 1.7966,
    .nu    = 0.0
};

orion::KeplerianElements mars {
    .a     = 2.279e8,    // km (1.524 AU)
    .e     = 0.0934,
    .i     = 0.03229,    // radians (~1.85 degrees)
    .raan  = 0.8650,     // radians
    .omega = 5.0004,     // radians
    .nu    = 0.0         // arbitrary start
};

double T = 365.25 * 24 * 3600; // one year in seconds

// At t=0, should match keplerian_to_state directly
auto sv0 = orion::propagate(earth, 0.0);
std::cout << "At t=0:\n";
std::cout << "  |r| = " << sv0.r.norm() << " km (expect 1.471e8)\n";
std::cout << "  |v| = " << sv0.v.norm() << " km/s (expect 30.29)\n";

// At t=T/2 (half orbit), Earth should be at aphelion
auto sv1 = orion::propagate(earth, T/2);
std::cout << "\nAt t=T/2 (aphelion):\n";
std::cout << "  |r| = " << sv1.r.norm() << " km (expect 1.521e8)\n";
std::cout << "  |v| = " << sv1.v.norm() << " km/s (expect ~29.3)\n";

// At t=T (full orbit), should be back near periapsis
auto sv2 = orion::propagate(earth, T);
std::cout << "\nAt t=T (full orbit):\n";
std::cout << "  |r| = " << sv2.r.norm() << " km (expect ~1.471e8)\n";
std::cout << "  |v| = " << sv2.v.norm() << " km/s (expect ~30.29)\n";

    return 0;
}