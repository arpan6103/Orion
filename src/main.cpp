#include <iostream>
#include "../include/orbit.h"

int main() {
    // Earth's orbit around the Sun (approximately)
    orion::KeplerianElements earth {
        .a     = 1.496e8,      // 1 AU in km
        .e     = 0.0167,       // nearly circular
        .i     = 0.0,          // reference plane
        .raan  = 0.0,
        .omega = 1.7966,       // ~102.9 degrees
        .nu    = 0.0           // at perihelion
    };

    auto sv = orion::keplerian_to_state(earth);

    std::cout << "Earth position (km):\n";
    std::cout << "  x = " << sv.r.x << "\n";
    std::cout << "  y = " << sv.r.y << "\n";
    std::cout << "  z = " << sv.r.z << "\n";

    std::cout << "Earth velocity (km/s):\n";
    std::cout << "  vx = " << sv.v.x << "\n";
    std::cout << "  vy = " << sv.v.y << "\n";
    std::cout << "  vz = " << sv.v.z << "\n";

    // Sanity check: orbital speed should be ~29.8 km/s
    std::cout << "\n|v| = " << sv.v.norm() << " km/s (expect ~29.8)\n";
    // Distance should be ~perihelion = a*(1-e) = 1.471e8 km
    std::cout << "|r| = " << sv.r.norm() << " km (expect ~1.471e8)\n";

    return 0;
}