#include <iostream>
#include "../include/orbit.h"

int main() {
    orion::KeplerianElements test {
    .a     = 1.0e8,
    .e     = 0.3,
    .i     = 0.5,     // ~28.6 degrees - nonzero, this is the key
    .raan  = 0.8,
    .omega = 1.2,
    .nu    = 0.6
};

auto sv2 = orion::keplerian_to_state(test);
auto back = orion::state_to_keplerian(sv2);

std::cout << "\nRound trip test (i != 0):\n";
std::cout << "  a:     " << back.a     << " (expect 1.0e8)\n";
std::cout << "  e:     " << back.e     << " (expect 0.3)\n";
std::cout << "  i:     " << back.i     << " (expect 0.5)\n";
std::cout << "  raan:  " << back.raan  << " (expect 0.8)\n";
std::cout << "  omega: " << back.omega << " (expect 1.2)\n";
std::cout << "  nu:    " << back.nu    << " (expect 0.6)\n";

    return 0;
}