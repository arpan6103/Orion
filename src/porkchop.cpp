#include "../include/orbit.h"
#include <iostream>
#include <fstream>
#include <cmath>

int main() {
    using namespace orion;

    KeplerianElements earth {
        .a=1.496e8, .e=0.0167, .i=0.0,
        .raan=0.0, .omega=1.7966, .nu=0.0
    };
    KeplerianElements mars {
        .a=2.279e8, .e=0.0934, .i=0.03229,
        .raan=0.8650, .omega=5.0004, .nu=0.0
    };

    const double DAY = 24.0 * 3600.0;

    // Sweep ranges (in days from epoch)
    double dep_start = 0;
    double dep_end   = 800;     // 400 days of possible departure dates
    double arr_start = 150;
    double arr_end   = 1100;     // arrival dates
    int    grid      = 100;     // 100x100 grid

    double dep_step = (dep_end - dep_start) / grid;
    double arr_step = (arr_end - arr_start) / grid;

    std::ofstream out("porkchop.csv");
    out << "departure_day,arrival_day,tof_day,total_dv\n";

    double min_dv = 1e9;
    double best_dep = 0, best_arr = 0;

    for (int i = 0; i < grid; i++) {
        double dep_day = dep_start + i * dep_step;
        double t_dep = dep_day * DAY;

        for (int j = 0; j < grid; j++) {
            double arr_day = arr_start + j * arr_step;
            double t_arr = arr_day * DAY;

            double tof = t_arr - t_dep;
            if (tof <= 0) continue;   // arrival must be after departure

            auto sv_earth = propagate(earth, t_dep);
            auto sv_mars  = propagate(mars,  t_arr);

            auto res = lambert(sv_earth.r, sv_mars.r, tof);
            if (!res.converged) continue;

            double dv1 = (res.v1 - sv_earth.v).norm();
            double dv2 = (sv_mars.v - res.v2).norm();
            double dv  = dv1 + dv2;

            out << dep_day << "," << arr_day << ","
                << tof/DAY << "," << dv << "\n";

            if (dv < min_dv) {
                min_dv = dv;
                best_dep = dep_day;
                best_arr = arr_day;
            }
        }
    }

    out.close();

    std::cout << "Porkchop sweep complete.\n";
    std::cout << "Grid: " << grid << "x" << grid << "\n";
    std::cout << "Minimum delta-V: " << min_dv << " km/s\n";
    std::cout << "  Departure day: " << best_dep << "\n";
    std::cout << "  Arrival day:   " << best_arr << "\n";
    std::cout << "  Transfer time: " << (best_arr - best_dep) << " days\n";
    std::cout << "Data written to porkchop.csv\n";

    return 0;
}