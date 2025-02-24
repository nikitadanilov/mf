#include <cstdint>
#include <cmath>
#include "stdio.h"

#include "rat.hpp"

int main() {
    using rat = rat_t<int64_t, uint8_t>;
    double epsilon = pow(2.0, -60.0);
    for (double i = -2.0; i < 2.0; i += 0.1) {
        for (double j = 0.1; j < 2.0; j += 0.1) {
            rat ri = rat::fromdouble(i);
            rat rj = rat::fromdouble(j);
            if (fabs(i * j - (ri * rj).todouble()) > epsilon) {
                double rr = (ri * rj).todouble();
                printf("%.1f * %.1f = %.20f != %.20f (%g)\n", i, j, i * j, rr, i * j - rr);
            }
            if (j != 0.0 && fabs(i / j - (ri / rj).todouble()) > epsilon) {
                double dd = (ri / rj).todouble();
                printf("%.1f / %.1f = %.20f != %.20f (%g)\n", i, j, i / j, dd, i / j - dd);
            }
        }
    }
    return 0;
}
