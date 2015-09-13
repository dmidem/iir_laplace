#include <stdio.h>

#include "filter.h"
#include "dumper.h"

int main()
{
    double snum[] = { 1 };
    double sden[] = { 1, 1 };

    double dT = 0.01;
    int L = 1;
    double tmax = 0.2;

    int nsnum = sizeof(snum) / sizeof(snum[0]);
    int nsden = sizeof(sden) / sizeof(sden[0]);
    double dt = dT / L;

    IIR_FILTER *iir = laplace_nd_filter(snum, nsnum, sden, nsden, dt, L);

    print_polynom_nd(snum, nsnum, sden, nsden, 's', "Laplace H(s)");
    printf("dT: %f, L: %d, dt: %f\n", dT, L, dt);
    print_polynom_nd(iir->b, iir->m + 1, iir->a, iir->n + 1, 'z', "IIR H(z)");

    printf("t\tx\ty\n");
    double t, x, y;
    for (t = 0; t < tmax; t += dt) {
	x = 1;
	y = iir_filter_next(iir, x);
	printf("%f\t%f\t%f\n", t, x, y);
    }

    iir_filter_free(iir);

    return 0;
}
