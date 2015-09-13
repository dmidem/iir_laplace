#include <stdio.h>
#include <math.h>

#include "filter.h"
#include "dumper.h"

double f(const double t)
{
    return sin(t) * sin(t * 0.1);
}

int main()
{
    double snum[] = { 1 };
    double sden[] = { 1 };

    double dT = 0.1;
    int L = 10;
    double tmax = 0.2;

    int nsnum = sizeof(snum) / sizeof(snum[0]);
    int nsden = sizeof(sden) / sizeof(sden[0]);
    double dt = dT / L;

    IIR_FILTER *iir = laplace_nd_filter(snum, nsnum, sden, nsden, dt, L);

    print_polynom_nd(snum, nsnum, sden, nsden, 's', "Laplace H(s)");
    printf("dT: %f, L: %d, dt: %f\n", dT, L, dt);
    print_polynom_nd(iir->b, iir->m + 1, iir->a, iir->n + 1, 'z', "IIR H(z)");

    printf("t\tx\ty\n");
    double t, x = 0, y;
    int l = 0;
    for (t = 0; t < tmax; t += dt) {
	if (--l <= 0) {
	    x = f(t);
	    l = L;
	}
	y = iir_filter_next(iir, x);
	printf("%f\t%f\t%f\n", t, x, y);
    }

    iir_filter_free(iir);

    return 0;
}
