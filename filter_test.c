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
    double sden[] = { 1, 1 };
    int nsnum = sizeof(snum) / sizeof(snum[0]);
    int nsden = sizeof(sden) / sizeof(sden[0]);

    complex szero[] = { { 1, 2 }, { 3, 4 } };
    complex spole[] = { { 5, 6 }, { 7, 8 } };
    complex sk = { 9, 0 };
    int nszero = sizeof(szero) / sizeof(szero[0]);
    int nspole = sizeof(spole) / sizeof(spole[0]);

    double dT = 0.1;
    int L = 10;
    double tmax = 0.2;

    double dt = dT / L;

    IIR_FILTER *iir_nd = laplace_nd_filter(snum, nsnum, sden, nsden, dt, L);
    IIR_FILTER *iir_zp = laplace_zp_filter(szero, nszero, spole, nspole, sk, dt, L);

    printf("dT: %f, L: %d, dt: %f\n", dT, L, dt);

    print_polynom_nd(snum, nsnum, sden, nsden, 's', "Laplace ND H(s)");
    print_polynom_nd(iir_nd->b, iir_nd->m + 1, iir_nd->a, iir_nd->n + 1, 'z', "IIR H(z)");

    printf("t\tx\ty\n");
    double t, x = 0, y;
    int l = 0;
    for (t = 0; t < tmax; t += dt) {
	if (--l <= 0) {
	    x = f(t);
	    l = L;
	}
	y = iir_filter_next(iir_nd, x);
	printf("%f\t%f\t%f\n", t, x, y);
    }

    print_polynom_zp_c(szero, nszero, spole, nspole, sk, 's', "Laplace ZP H(s)");
    print_polynom_nd(iir_zp->b, iir_zp->m + 1, iir_zp->a, iir_zp->n + 1, 'z', "IIR H(z)");

    printf("t\tx\ty\n");
    x = 0;
    l = 0;
    for (t = 0; t < tmax; t += dt) {
	if (--l <= 0) {
	    x = f(t);
	    l = L;
	}
	y = iir_filter_next(iir_zp, x);
	printf("%f\t%f\t%f\n", t, x, y);
    }

    iir_filter_free(iir_nd);
    iir_filter_free(iir_zp);

    return 0;
}
