#include <string.h>
#include <stdlib.h>

#include "filter.h"


void bilinear_transform_n(const double *s, const int ns, const double t, double *znum, double *zden)
{
    const double z1num[2] = { 2 / t, - 2 / t };
    const double z1den[2] = { 1, 1 };

    double z[ns];
    int nz;

    double h = 1;
    int i;

    if (ns == 0)
	return;

    clear_polynom(znum, ns);
    clear_polynom(zden, ns);
    znum[0] = 1;
    zden[0] = 1;
    nz = 1;

    for (i = ns - 1; i > 0; --i) {
	h = s[i];
	if (h != 0)
	    break;
    }

    if (i == 0) {
	znum[0] = s[0];
	return;
    }

    for (i = i - 1; i >= 0; --i) {
	expand_polynom(znum, nz, z1num, 2, z);
	copy_polynom(z, nz + 1, znum);
	expand_polynom(zden, nz, z1den, 2, z);
	copy_polynom(z, nz + 1, zden);
	++nz;
	add_polynom(znum, nz, s[i] / h, zden);
    }
    mul_polynom(znum, ns, h);
}


int bilinear_transform_nd(const double *snum, const int nsnum, const double *sden, const int nsden,
    const double t, double *znum, double *zden)
{
    double znumnum[nsnum], znumden[nsnum], zdennum[nsden], zdenden[nsden];
    int nz = (nsnum > nsden ? nsnum : nsden);
    int nznumden, nzdenden;
    double l;

    const double z1[2] = { 1, 1 };

    bilinear_transform_n(snum, nsnum, t, znumnum, znumden);
    bilinear_transform_n(sden, nsden, t, zdennum, zdenden);

/*
    print_polynom_nd(znumnum, nsnum, znumden, nsnum, 'z', "zn");
    print_polynom_nd(zdennum, nsden, zdenden, nsden, 'z', "zd");
*/

    clear_polynom(znumden, nsnum);
    clear_polynom(zdenden, nsden);
    nzdenden = 0;
    nznumden = 0;
    if (nsden > 0) {
	zdenden[0] = 1;
	nzdenden = 1;
    }
    if (nsnum > 0) {
	znumden[0] = 1;
	nznumden = 1;
    }
    if (nsnum > nsden) {
	pow_polynom(z1, 2, nsnum - nsden, znumden);
	nznumden = (nsnum - nsden) + 1;
    }
    if (nsden > nsnum) {
	pow_polynom(z1, 2, nsden - nsnum, zdenden);
	nzdenden = (nsden - nsnum) + 1;
    }

    expand_polynom(znumnum, nsnum, zdenden, nzdenden, znum);
    expand_polynom(zdennum, nsden, znumden, nznumden, zden);

    if (nz > 0) {
	l = zden[0];
	if (l != 0) {
	    mul_polynom(zden, nz, 1 / l);
	    mul_polynom(znum, nz, 1 / l);
	}
    }

    return nz;
}


IIR_FILTER *iir_filter_init(const double *b, const int m, const double *a, const int n)
{
    IIR_FILTER *iir;

    if ((iir = (IIR_FILTER *)malloc(sizeof(IIR_FILTER))) == NULL)
	return NULL;

    memset(iir, 0, sizeof(IIR_FILTER));

    iir->b = (double *)calloc(m + 1, sizeof(double));
    iir->m = m;
    iir->a = (double *)calloc(n + 1, sizeof(double));
    iir->n = n;

    iir->x = (double *)calloc(m + 1, sizeof(double));
    iir->y = (double *)calloc(n + 1, sizeof(double));

    memcpy(iir->b, b, (m + 1) * sizeof(double));
    memcpy(iir->a, a, (n + 1) * sizeof(double));

    return iir;
}


void iir_filter_free(IIR_FILTER *iir)
{
    if (iir) {
	free(iir->b);
	free(iir->a);
	free(iir->x);
	free(iir->y);
    }
    free(iir);
}


double iir_filter_next(IIR_FILTER *iir, const double x)
{
    double sx = 0, sy = 0;
    double y;
    int i;

    for (i = iir->m; i > 0; --i)
        iir->x[i] = iir->x[i - 1];
    iir->x[0] = x;
    for (i = 0; i <= iir->m; ++i)
        sx += iir->b[i] * iir->x[i];

    for (i = iir->n; i > 0; --i)
        iir->y[i] = iir->y[i - 1];
    for (i = 1; i <= iir->n; ++i)
        sy += iir->a[i] * iir->y[i];
    y = sx - sy;
    iir->y[0] = y;

    return y;
}


IIR_FILTER *laplace_nd_filter(const double *num, const int nnum,
    const double *den, const int nden,
    const double dt, const int L)
{
    double interp_znum[] = { 0.0221824, -0.0134342, -0.0134342, 0.0221824 };
    int ninterp_znum = sizeof(interp_znum) / sizeof(interp_znum[0]);
    double interp_zden[] = { 1, -2.5885253, 2.3163762, -0.7103546 };
    int ninterp_zden = sizeof(interp_zden) / sizeof(interp_zden[0]);

    int nz = (nnum > nden ? nnum : nden);
    double laplace_znum[nz];
    double laplace_zden[nz];

    int nznum = nz + ninterp_znum - 1;
    int nzden = nz + ninterp_zden - 1;
    double znum[nznum];
    double zden[nzden];

    clear_polynom(laplace_znum, nz);
    clear_polynom(laplace_zden, nz);
    bilinear_transform_nd(num, nnum, den, nden, dt, laplace_znum, laplace_zden);

    if (L > 1) {
	expand_polynom(laplace_znum, nz, interp_znum, ninterp_znum, znum);
	expand_polynom(laplace_zden, nz, interp_zden, ninterp_zden, zden);
    } else {
	copy_polynom(laplace_znum, nz, znum);
	copy_polynom(laplace_zden, nz, zden);
	nznum = nz;
	nzden = nz;
    }

    IIR_FILTER *iir = iir_filter_init(znum, nznum - 1, zden, nzden - 1);

    return iir;
}


IIR_FILTER *laplace_zp_filter(const complex *zlc, const int nzlc,
    const complex *plc, const int nplc,
    const complex zkc, const double dt, const int L)
{
    int nnum = nzlc + 1;
    complex num_c[nnum];
    double num[nnum];

    complex pkc = { 1, 0 };
    int nden = nplc + 1;
    complex den_c[nden];
    double den[nden];

    expand_zerolist_c(zlc, nzlc, zkc, num_c);
    copy_polynom_cr(num_c, nnum, num);
    expand_zerolist_c(plc, nplc, pkc, den_c);
    copy_polynom_cr(den_c, nden, den);

    //print_polynom_nd(num, nnum, den, nden, 's', "zp nd");

    return laplace_nd_filter(num, nnum, den, nden, dt, L);
}
