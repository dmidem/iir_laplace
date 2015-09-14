#ifndef FILTER_H
#define FILTER_H

#include "polynom.h"

typedef struct {
    double *b;
    int m;

    double *a;
    int n;

    double *x;
    double *y;
} IIR_FILTER;

void bilinear_transform_n(const double *s, const int ns, const double t, double *znum, double *zden);
int bilinear_transform_nd(const double *snum, const int nsnum, const double *sden, const int nsden,
    const double t, double *znum, double *zden);

IIR_FILTER *iir_filter_init(const double *b, const int m, const double *a, const int n);
void iir_filter_free(IIR_FILTER *iir);
double iir_filter_next(IIR_FILTER *iir, const double x);

IIR_FILTER *laplace_nd_filter(const double *num, const int nnum,
    const double *den, const int nden,
    const double dt, const int L);
IIR_FILTER *laplace_zp_filter(const complex *zlc, const int nzlc,
    const complex *plc, const int nplc,
    const complex zkc, const double dt, const int L);

#endif // FILTER_H
