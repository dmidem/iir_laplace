/*
    IIR filter-related functions.

*/

#ifndef FILTER_H
#define FILTER_H

#include "polynom.h"

// IIR-filter
typedef struct {
    // Feedforward filter coefficients
    double *b;
    int m;

    // Feedback filter coefficients
    double *a;
    int n;

    // Filter state (input and output signal values history)
    double *x;
    double *y;
} IIR_FILTER;

// Bilinear transform implementations (https://en.wikipedia.org/wiki/Bilinear_transform)
void bilinear_transform_n(const double *s, const int ns, const double t, double *znum, double *zden);
int bilinear_transform_nd(const double *snum, const int nsnum, const double *sden, const int nsden,
    const double t, double *znum, double *zden);

// Create IIR-filter based on given difference equation 
IIR_FILTER *iir_filter_init(const double *b, const int m, const double *a, const int n);

// Free (destroy) IIR-filter
void iir_filter_free(IIR_FILTER *iir);

// Get next output value of IIR-filter
double iir_filter_next(IIR_FILTER *iir, const double x);

// Zero-Denominator Laplace transform IIR-filter
IIR_FILTER *laplace_nd_filter(const double *num, const int nnum,
    const double *den, const int nden,
    const double dt, const int L);

// Zero-Pole Laplace transform IIR-filter
IIR_FILTER *laplace_zp_filter(const complex *zlc, const int nzlc,
    const complex *plc, const int nplc,
    const complex zkc, const double dt, const int L);

#endif // FILTER_H
