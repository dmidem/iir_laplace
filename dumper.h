/*
    Routines to print polynoms and filter transform functions (used for testing purposes).
*/

#ifndef DUMPER_H
#define DUMPER_H

#include "polynom.h"

void print_polynom(const double *p, const int np, const char v);
void print_polynom_c(const complex *p, const int np, const char v);
void print_polynom_nd(const double *num, const int nnum, const double *den, const int nden,
    const char v, const char *title);
void print_polynom_nd_c(const complex *num, const int nnum, const complex *den, const int nden,
    const char v, const char *title);
void print_zerolist(const double *z, const int nz, const char v);
void print_zerolist_c(const complex *z, const int nz, const char v);
void print_polynom_zp(const double *zero, const int nzero, const double *pole, const int npole,
    const double k, const char v, const char *title);
void print_polynom_zp_c(const complex *zero, const int nzero, const complex *pole, const int npole,
    const complex k, const char v, const char *title);

#endif // DUMPER_H
