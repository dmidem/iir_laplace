#ifndef POLYNOM_H
#define POLYNOM_H

typedef struct { double r; double i; } complex;

void clear_polynom(double *a, const int na);
void clear_polynom_c(complex *a, const int na);
void copy_polynom(const double *a, const int nb, double *b);
void copy_polynom_c(const complex *a, const int nb, complex *b);
void copy_polynom_cr(const complex *a, const int nb, double *b);
void copy_polynom_ci(const complex *a, const int nb, double *b);
void expand_polynom(const double *a, const int na, const double *b, const int nb, double *c);
void expand_polynom_c(const complex *a, const int na, const complex *b, const int nb, complex *c);
void expand_zerolist(const double *zl, const int nzl, const double zk, double *r);
void expand_zerolist_c(const complex *zl, const int nzl, const complex zk, complex *r);
void pow_polynom(const double *a, const int na, const int m, double *c);
void pow_polynom_c(const complex *a, const int na, const int m, complex *c);
void add_polynom(double *a, const int na, const double k, const double *b);
void add_polynom_c(complex *a, const int na, const double k, const complex *b);
void add_polynom_cc(complex *a, const int na, const complex k, const complex *b);
void mul_polynom(double *a, const int na, const double k);
void mul_polynom_c(complex *a, const int na, const double k);
void mul_polynom_cc(complex *a, const int na, const complex k);

#endif // POLYNOM_H
