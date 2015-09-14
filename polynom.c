#include <string.h>
#include <stdlib.h>

#include "polynom.h"

void clear_polynom(double *a, const int na)
{
    memset(a, 0, na * sizeof(a[0]));
}

void clear_polynom_c(complex *a, const int na)
{
    memset(a, 0, na * sizeof(a[0]));
}

void copy_polynom(const double *a, const int nb, double *b)
{
    memcpy(b, a, nb * sizeof(b[0]));
}

void copy_polynom_c(const complex *a, const int nb, complex *b)
{
    memcpy(b, a, nb * sizeof(b[0]));
}

void copy_polynom_cr(const complex *a, const int nb, double *b)
{
    int i;
    for (i = 0; i < nb; ++i)
	b[i] = a[i].r;
}

void copy_polynom_ci(const complex *a, const int nb, double *b)
{
    int i;
    for (i = 0; i < nb; ++i)
	b[i] = a[i].i;
}

void expand_polynom(const double *a, const int na, const double *b, const int nb, double *c)
{
    int i, j, nc = na + nb - 1;

    if (nc > 0)
	clear_polynom(c, nc);
    for (i = 0; i < na; ++i)
	for (j = 0; j < nb; ++j)
	    c[i + j] += a[i] * b[j];
}

void expand_polynom_c(const complex *a, const int na, const complex *b, const int nb, complex *c)
{
    int i, j, nc = na + nb - 1;

    if (nc > 0)
	clear_polynom_c(c, nc);
    for (i = 0; i < na; ++i)
	for (j = 0; j < nb; ++j) {
	    c[i + j].r += a[i].r * b[j].r - a[i].i * b[j].i;
	    c[i + j].i += a[i].r * b[j].i + a[i].i * b[j].r;
	}
}

void expand_zerolist(const double *zl, const int nzl, const double zk, double *r)
{
    double ri[nzl + 1];
    double zi[2] = { 0, 1 };
    int i;

    clear_polynom(r, nzl + 1);
    r[0] = zk;
    for (i = 0; i < nzl; ++i) {
        zi[0] = zl[i];
        copy_polynom(r, i + 1, ri);
        expand_polynom(ri, i + 1, zi, 2, r);
    }
}

void expand_zerolist_c(const complex *zl, const int nzl, const complex zk, complex *r)
{
    complex ri[nzl + 1];
    complex zi[2] = { { 0, 0 }, { 1, 0 } };
    int i;

    clear_polynom_c(r, nzl + 1);
    r[0] = zk;
    for (i = 0; i < nzl; ++i) {
        zi[0] = zl[i];
        copy_polynom_c(r, i + 1, ri);
        expand_polynom_c(ri, i + 1, zi, 2, r);
    }
}

void pow_polynom(const double *a, const int na, const int m, double *c)
{
    double b[na * m];
    int i, nb = na;

    if (m < 0)
	return;

    if (m == 0) {
	c[0] = 1;
	return;
    }

    copy_polynom(a, nb, c);
    for (i = 1; i < m; ++i) {
	copy_polynom(c, nb, b);
	expand_polynom(a, na, b, nb, c);
	nb += na;
    }
}

void pow_polynom_c(const complex *a, const int na, const int m, complex *c)
{
    complex b[na * m];
    int i, nb = na;

    if (m < 0)
	return;

    if (m == 0) {
	c[0].r = 1;
	c[0].i = 0;
	return;
    }

    copy_polynom_c(a, nb, c);
    for (i = 1; i < m; ++i) {
	copy_polynom_c(c, nb, b);
	expand_polynom_c(a, na, b, nb, c);
	nb += na;
    }
}

void add_polynom(double *a, const int na, const double k, const double *b)
{
    int i;
    for (i = 0; i < na; ++i)
	a[i] += b[i] * k;
}

void add_polynom_c(complex *a, const int na, const double k, const complex *b)
{
    int i;
    for (i = 0; i < na; ++i) {
	a[i].r += b[i].r * k;
	a[i].i += b[i].i * k;
    }
}

void add_polynom_cc(complex *a, const int na, const complex k, const complex *b)
{
    int i;
    for (i = 0; i < na; ++i) {
	a[i].r += b[i].r * k.r - b[i].i * k.i;
	a[i].i += b[i].r * k.i + b[i].i * k.r;
    }
}

void mul_polynom(double *a, const int na, const double k)
{
    int i;
    for (i = 0; i < na; ++i)
	a[i] = a[i] * k;
}

void mul_polynom_c(complex *a, const int na, const double k)
{
    int i;
    for (i = 0; i < na; ++i) {
	a[i].r = a[i].r * k;
	a[i].i = a[i].i * k;
    }
}

void mul_polynom_cc(complex *a, const int na, const complex k)
{
    int i;
    for (i = 0; i < na; ++i) {
	a[i].r = a[i].r * k.r - a[i].i * k.i;
	a[i].i = a[i].r * k.i + a[i].i * k.r;
    }
}
