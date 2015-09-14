#include <stdio.h>

#include "dumper.h"

void print_polynom(const double *p, const int np, const char v)
{
    int i, j = 0;

    for (i = 0; i < np; ++i) {
	double pi = p[i];
	if (pi == 0)
	    continue;
	if (j > 0) {
	    if (pi < 0)
		printf(" - ");
	    else
		printf(" + ");
	} else {
	    if (pi < 0)
		printf("- ");
	}
	if (pi < 0)
	    pi = -pi;
	printf("%f", pi);
	if (i > 0)
	    printf("%c", v);
	if (i > 1)
	    printf("^%d", i);
	++j;
    }

    printf("\n");
}

void print_polynom_c(const complex *p, const int np, const char v)
{
    int i, j = 0;

    for (i = 0; i < np; ++i) {
	complex pi = p[i];
	if ((pi.i == 0) && (pi.r == 0))
	    continue;
	if (j > 0)
	    printf(" + ");
	printf("[%f,%f]", pi.r, pi.i);
	if (i > 0)
	    printf("%c", v);
	if (i > 1)
    	    printf("^%d", i);
	++j;
    }

    printf("\n");
}

void print_polynom_nd(const double *num, const int nnum, const double *den, const int nden,
    const char v, const char *title)
{
    if (title)
	printf("\n%s:\n", title);
    print_polynom(num, nnum, v);
    printf("----------\n");
    print_polynom(den, nden, v);
    printf("\n");
}

void print_polynom_nd_c(const complex *num, const int nnum, const complex *den, const int nden,
    const char v, const char *title)
{
    if (title)
	printf("\n%s:\n", title);
    print_polynom_c(num, nnum, v);
    printf("----------\n");
    print_polynom_c(den, nden, v);
    printf("\n");
}

void print_zerolist(const double *z, const int nz, const char v)
{
    int i;

    for (i = 0; i < nz; ++i) {
	if (i > 0)
	    printf(" * ");
	printf("(%f + %c", z[i], v);
	if (i > 1)
    	    printf("^%d", i);
	printf(")");
    }

    printf("\n");
}

void print_zerolist_c(const complex *z, const int nz, const char v)
{
    int i;

    for (i = 0; i < nz; ++i) {
	if (i > 0)
	    printf(" * ");
	printf("([%f,%f] + %c", z[i].r, z[i].i, v);
	if (i > 1)
    	    printf("^%d", i);
	printf(")");
    }

    printf("\n");
}

void print_polynom_zp(const double *zero, const int nzero, const double *pole, const int npole,
    const double k, const char v, const char *title)
{
    if (title)
	printf("\n%s:\n", title);
    printf("%f *", k);
    print_zerolist(zero, nzero, v);
    printf("----------\n");
    print_zerolist(pole, npole, v);
    printf("\n");
}

void print_polynom_zp_c(const complex *zero, const int nzero, const complex *pole, const int npole,
    const complex k, const char v, const char *title)
{
    if (title)
	printf("\n%s:\n", title);
    printf("[%f,%f] *", k.r, k.i);
    print_zerolist_c(zero, nzero, v);
    printf("----------\n");
    print_zerolist_c(pole, npole, v);
    printf("\n");
}
