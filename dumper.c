#include <stdio.h>

#include "dumper.h"

void print_polynom(const double *p, const int np, const char v)
{
    int i, j = 0;
    for (i = 0; i < np; ++i) {
	double pi = p[i];
	if (pi == 0)
	    continue;
	if (j == 0) {
	    if (pi < 0)
		printf("- ");
	} else {
	    if (pi < 0)
		printf(" - ");
	    else
		printf(" + ");
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

void print_polynom_nd(const double *num, const int nnum, const double *den, const int nden, const char v, const char *title)
{
    if (title)
	printf("\n%s:\n", title);
    print_polynom(num, nnum, v);
    printf("----------\n");
    print_polynom(den, nden, v);
    printf("\n");
}
