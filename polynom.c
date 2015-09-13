#include <string.h>
#include <stdlib.h>

void clear_polynom(double *a, const int na)
{
    memset(a, 0, na * sizeof(a[0]));
}

void copy_polynom(const double *a, const int nb, double *b)
{
    memcpy(b, a, nb * sizeof(b[0]));
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

void pow_polynom(const double *a, const int na, int m, double *c)
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


void add_polynom(double *a, const int na, const double k, double *b)
{
    int i;
    for (i = 0; i < na; ++i)
	a[i] += b[i] * k;
}


void mul_polynom(double *a, const int na, const double k)
{
    int i;
    for (i = 0; i < na; ++i)
	a[i] *= k;
}
