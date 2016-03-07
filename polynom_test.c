/*
    Main module for polynom operations testing.
*/

#include <stdio.h>

#include "polynom.h"
#include "dumper.h"

int main()
{
    double a[] = { 1, 2, 3, 4 };
    int na = sizeof(a) / sizeof(a[0]);
    double b[] = { 5, 6, 7 };
    int nb = sizeof(b) / sizeof(b[0]);
    int nc = 20;
    double c[nc];

    double zl[] = { 1, 2, 3, 4 };
    int nzl = sizeof(zl) / sizeof(zl[0]);
    double zk = 1;
    int nr = nzl + 1;
    double r[nr];

    complex szero[] = { { 1, 2 }, { 3, 4 } };
    complex spole[] = { { 5, 6 }, { 7, 8 } };
    complex szerok = { 9, 0 };
    complex spolek = { 1, 0 };
    int nszero = sizeof(szero) / sizeof(szero[0]);
    int nspole = sizeof(spole) / sizeof(spole[0]);
    int nrzero = nszero + 1;
    complex rzero[nrzero];
    int nrpole = nspole + 1;
    complex rpole[nrpole];

    printf("a (%d):\n", na);
    print_polynom(a, na, 'x');
    printf("b (%d):\n", nb);
    print_polynom(b, nb, 'x');

    clear_polynom(c, nc);
    expand_polynom(a, na, b, nb, c);
    printf("c (%d):\n", nc);
    print_polynom(c, nc, 'x');

    clear_polynom(c, nc);
    pow_polynom(a, na, 2, c);
    printf("c (%d):\n", nc);
    print_polynom(c, nc, 'x');

    expand_zerolist(zl, nzl, zk, r);
    printf("r (%d):\n", nr);
    print_polynom(r, nr, 'x');

    expand_zerolist_c(szero, nszero, szerok, rzero);
    expand_zerolist_c(spole, nspole, spolek, rpole);
    print_polynom_nd_c(rzero, nrzero, rpole, nrpole, 's', "zp nd");

    return 0;
}
