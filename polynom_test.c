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
}
