#include <cmath>

namespace Statistics {

// Computes Pearson's r (correlation coefficient) between two variables x and y.
//
// Algorithm from "Numerical Recipes in C++", second edition, p643.
//
double correlation_r(double* x, double* y, int n)
{
    const double TINY = 1.0e-20;
    int j;
    double yt, xt;
    double syy = 0, sxy = 0, sxx = 0, ay = 0, ax = 0;

    // Compute averages
    for (j = 0; j < n; ++j) {
	ax += x[j];
	ay += y[j];
    }

    ax /= (double)n;
    ay /= (double)n;

    // Compute R
    for (j = 0; j < n; ++j) {
	xt = x[j] - ax;
	yt = y[j] - ay;
	sxx += xt * xt;
	syy += yt * yt;
	sxy += xt * yt;
    }

    return sxy / (sqrt(sxx * syy) + TINY);
}

// Given two sets of data x and y - which we suspect to have a linear correlation - 
// compute the coefficients a & b for the line y = a + bx, using least-squares
// as the fitting condition. We do not have weighting data, all points are considered equal.
//
// a & b are output parameters. No return value as the algorithm cannot fail, but as always,
// GIGO applies.
//
// Algorithm from "Numerical Recipes in C++", second edition, p670.
//
void linear_fit(double* x, double* y, int n, double& a, double& b)
{
    int i;
    double t, sxoss, ss;
    double sx = 0, sy = 0, st2 = 0;

    b = 0.0;

    // Accumulate sums
    for (i = 0; i < n; ++i) {
	sx += x[i];
	sy += y[i];
    }
    ss = (double)n;

    sxoss = sx / ss;

    for (i = 0; i < n; ++i) {
	t = x[i] - sxoss;
	st2 += t * t;
	b += t * y[i];
    }

    b /= st2;

    a = (sy - (sx * b)) / ss;

}

}
