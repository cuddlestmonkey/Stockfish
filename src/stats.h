#ifndef SF_STATS_H
#define SF_STATS_H

namespace Statistics {
double correlation_r(double* x, double* y, int n);
void linear_fit(double* x, double* y, int n, double& a, double& b);
}

#endif // SF_STATS_H
