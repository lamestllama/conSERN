//
//  edgeprobfuncs.c
//  conSERN
//
//  Created by Eric Parsonage on 8/4/15.
//
//


#include "edgeprobfuncs.h"

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


double waxman(double s, __attribute__((unused))double unused, double d)
{
    return exp(-s * d);
}

double clipped_waxman(double s, double r, double d)
{
    return (d > r) ? 0.0 : exp(-s * d);
}

double constant(__attribute__((unused))double unused1 ,
                __attribute__((unused))double unused2,
                __attribute__((unused))double d)
{
    return 1.0;
}

double threshold(double r, __attribute__((unused))double unused,
                 __attribute__((unused))double d)
{
    return (d < r) ? 1.0 : 0.0;
}

double powerlaw(double theta1, double theta2, double d)
{
    return pow(1 + theta1 * d, -theta2);
}

double cauchy(double theta1, double __attribute__((unused))unused, double d)
{
    return pow(1 + theta1 * d * d, -1);
}

double exponential(double L, __attribute__((unused))double unused, double d)
{
    return exp(-d / (L - d));
}

double maxentropy(double s, double q2, double d)
{
   // to do test this against a version with a
   // temporary variable for speed
    return (q2 * exp(-s * d))/ (1 + q2 * exp(-s * d));
}

double clipped_waxman2(double s, double h, double d)
{
    return  min(exp(-s * d),h);
}

