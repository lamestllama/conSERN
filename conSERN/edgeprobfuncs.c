//
//  edgeprobfuncs.c
//  conSERN
//
//  Created by Eric Parsonage on 8/4/15.
//
//


#include "edgeprobfuncs.h"




double waxman(double s, double unused, double d)
{
    return exp(-s * d);
}

double clipped_waxman(double s, double r, double d)
{
    return (d > r) ? 0.0 : exp(-s * d);
}

double constant(double unused1 , double unused2, double d)
{
    return 1.0;
}

double threshold(double r, double unused, double d)
{
    return (d < r) ? 1.0 : 0.0;
}

double powerlaw(double theta1, double theta2, double d)
{
    return pow(1 + theta1 * d, -theta2);
}

double cauchy(double theta1, double unused, double d)
{
    return pow(1 + theta1 * d * d, -1);
}

double exponential(double L, double unused, double d)
{
    return exp(-d / (L - d));
}

double maxentropy(double s, double unused, double d)
{
   // to do test this against a version with a
   // temporary variable for speed
    return exp(-s * d)/ (1 + exp(-s * d));
}

