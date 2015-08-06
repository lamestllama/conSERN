//
//  metrics.c
//  conSERN
//
//  Created by Eric Parsonage on 8/4/15.
//
//

#include "metrics.h"
#include <math.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

double euclidean(double xdiff, double ydiff)
{
    return sqrt(xdiff * xdiff + ydiff * ydiff);
}
double manhattan(double xdiff, double ydiff)
{
    return xdiff + ydiff;
}

double maxdist(double xdiff, double ydiff)
{
    return max(xdiff, ydiff);
}

double mindist(double xdiff, double ydiff)
{
    return min(xdiff, ydiff);
}
