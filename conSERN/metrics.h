//
//  metrics.h
//  conSERN
//
//  Created by Eric Parsonage on 8/4/15.
//
//

#ifndef __conSERN__metrics__
#define __conSERN__metrics__

#include <math.h>

typedef double (*DistanceFunction)(double xdiff, double ydiff);


extern double euclidean(double xdiff, double ydiff);
extern double manhattan(double xdiff, double ydiff);
extern double max(double xdiff, double ydiff);


static DistanceFunction  metrics[] =
{   euclidean,
    manhattan,
    max
};

#endif /* defined(__conSERN__metrics__) */
