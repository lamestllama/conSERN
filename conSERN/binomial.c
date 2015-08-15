//
//  binomial.c
//  NodePlacementTesting
//
//  Created by Eric Parsonage on 7/24/15.
//  Copyright 2015. All rights reserved.
//
//

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

#include "binomial.h"


static double GammaGetRandSmallShapeParameter( double shape );

double GetGammaRand( double shape )
{
    double u, v;
    double w;
    double x;
    double y;
    double z;
    double a = shape - 1;
    
    if (shape < 1.0) return GammaGetRandSmallShapeParameter(shape);
    if (shape == 1.0)
    {
        u = (double)rand() / (double) RAND_MAX;
        if (u == 0.0) return DBL_MAX;
        if (u == 1.0) return 0.0;
        return -log(u);
        
    }
    
    for ( ; ;)
    {
        u = (double)rand() / (double) RAND_MAX;
        v = (double)rand() / (double) RAND_MAX;
        w = u * (1.0 - u);
        y = sqrt( 3.0 * (shape - 0.25) / w) * (u - 0.5);
        x = y + a;
        if ( v == 0.0 ) return x;
        w *= 4.0;
        v *= w;
        z = w * v * v;
        if (log(z) <= 2.0 * (a * log(x/a) - y) ) return x;
    }
}


static  double GammaGetRandSmallShapeParameter( double shape )
{
    double u;
    
    if (shape >= 1.0) return GetGammaRand(shape);
    u = (double)rand() / (double) RAND_MAX;
    return (u == 0.0) ? 0.0 : GetGammaRand(shape+1.0)*pow(u, 1.0/shape);
}


double GetBetaRnd( double a, double b)
{
    double ga = GetGammaRand(a);
    double gb = GetGammaRand(b);
    
    return ga / (ga + gb);
}



/* Binomial algorithm by Knuth */
uint32_t  GetBinomialRand ( double p, uint32_t n)
{
    uint32_t i, a, b, k = 0;
    
    while (n > 25)      // while this is quicker than direct calculation
    {
        double x;
        a = 1 + (n / 2);
        b = 1 + n - a;
        /* beta function */
        x = GetBetaRnd(a, b);
        
        if (x >= p)
        {
            n = a - 1;
            p /= x;
        }
        else
        {
            k += a;
            n = b - 1;
            p = (p - x) / (1 - x);
        }
    }
    
    for (i = 0; i < n; i++)  // direct calculation
    {
        double u = (double)rand() / (double) RAND_MAX;
        if (u < p) k++;
    }
    
    return k;
}
