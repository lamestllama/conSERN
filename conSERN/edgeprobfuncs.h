//
//  edgeprobfuncs.h
//  conSERN
//
//  Created by Eric Parsonage on 8/4/15.
//
//

#ifndef __conSERN__edgeprobfuncs__
#define __conSERN__edgeprobfuncs__


#include <math.h>

typedef double (*ProbabilityFunction)(double p1, double p2, double d);


extern double waxman(double s, double unused, double d);


extern double clipped_waxman(double s, double r, double d);


extern double constant(double unused1 , double unused2, double d);


extern double threshold(double r, double unused, double d);


extern double powerlaw(double theta1, double theta2, double d);


extern double cauchy(double theta1, double unused, double d);


extern double exponential(double L, double unused, double d);


extern double maxentropy(double s, double unused, double d);


static ProbabilityFunction  distfunc[] =
    {   waxman,
        clipped_waxman,
        constant,
        threshold,
        powerlaw,
        cauchy,
        exponential,
        maxentropy
    };


#endif /* defined(__conSERN__edgeprobfuncs__) */
