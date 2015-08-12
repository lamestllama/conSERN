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

typedef double (*ProbabilityFunction)(double p1, double p2, double q, double d);


extern double waxman(double s, double unused, double q, double d);


extern double clipped_waxman(double s, double r, double q, double d);


extern double waxman_transition_threshold(double s, double h, double q, double d);


extern double threshold(double r, double unused, double q, double d);


extern double constant(double unused1 , double unused2, double q, double d);


extern double powerlaw(double theta1, double theta2, double q, double d);


extern double cauchy(double theta1, double unused, double q, double d);


extern double exponential(double L, double unused, double q, double d);


extern double maxentropy(double s, double unused, double q, double d);



static ProbabilityFunction  probabilityFunctions[] =
    {   waxman,
        clipped_waxman,
        waxman_transition_threshold,
        threshold,
        constant,
        powerlaw,
        cauchy,
        exponential,
        maxentropy
    };


#endif /* defined(__conSERN__edgeprobfuncs__) */
