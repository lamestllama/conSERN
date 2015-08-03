//
//  geometric.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#ifndef __conSERN__geometric__
#define __conSERN__geometric__

#include "uniform.h"
#include <math.h>


static inline uint32_t geom_rand(double p, uint32_t thread)
{
    /* should probably check that 0<=p<1 */
    return(1 + floor(log((GetUint(thread) + 1) * INTP_TO_DOUBLEP) / log(1 - p)));//
}



/* lambda = -log2(1-p) */
static inline uint32_t geom_rand2(double lambda, uint32_t thread)
{
    /* should probably check that 0<=p<1 */
    return(-floor((log2(GetUint(thread) + 1) - 32) / lambda));
}

#endif /* defined(__conSERN__geometric__) */
