//
//  Header.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#ifndef conSERN_Header_h
#define conSERN_Header_h
#include "edgeprobfuncs.h"
#include "metrics.h"

typedef struct
{
    //uint32_t probabilityFunction;
    double s1;
    double s2;
    double q;
    uint32_t N;
    uint32_t M;
    //uint32_t distanceFunction;
    
    uint32_t connected;
    uint32_t algorithm;
    uint32_t ThreadCount;
    uint32_t BufferSize;
    
    
    uint32_t seedval;
    uint32_t components_enabled;
    uint32_t weights_enabled;
    ProbabilityFunction probGivenDistance;
    DistanceFunction distance;
   
    

} Options;

#endif
