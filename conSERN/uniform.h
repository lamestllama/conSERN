//
//  uniform.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#ifndef __conSERN__uniform__
#define __conSERN__uniform__

#include <stdlib.h>
#include <inttypes.h>

// NOTE: this stuff is in a header so it can remain inline
// do not move it to a seperate compilation unit

#define INTP_TO_DOUBLEP 2.328306435454494e-10

extern  uint32_t *uint_m_z;
extern  uint32_t *uint_m_w;

extern void FreeRandom();


static inline void AllocRandom(uint32_t u, uint32_t ThreadCount)
{
    
    uint32_t i;
    
    srand(u);
    uint_m_z = malloc(sizeof(uint32_t) * ThreadCount);
    uint_m_w = malloc(sizeof(uint32_t) * ThreadCount);
    
    for(i = 0; i < ThreadCount; i++)
    {
        uint_m_z[i] = rand();
        uint_m_w[i] = rand();
    }
    
}

/* inline to avoid overhead of function call */
static inline uint32_t GetUint(uint32_t thread)
{
    uint32_t t;
    uint_m_z[thread] =
        36969 * (uint_m_z[thread] & 65535) + (uint_m_z[thread] >> 16);
    uint_m_w[thread] =
        18000 * (uint_m_w[thread] & 65535) + (uint_m_w[thread] >> 16);
    t = (uint_m_z[thread] << 16) + uint_m_w[thread];
    return t;
}

/* inline to avoid overhead of function call */
static inline double GetUniform(uint32_t thread)
{
    /* // 0 <= u < 2^32 */
    uint32_t u = GetUint(thread);
    /* // The magic number below is 1/(2^32 + 2). */
    /* // The result is strictly between 0 and 1. */
    return (u + 1.0) * INTP_TO_DOUBLEP;
}



#endif /* defined(__conSERN__uniform__) */
