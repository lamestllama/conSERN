//
//  uniform.c
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#include "uniform.h"
#include <inttypes.h>
#include <stdlib.h>

uint32_t *uint_m_z;
uint32_t *uint_m_w;


void FreeUint()
{
    free(uint_m_z);
    free(uint_m_w);
};
