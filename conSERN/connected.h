//
//  connected.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#ifndef __conSERN__connected__
#define __conSERN__connected__

#include "edgelist.h"
#include "nodegen.h"

int Components(Options *options,  NodeList *nodes, EdgeList* edges);

void MakeConnected(int32_t N, NodeList *nodes, EdgeList *edges, uint32_t BufferSize,
                   float* x, float* y, uint32_t component_count);

#endif /* defined(__conSERN__connected__) */
