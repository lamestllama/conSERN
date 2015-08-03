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

void Connect(uint64_t node, uint32_t component, EdgeList* edges, uint32_t* adjacency_list, uint64_t* node_offsets, uint32_t* node_degree);

int MarkConnectedComponents(int32_t N,  EdgeList* edges);

void MakeConnected(int32_t N,  EdgeList *edges, uint32_t BufferSize,
                   float* x, float* y, uint32_t component_count);
#endif /* defined(__conSERN__connected__) */
