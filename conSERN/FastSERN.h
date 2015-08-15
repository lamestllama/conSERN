//
//  FastSERN.h
//  conSERN
//
//  Created by Eric Parsonage on 3/30/15.
//  Copyright 2015. All rights reserved.
//
//

#ifndef __FastSERN__
#define __FastSERN__

#include <stdint.h>
#include "edgelist.h"
#include "options.h"
#include "nodegen.h"


extern int GenSERN(NodeList* nodes, EdgeList* edges, Options* options, GeometryStruct *geometry);


#endif /* defined(__FastSERN__) */
