//
//  FastSERN.h
//  MexFastSERN
//
//  Created by Eric Parsonage on 3/30/15.
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
