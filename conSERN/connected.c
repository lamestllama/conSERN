//
//  connected.c
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#include "connected.h"
#include "options.h"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#define LOCALE 5



uint64_t find(uint32_t x, uint32_t *roots, int8_t *signs)
{
    // find the root
    int64_t root = x;
    while (signs[root] >= 0) root = roots[root];
    
    // compress paths
    while (x != root)
    {
        int64_t old = x;
        x = roots[x];
        roots[old] =  root ;
    }
    return root ;
}

// combines the sets that contain x and y
// false if x and y already in the same set
bool merge(uint32_t x, uint32_t y, uint32_t *roots, int8_t *signs)
{
    uint32_t tmp;
    x = find(x, roots, signs);
    y = find(y, roots, signs);
    if (x == y) return false;
    // roots[x] is the  size of the set x is in
    // make sure x is the larger set
    if (roots[x] < roots[y]) {tmp = x; x = y; y = tmp;}
    roots[x] += roots[y];
    // move the y tree under x
    roots[y] =  x;
    signs[y] =  1;
    return true ;
}


// Union Set algorithm to find the connected components and then
// a heuristic to link the components together in reasonable time

int Components(Options *options,  NodeList *nodes, EdgeList* edges)
{
    
    uint32_t *roots;
    int8_t *signs;
    
    int64_t i, j, d, e;
    uint32_t c;
    uint32_t N;
    uint32_t sLargest;
    uint32_t iLargest;
    double distance;
    EdgeList edge_buffer;
    
    // shorthand
    N = options->N;
    
    // if we want to return the connected components
    // to the caller then use memory allocated with
    // the callers memory allocation routines
    if (options->connected)
    {
        nodes->component = options->calloc(options->N, sizeof(uint32_t));
        roots = nodes->component;
    }
    else  // otherwise use the stdlib allocation routines
        roots = calloc(options->N, sizeof(uint32_t));
    
    // use this to indicate the roots of connected components
    signs = (int8_t *) calloc(N, sizeof(int8_t));
    
    // set all nodes to be their own root
    for(i = 0; i < N; i++)
    {
        roots[i] = 1;
        signs[i] = -1;  // -1 inicates this is a root
    }
    
    // merge into components
    for(i = 0; i < edges->count; i++)
        merge(edges->from[i] - 1, edges->to[i] - 1, roots, signs);
    
    
    // compress the paths and find the largest component
    sLargest = 0;
    for(i = 0; i < N; i++)
    {
        find(i, roots, signs);
        if (signs[i] < 0 && roots[i] > sLargest)
        {
            sLargest = roots[i];
            iLargest = i;
        }
    
    }

    // for convenience we make the largest
    // component have the label 1 the rest
    // are labeled arbitrarily
    // first we label the roots
    roots[iLargest] = 1;
    for (i = 0, c = 2; i < N; i++)
    {
        if (i == iLargest) continue;
        if (signs[i] < 0) roots[i] = c++;
    }
    c--;
    
    
    // now labels the rest of the nodes
    // using the labels of their roots
    for (i = 0; i < N; i++)
        if (signs[i] > 0)
            roots[i] = roots[roots[i]];
    
    // only components required
    if (!options->connected)
    {
        free(signs);
        return c;
    }
    
    // we will only add #components - 1 edges to connect
    // the graph now we can allocate the right sized buffer
    // prepare for allocation
    AllocateEdgeBuffer(&edge_buffer, c, edges->weights_enabled);
    
    // this is dirty and uses the fact that nodes are numbered
    // close together so we will mainly be adding short links
    // in order to make the graph a connected graph the more
    // buckets we have the shorter the links will be so not ideal
    for (i = 0, e = c - 1; i < N && e > 0; i++)
    {
        // if we find the root of a component that is not the largest one
        if (signs[i] < 0 &&  i != iLargest)
        {
            d = (i < iLargest) ? 1 : -1;
            --e; // keeping track of how many left to do
            
            // look at the next few in the direction of the largest
            // components root to see if a close node that belongs
            // to the largest component can be connected too otherwise
            // connect too the root of the largest component
            for (j = i + d; abs((int)(j - i)) <= LOCALE; j = j + d)
                if (roots[j] == 1)
                    break;
                    
            j = (roots[j] == 1) ? j : iLargest;
            
            // use the correct distance function
            distance = options->distance(nodes->x[i]-nodes->x[j],
                                         nodes->y[i]-nodes->y[j]);
            
            AddEdgeToBuffer(options, edges, &edge_buffer,
                            (uint32_t)i, (uint32_t)j, distance );
        }
        
    }
    
    if (!options->connected) free(roots);
    
    free(signs);
    CopyEdgeList(options, edges, &edge_buffer);
    
    free(edge_buffer.from);
    free(edge_buffer.to);
    
    if (edges->weights_enabled)
        free(edge_buffer.weight);

    return c;
    
}




