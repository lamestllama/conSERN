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
    uint32_t *componentRoots;
    int8_t *signs;
    
    int64_t i, j, d;
    uint32_t c;
    uint32_t N;
    uint32_t sLargest;
    uint32_t iLargest;
    double distance;
    EdgeList edge_buffer;
    
    // shorthand
    N = options->N;
    
    // if we want to return the  components
    // to the caller then use memory allocated with
    // the callers memory allocation routines
    if (options->components_enabled)
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
    
    
    // compress the paths, count the components
    // and find the largest component
    sLargest = 0;
    for(i = 0, c = 0; i < N; i++)
    {
        find(i, roots, signs);
        if (signs[i] < 0)
        {
            c++;
            if(roots[i] > sLargest)
            {
                sLargest = roots[i];
                iLargest = i;
            }
            
        }
    
    }
    
    
    // simple form of assigning component labels can be
    // used if only the connected components and not a
    // connected graph are required
    if (!options->connected)
    {
         // give the largest components root label 1
        roots[iLargest] = 1;
        for (i = 0, c = 2; i < N; i++)
        {
            if (i == iLargest) continue;
            // label the remaining roots arbitrarly
            if (signs[i] < 0) roots[i] = c++;
        }
        c--;
        
        
        // labels the rest of the nodes using
        // the labels assigned their roots
        for (i = 0; i < N; i++)
            if (signs[i] > 0)
                roots[i] = roots[roots[i]];
        
        
        free(signs);
        
        return c;
    }
    
    
    // More complex form of assigning component
    // labels keeps track of the location of
    // component roots to speed subsequent itteration
    
    componentRoots = (uint32_t *) calloc(c, sizeof(uint32_t));
    
    // largest component assigned label 1 the
    // others are labeled arbitrarily labeled
    roots[iLargest] = 1;
    componentRoots[0] = iLargest;
    for (i = 0, c = 2; i < N; i++)
    {
        if (signs[i] < 0 && i != iLargest)
        {
           componentRoots[c - 1] = i;
           roots[i] = c++;
        }
    }
    c--;
    
    
    // now labels the rest of the nodes
    // using the labels of their roots
    for (i = 0; i < N; i++)
        if (signs[i] > 0)
            roots[i] = roots[roots[i]];
    
    // we will only add #components - 1 edges to connect
    // the graph now we can allocate the right sized buffer
    // prepare for allocation
    AllocateEdgeBuffer(&edge_buffer, c, edges->weights_enabled);
    
    // this is dirty and uses the fact that nodes are numbered
    // close together so we will mainly be adding short links
    // in order to make the graph a connected graph the more
    // buckets we have the shorter the links will be so not ideal
    for (i = 1; i < c; i++)
    {
        
        d = (componentRoots[i] < iLargest) ? 1 : -1;
        
        // look in  the direction of the largest components root
        // for a nearby node that belongs to the largest component.
        for (j = componentRoots[i]  + d;
             abs((int)(j - componentRoots[i])) <= LOCALE; j = j + d)
            if (roots[j] == 1) break;
        
        // if such a node is found within LOCALE steps connect to
        // it otherwise connect to the root of the largest component
        j = (roots[j] == 1) ? j : iLargest;
        
        // use the correct distance function
        distance = options->distance(nodes->x[componentRoots[i]]-nodes->x[j],
                                     nodes->y[componentRoots[i]]-nodes->y[j]);
        
        AddEdgeToBuffer(options, edges, &edge_buffer,
                        (uint32_t)componentRoots[i], (uint32_t)j, distance);
    }
    
    // if we are not returning the components
    // to the caller then free the roots
    if (!options->components_enabled) free(roots);
    
    free(componentRoots);
    free(signs);
    
    CopyEdgeList(options, edges, &edge_buffer);
    
    // now we have finished with the edge buffer
    // TODO write a free edgebuffer function 
    free(edge_buffer.from);
    free(edge_buffer.to);
    if (edges->weights_enabled)
        free(edge_buffer.weight);

    return c;
    
}




