//
//  connected.c
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#include "connected.h"
#include "fastSERN.h"
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>



// Union Set algorithm
// Returns the root node of x
int64_t find(int64_t x, int64_t *roots)
{
    // find the root
    int64_t root = x;
    while (roots[root] >= 0) root = roots[root];
    
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
bool merge(int64_t x, int64_t y, int64_t *roots)
{
    int64_t tmp;
    x = find(x, roots);
    y = find(y, roots);
    if (x == y) return false;
    // roots[x] is the negated size of the set x is in
    // make sure x is the larger set
    if (roots[x] > roots[y]) {tmp = x; x = y; y = tmp;}
    roots[x] += roots[y];
    // move the y tree under x
    roots[y] =  x;
    return true ;
}



int Components(int32_t N,  EdgeList* edges)
{
   
    int64_t *roots;
    int64_t i;
    uint32_t c;
    
    roots = (int64_t *) calloc(N, sizeof(int64_t));
    
    for(i = 0; i < N; i++) roots[i] = -1;
    
    for(i = 0; i < edges->count; i++)
        merge(edges->from[i] - 1, edges->to[i] - 1, roots);
    
    for(i = 0; i < N; i++) find(i, roots);
   
    for (i = 0, c = 1; i < N; i++)
        if (roots[i] < 0)
            edges->component[i] = c++;
    
    for (i = 0; i < N; i++)
        if (roots[i] > 0)
            edges->component[i] = edges->component[roots[i]];
    
    free(roots);
    return c - 1;
}




void MakeConnected(int32_t N,  EdgeList *edges, uint32_t BufferSize,
                   float* x, float* y, uint32_t component_count)
{
    
    uint64_t i, j;
    
    uint32_t *component_sizes;
    uint32_t *massive_component;
    uint32_t massive_size;
    uint32_t massive;
    double  distance = 0;
    EdgeList edge_buffer = {NULL, NULL, NULL, 0, 0, 0, 0};
    
    
    // TODO memory allocation error handling
    // components are numbered from 1
    // but we store there sizes in an array
    // indexed from zero
    component_sizes = calloc(component_count, sizeof(uint32_t));
    
    // the largest component is initially
    // set to the first component
    massive = 0;
    // for each node
    for(i = 0; i < N; i++)
    {
        // increment the component size
        // of the component that this
        // node belongs too
        component_sizes[edges->component[i]-1] += 1;
        // update which component is the largest
        if (component_sizes[edges->component[i]-1] > component_sizes[massive])
            massive = edges->component[i]-1;
        
    }
    
    
    // we know what component is the largest one and
    // how much space we need to accomodate it
    // TODO memory allocation error handling
    massive_component = calloc(component_sizes[massive], sizeof(uint32_t));
    
    // store the nodes that make up the largest component
    for (i = j = 0; i < N; i++)
    {
        if (edges->component[i]-1 == massive)
        {
            massive_component[j] = (uint32_t)i;
            j++;
            
        }continue;
    }
    massive_size = component_sizes[massive];
    
    // almost randomly choose the node to connect for each
    // component not in the largest component Note: slight bias to
    // smaller values because of the use of mod TODO fix randomness
    
    for(i = 0; i < component_count; i++)
    {
        component_sizes[i] = (rand() % component_sizes[i]) + 1;
    }
    
    // ensure nodes in largest component are ignored;
    component_sizes[massive] = 0;
    
    AllocateEdgeBuffer(&edge_buffer, BufferSize, edges->weights_enabled);
    
    // we should only need to add component - 1 more edges
    // for every node
    for(i = j= 0; i < N; i++)
    {
        switch(component_sizes[edges->component[i] - 1])
        {
                // the component that this node part of
                // has already been connected to the largest component
            case 0: break;
                
            case 1:
                // connect the component that this node is part of
                // to the largest component by selecting  random
                // node in the largest component
                j = rand() % massive_size;
                if(edges->weights_enabled)
                {
                    double xdiff = x[i] - x[massive_component[j]];
                    double ydiff = y[i] - y[massive_component[j]];
                    
                    distance = sqrt((xdiff * xdiff) + (ydiff * ydiff));
                }
                AddEdgeToBuffer(edges, &edge_buffer, (uint32_t)i, massive_component[j], distance);
                // case 1 purposely drops through to default processing
            default: component_sizes[edges->component[i] - 1] -= 1;
                
        }
    }
    
    free(component_sizes);
    free(massive_component);
    
    // flush out our buffer
    CopyEdgeList(edges, &edge_buffer);
}
