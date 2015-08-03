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


// This is a dog that uses far too much stack and will crash
void Connect(uint64_t node, uint32_t component, EdgeList* edges, uint32_t* adjacency_list, uint64_t* node_offsets, uint32_t* node_degree)
{
    uint64_t i;
    
    // asign this node a component
    edges->component[node] = component;
    
    // for every node that this node is connected to
    for(i = node_degree[node]; i > 0; i--)
    {
        // if it isnt already assigned a component call this method on it
        if (edges->component[adjacency_list[node_offsets[node] - i]] == 0)
            Connect( adjacency_list[node_offsets[node] - i],
                    component, edges, adjacency_list,
                    node_offsets, node_degree);
    }
    
}


int MarkConnectedComponents(int32_t N,  EdgeList* edges)
{
    int64_t i;
    uint32_t component;
    uint32_t *node_degree;
    uint64_t *node_offsets;
    uint32_t *adjacency_list;
    
    // TODO handle memory allocation error
    node_degree = calloc(N, sizeof(uint32_t));
    if (node_degree==NULL) exit (1);
    
    // calculate the degree of each node.
    // the edges are assuming the nodes are
    // indexed from 1
    for (i = 0; i < edges->count; i++)
    {
        node_degree[edges->from[i] - 1] += 1;
        node_degree[edges->to[i] - 1] += 1;
    }
    
    // TODO handle memory allocation error
    node_offsets = calloc(N, sizeof(uint64_t));
    if (node_offsets==NULL) exit (2);
    // calculate the index of each node in an
    // array representing an adjacencly list;
    // the index stored is actually the begining
    // of the next nodes adjacency list.
    // i.e
    // for a square (4 nodes 4 edges)
    // node_degree = [2 2 2 2]
    // node_offset = [2 4 6 8]
    
    node_offsets[0] = node_degree[0];
    for (i = 1; i < N; i++)
        node_offsets[i] = node_degree[i] + node_offsets[i - 1];
    
    // TODO handle memory allocation error
    adjacency_list = calloc(node_offsets[N - 1], sizeof(uint32_t));
    
    // fill in the adjacency list
    for (i = 0; i < edges->count; i++)
    {
        // update the from node in the adjacency list
        // with the number of the too node
        // note nodes in the adjacency list are
        // indexed from zero
        adjacency_list[node_offsets[edges->from[i] - 1]
                       - node_degree[edges->from[i] - 1]]
        = edges->to[i] - 1;
        
        // using node_degree as a count of how
        // many edges left to store for a node
        node_degree[edges->from[i]-1] -= 1;
        
        
        adjacency_list[node_offsets[edges->to[i] - 1]
                       - node_degree[edges->to[i] - 1]]
        = edges->from[i] - 1;
        
        node_degree[edges->to[i] - 1] -= 1;
    }
    
    
    // so using are above example assuming clockwise from top left
    // numbering adjacency_list = [2 4 1 3 2 4 3 1]
    
    // we are going to need node_degrees again so we rebuild them
    // from the node_offsets;
    node_degree[0] = (uint32_t)node_offsets[0];
    for (i = 1; i < N; i++)
        node_degree[i] = (uint32_t)(node_offsets[i] - node_offsets[i - 1]);
    
    
    // now we are ready for a standard connected components algorithm
    component = 0;
    for(i = 0; i < N; i++)
    {
        // this node has already been visited
        if (edges->component[i] > 0) continue;
        // a new component has been discovered process it
        component++;
        edges->component[i] = component;
        // TODO write an iterative version of this
        Connect(i, component, edges, adjacency_list, node_offsets, node_degree);
    }
    
    free(adjacency_list);
    free(node_offsets);
    free(node_degree);
    
    return component;
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
