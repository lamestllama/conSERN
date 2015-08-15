//
//  FastSERN.c
//  conSERN
//
//  Created by Eric Parsonage on 3/30/15.
//  Copyright 2015. All rights reserved.
//
//

#include "FastSERN.h"
#include "binomial.h"
#include "uniform.h"
#include "geometric.h"
#include "edgelist.h"
#include "connected.h"
#include "nodegen.h"
#include "threads.h"
#include <assert.h>
#include <pthread.h>



void * BusyWork(void *t)
{
    
    int64_t i, j;
    int64_t k, S;
    uint32_t bucket_a, bucket_b;
    double p, lambda, distance, x_diff, y_diff;
    BucketStruct bucket_A, bucket_B;
    
    ThreadDataStruct *thread_data;
    uint32_t     thread_id;
    uint32_t     thread_count;
    uint32_t         Mx;
    uint32_t         My;
    double       s1, s2, q;
    
    EdgeList     *edges;
    BucketStruct *buckets;
    double       *Q;
    ProbabilityFunction prob;
    DistanceFunction dist;
    
    EdgeList edge_buffer;
    
    thread_data = (ThreadDataStruct *)t;
    
    // store all the frequently used values locally
    thread_id = thread_data->thread_id;
    thread_count = thread_data->thread_count;
    
    // geometry
    Mx = thread_data->geometry->Mx;
    My = thread_data->geometry->My;
    
    // functions for calculating probability and distance
    // based on chosen model
    prob = thread_data->options->probGivenDistance;
    dist = thread_data->options->distance;
    
    // shape parameters for model
    s1 = thread_data->options->s1;
    s2 = thread_data->options->s1;
    
    // thining parameter
    q = thread_data->options->q;
    
    // pointer to storage for the edges generated
    edges = thread_data->edges;
    
    // the set of buckets that cover the region our graph is on
    buckets = thread_data->buckets;
    
    // lookup table of maximum probabilities between buckets
    Q = thread_data->Q;
    
    
    // create some space to store the edges
    AllocateEdgeBuffer(thread_data->options, &edge_buffer,
                       thread_data->BufferSize,
                       edges->weights_enabled);
    
    if ((thread_data->options->algorithm == 0) && (1.0 - q > 0.0))
    {
        
        /* Generate intra bucket edges in bucket_a */
        // this is the fast algorithm
        lambda = -log2(1 - q);
        
        for (bucket_a = thread_id; bucket_a < Mx * My; bucket_a += thread_count)
        {
            
            k = -1;
            bucket_A = buckets[bucket_a];
            S = ((int64_t)bucket_A.count * ((int64_t)bucket_A.count - 1))/2;
            
            // find the next link that might exist
            // geom_rand2 returns an int64 in the range 0 ... INT64_MAX
            // loop terminates when either k is greater than the number
            // of edges possible between nodes in the bucket or k
            // has wrapped meaning it has skipped past the 2^63 possible
            // edges we support and is now negative.
            while ((k += geom_rand2(lambda, thread_id)) < S && k >= 0)
            {
                // can we do an integer square root here ?
                j = 1 + ((((int64_t)sqrtl(8 * k + 1)) - 1) / 2);
                i = k - j * (j - 1) / 2;
                
                x_diff = bucket_A.x[i] - bucket_A.x[j];
                y_diff = bucket_A.y[i] - bucket_A.y[j];
                
                distance = dist(x_diff, y_diff);
                
                if (GetUniform(thread_id) < prob(s1, s2, q, distance))
                    AddEdgeToBuffer(thread_data->options, edges, &edge_buffer,
                                    bucket_A.start + (uint32_t)i,
                                    bucket_A.start + (uint32_t)j,
                                    distance);
                
            }
            
        }
        
        
        // Generate inter bucket edges between bucket_a and bucket_b
        for (bucket_a = thread_id; bucket_a < Mx * My; bucket_a += thread_count)
        {
            bucket_A = buckets[bucket_a];
            
            for (bucket_b = bucket_a + 1; bucket_b < Mx * My; bucket_b++)
            {
                k = -1;
                bucket_B = buckets[bucket_b];
                
                p =  Q[abs((int32_t)bucket_A.i - (int32_t)bucket_B.i) +
                       (int32_t)Mx *
                       abs((int32_t)bucket_A.j - (int32_t)bucket_B.j)];
                
                lambda = -log2(1 - p * q);
                
                S = (int64_t)bucket_A.count * (int64_t)bucket_B.count;
                
                // find the next link that might exist
                while ((k += geom_rand2(lambda, thread_id)) < S && k >= 0)
                {
                    i = k % bucket_A.count;
                    j = k / bucket_A.count;
                    
                    x_diff = bucket_A.x[i] - bucket_B.x[j];
                    y_diff = bucket_A.y[i] - bucket_B.y[j];
                    
                    distance = dist(x_diff, y_diff);
                    
                    if (GetUniform(thread_id) * p < prob(s1, s2, q, distance))
                        AddEdgeToBuffer(thread_data->options, edges, &edge_buffer,
                                        bucket_A.start + (uint32_t)i,
                                        bucket_B.start + (uint32_t)j,
                                        distance);
                }
            }
        }
        
    }
    else
    {
        // for now any other option means use the N^2 algorithm
        // but at least in a parallel fashion
        /* Generate intra bucket edges in bucket_a */
        
        for (bucket_a = thread_id; bucket_a < Mx * My;
             bucket_a += thread_data->thread_count)
        {
            bucket_A = buckets[bucket_a];
            
            for (i = 0; i < bucket_A.count; i++)
            {
                for(j = i + 1; j < bucket_A.count; j++)
                {
                    x_diff = bucket_A.x[i] - bucket_A.x[j];
                    y_diff = bucket_A.y[i] - bucket_A.y[j];
                    
                    distance = dist(x_diff, y_diff);
                    
                    if (GetUniform(thread_id) < q * prob(s1, s2, q, distance))
                        AddEdgeToBuffer(thread_data->options, edges, &edge_buffer,
                                        bucket_A.start + (uint32_t)i,
                                        bucket_A.start + (uint32_t)j,
                                        distance);
                }
            }
        }
        
        
        // Generate inter bucket edges between bucket_a and bucket_b
        for (bucket_a = thread_id; bucket_a < Mx * My;
             bucket_a += thread_data->thread_count)
        {
            bucket_A = buckets[bucket_a];
            
            for (bucket_b = bucket_a + 1; bucket_b < Mx * My; bucket_b++)
            {
                bucket_B = buckets[bucket_b];
                for ( i = 0; i < bucket_A.count; i++)
                {
                    for(j = 0; j < bucket_B.count; j++)
                    {
                        
                        x_diff = bucket_A.x[i] - bucket_B.x[j];
                        y_diff = bucket_A.y[i] - bucket_B.y[j];
                        
                        distance = dist(x_diff, y_diff);
                        
                        if (GetUniform(thread_id) < q * prob(s1, s2, q, distance))
                            AddEdgeToBuffer(thread_data->options,
                                            edges, &edge_buffer,
                                            bucket_A.start + (uint32_t)i,
                                            bucket_B.start + (uint32_t)j,
                                            distance);
                    }
                }
            }
        }
    }
    
    // flush out our buffer
    CopyEdgeList(thread_data->options, edges, &edge_buffer);
    
    free(edge_buffer.from);
    free(edge_buffer.to);
    
    if (edges->weights_enabled)
        free(edge_buffer.weight);
    
    
     pthread_exit(NULL);
}



// no longer relies on the surface being a unit square
double *CreateQ(const GeometryStruct *g, const Options *options)
{
    double distance;
    double *Q;
    uint32_t i;
    uint32_t j;
    
    /* calloc initialises with zeros */
    double *t = calloc(options->M, sizeof(double));
    
    if (t == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    
    for (i = 1; i < options->M; i++) t[i] = (i - 1) * g->bucketSize;
    
    Q = malloc((size_t) sizeof(double) * g->Mx * g->My);
    
    if (Q == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    for(j = 0; j < g->My; j++)
    {
        for(i = 0; i < g->Mx; i++)
        {
            
            // distance = sqrt(t[i] * t[i] + t[j] * t[j]);
            distance = options->distance(t[i], t[j]);
            Q[i + g->Mx * j] =
            options->probGivenDistance(options->s1, options->s2,
                                       options->q, distance);
        }
    }
    
    free(t);
    return Q;
}



//------------------------------------------------------------------------------
//  ENTRYPOINT
//
//------------------------------------------------------------------------------

int GenSERN(NodeList* nodes, EdgeList* edges,
            Options* options, GeometryStruct* geometry)
{
    pthread_t *threads;
    pthread_attr_t attr;
    ThreadDataStruct *thread_data;
    int rc;
    uint32_t t;
    void *status;
    BucketStruct *buckets;
    double *Q;
    
    
    // see if we have alternatives for realloc and calloc
    // if not use the default ones provided by stdlib
    if (!options->realloc) options->realloc = realloc;
    if (!options->calloc) options->calloc = calloc;
    
    
    // start up the random number generator
    AllocRandom(options->seedval, options->ThreadCount);
    
    Q = CreateQ(geometry, options);
    
    //--------------------------------------------------------------------------
    //  NODE GENERATION
    //
    //--------------------------------------------------------------------------
    
    
    // some memory to return to the caller thus allocated with
    // provided memory allocation function
    nodes->x = options->calloc(options->N, sizeof(float));
    if (nodes->x == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    nodes->y = options->calloc(options->N, sizeof(float));
    if (nodes->y == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    
    // creates bucket structures with pre-initialised node counts
    // and their offsets into the node list.
    buckets = GenerateBuckets(options, geometry, nodes);
    
    
    // create an array of thread structures pre-initialised
    // with everything a worker thread needs to create nodes
    // or edges.
    
    threads = calloc(options->ThreadCount, sizeof(pthread_t));
    if (threads == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    thread_data = calloc(options->ThreadCount, sizeof(ThreadDataStruct));
    if (thread_data == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    for(t=0; t < options->ThreadCount; t++)
    {
        thread_data[t].thread_id = t;
        thread_data[t].thread_count = options->ThreadCount;
        thread_data[t].nodes = nodes;
        thread_data[t].edges = edges;
        thread_data[t].options = options;
        thread_data[t].buckets = buckets;
        thread_data[t].geometry = geometry;
        thread_data[t].Q = Q;
        thread_data[t].BufferSize = options->BufferSize;
        
    }
    
    
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    for(t=0; t < options->ThreadCount; t++)
    {
        
        rc = pthread_create(&threads[t], &attr,
                            (void *) &GenerateNodes,(void *) &thread_data[t]);
        if (rc)
        {
            
            options->errIdAndTxt(__FILE__, __LINE__,
                               "Error creating threads");

        }
    }
    
    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    
    for(t = 0; t < options->ThreadCount; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            
            options->errIdAndTxt(__FILE__, __LINE__,
                                 "Error joining threads");
        }
        
    }
   
    
    //--------------------------------------------------------------------------
    //  EDGE GENERATION
    //
    //--------------------------------------------------------------------------
    
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    
    openEdgeList(edges, options);
    
    for(t=0; t < options->ThreadCount; t++)
    {
        
        
        rc = pthread_create(&threads[t],
                            &attr, BusyWork,(void *)(&thread_data[t]));
        if (rc)
        {
            
            options->errIdAndTxt(__FILE__, __LINE__,
                                 "Error creating threads");
        }
    }
    
    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    
    for(t = 0; t < options->ThreadCount; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            options->errIdAndTxt(__FILE__, __LINE__,
                                 "Error joining threads");
        }
        
    }
    
    // if labeling of connected components is enabled then
    // we will return the labels of the components as they where
    // before we added links to connect them up as overwise use of this
    // option while ensuring the graph is connected would
    // be redundant
    
    if (options->components_enabled || options->connected)
        Components(options, nodes, edges);

    
    /* cleanup */
    free(buckets);
    free(Q);
    free(threads);
    free(thread_data);
    closeEdgeList();
    FreeRandom();
    
    return 0;
}

