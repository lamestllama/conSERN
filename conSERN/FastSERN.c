//
//  FastSERN.c
//  MexFastSERN
//
//  Created by Eric Parsonage on 3/30/15.
//
//

#include "FastSERN.h"
#include "binomial.h"
#include "uniform.h"
#include "geometric.h"
#include "edgelist.h"
#include "connected.h"
#include <pthread.h>



#ifdef MATLAB
#include "mex.h"
#define caller_realloc mxRealloc
#define caller_calloc mxCalloc
#else
#define caller_realloc realloc
#define caller_calloc  calloc
#endif





typedef struct
{
    float    *x;
    float    *y;
    uint64_t count;
    uint64_t start;
    int32_t i;
    int32_t j;
} BucketStruct;

typedef struct
{
    uint32_t     thread_id;
    uint32_t     thread_count;
    uint32_t      N;
    int          M;
    double       s;
    double       q;
    float        *x;
    float        *y;
    EdgeList     *edges;
    Options      *options;
    BucketStruct *buckets;
    double       *Q;
    uint32_t     BufferSize;
    
} ThreadDataStruct;



void SetBucketSizes (const uint32_t M, const uint32_t N, BucketStruct * buckets)
{
    uint32_t k;
    double p;
    double sum_p = 0.0;
    uint32_t sum_n = 0;
    
    
    /* the probability a randomly generated node
     will end up in a particular bucket t */
    p = 1.0/(double)(M * M);
    
    /* multinomial.  we have N nodes to distribute amongst
     M x M  boxes each of the boxes has equal probability
     of being chosen. Implemented one box at a time using
     binomial and conditioning on the distribution of the
     nodes in the previous boxes */
    for (k = 0; k < M * M ; k++)
    {
        buckets[k].count = GetBinomialRand(p / (1.0 - sum_p), N - sum_n);
        sum_p += p;
        sum_n += buckets[k].count;
    }
    
}




void GenerateNodes(uint32_t M, uint32_t N, BucketStruct *buckets)
{
    
    double x_offset, y_offset;
    uint32_t bucket_x = 0, bucket_y = 0, node_number;
    uint32_t i, j;
    // TODO this is what we need to thread
    for(i = 0; i < M * M; i++)
    {
        for(j = 0; j < buckets[i].count; j++)
        {
            buckets[i].x[j] = ((double)buckets[i].i + GetUniform(0)) / (double)M;
            buckets[i].y[j] = ((double)buckets[i].j + GetUniform(0)) / (double)M;
        }
    }
}


BucketStruct *GenerateBuckets(uint32_t M, uint32_t N, float *x, float *y)
{
    uint32_t i;
    uint32_t offset;
    
    /* calloc used to zero all fields */
    BucketStruct * buckets = calloc(M * M, sizeof(BucketStruct));
    
    
    SetBucketSizes (M,  N, buckets);
    
    
    /* Set each bucket to point at the appropriate part of the node arrays */
    offset = 0;
    buckets[i].start = 1;
    buckets[0].x = x;
    buckets[0].y = y;
    buckets[0].i = 0;
    buckets[0].j = 0;
    
    for(i = 1; i < M * M; i++)
    {
        buckets[i].start = buckets[i -1 ].start + buckets[i - 1].count ;
        buckets[i].x = x + buckets[i].start - 1;
        buckets[i].y = y + buckets[i].start - 1;
        buckets[i].i = i / M;
        buckets[i].j = i % M;
    
    }
    
    return buckets;
}



/* relies on the surface being a unit square */
double *CreateQ(uint32_t M, double s)
{
    double *Q;
    uint64_t i, j;
    double distance;
    
    /* calloc initialises with zeros */
    double *t = calloc(M, sizeof(double));
    
    for (i = 1; i < M; i++) t[i] = i - 1;
    
    Q = malloc((size_t) sizeof(double) * M * M);
    
    for(i = 0; i < M; i++)
    {
        for(j = 0; j < M; j++)
        {
            distance = sqrt(t[i] * t[i] + t[j] * t[j]) / M ;
            Q[i * M + j] =  exp(-s*distance);
        }
    }
    
    free(t);
    return Q;
}



void * BusyWork(void *t)
{
    
    int64_t i, j, k, S;
    uint32_t bucket_a, bucket_b;
    double p, lambda, distance, x_diff, y_diff;//, *Q;
    BucketStruct bucket_A, bucket_B;//, *buckets;
    
    ThreadDataStruct *thread_data;
    //int64_t     N;
    int         M;
    double      s;
    double      q;
    float       *x;
    float       *y;
    EdgeList    *edges;
    Options     *options;
    BucketStruct *buckets;
    double *Q;
    EdgeList edge_buffer = {NULL, NULL, NULL, NULL, 0, 0, 0, 0};
    
    
    thread_data = (ThreadDataStruct *)t;
    
    //N = thread_data->N;
    M = thread_data->M;
    s = thread_data->s;
    q = thread_data->q;
    x = thread_data->x;
    y = thread_data->y;
    edges = thread_data->edges;
    options =thread_data->options;
    buckets = thread_data->buckets;
    Q = thread_data->Q;
    
    
    
    AllocateEdgeBuffer(&edge_buffer, thread_data->BufferSize,
                       edges->weights_enabled);
    
    if (options->algorithm == 0)
    {
    
        /* Generate intra bucket edges in bucket_a */
        lambda = -log2(1 - q);
        	

        
        for (bucket_a = thread_data->thread_id; bucket_a < M * M;
             bucket_a += thread_data->thread_count)
        {
	    
            k = -1;
            bucket_A = buckets[bucket_a];
            
            S = (bucket_A.count * (bucket_A.count - 1))/2;
	    // fprintf(stdout, " q = %.8f, lambda = %.6f, bucket_A.count = %ld, S = %ld \n", q, lambda, bucket_A.count, S);  
            
            // find the next link that might exist
            while ((k += geom_rand2(lambda, thread_data->thread_id)) < S)
	    {
                // can we do an integer square root here ?
                j = 1 + ((((int64_t)sqrtl(8 * k + 1)) - 1) / 2);
                i = k - j * (j - 1) / 2; 

                x_diff = bucket_A.x[i] - bucket_A.x[j];
                y_diff = bucket_A.y[i] - bucket_A.y[j];
                distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                
                
                if (GetUniform(thread_data->thread_id) < exp(-s * distance))
                    AddEdgeToBuffer(edges, &edge_buffer, bucket_A.start + (uint32_t)i,
                                    bucket_A.start + (uint32_t) j, distance);
                
            }
            
        }
        
        
        // Generate inter bucket edges between bucket_a and bucket_b
        for (bucket_a = thread_data->thread_id; bucket_a < M * M;
             bucket_a += thread_data->thread_count)
        {
            bucket_A = buckets[bucket_a];
            
            for (bucket_b = bucket_a + 1; bucket_b < M * M; bucket_b++)
            {
                k = -1;
                bucket_B = buckets[bucket_b];
                
                p =  Q[abs(bucket_A.i - bucket_B.i) * M +
                       abs(bucket_A.j - bucket_B.j)];
                lambda = -log2(1 - p * q);
                
                S = bucket_A.count * bucket_B.count;
                
                // find the next link that might exist
                while ((k += geom_rand2(lambda, thread_data->thread_id)) < S)
                {
                    
                    i = k % bucket_A.count;
                    j = k / bucket_A.count;
                    
                    x_diff = bucket_A.x[i] - bucket_B.x[j];
                    y_diff = bucket_A.y[i] - bucket_B.y[j];
                    distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                    
                    if (GetUniform(thread_data->thread_id) * p < exp(-s * distance))
                        AddEdgeToBuffer(edges, &edge_buffer, bucket_A.start + (uint32_t) i,
                                        bucket_B.start + (uint32_t)j, distance);
                }
            }
        }
        
    }
    else
    {
        // for now any other option means use the N^2 algorithm
         /* Generate intra bucket edges in bucket_a */
        
        for (bucket_a = thread_data->thread_id; bucket_a < M * M;
             bucket_a += thread_data->thread_count)
        {

            bucket_A = buckets[bucket_a];
            
            for (i = 0; i < bucket_A.count; i++)
            {
                for(j = i + 1; j < bucket_A.count; j++)
                {
                    x_diff = bucket_A.x[i] - bucket_A.x[j];
                    y_diff = bucket_A.y[i] - bucket_A.y[j];
                    distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                    
                    if (GetUniform(thread_data->thread_id) < q * exp(-s * distance))
                        AddEdgeToBuffer(edges, &edge_buffer, bucket_A.start + (uint32_t)i,
                                        bucket_A.start + (uint32_t) j, distance);
                }
            }
        }
        
            
        // Generate inter bucket edges between bucket_a and bucket_b
        for (bucket_a = thread_data->thread_id; bucket_a < M * M;
             bucket_a += thread_data->thread_count)
        {
            bucket_A = buckets[bucket_a];
            
            for (bucket_b = bucket_a + 1; bucket_b < M * M; bucket_b++)
            {
                bucket_B = buckets[bucket_b];
                for ( i = 0; i < bucket_A.count; i++)
                {
                    for(j = 0; j < bucket_B.count; j++)
                    {
                        
                        x_diff = bucket_A.x[i] - bucket_B.x[j];
                        y_diff = bucket_A.y[i] - bucket_B.y[j];
                        distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                        
                        if (GetUniform(thread_data->thread_id) < q * exp(-s * distance))
                            AddEdgeToBuffer(edges, &edge_buffer, bucket_A.start + (uint32_t) i,
                                            bucket_B.start + (uint32_t)j, distance);
                    }
                }
            }
        }
    }
    
    // flush out our buffer
    CopyEdgeList(edges, &edge_buffer);
    
    free(edge_buffer.from);
    free(edge_buffer.to);
    
    if (edges->weights_enabled)
        free(edge_buffer.weight);
    
    
    pthread_exit((void*) 7);
}





int WaxmanGen(float* x, float* y, EdgeList* edges, Options* options)
{
    pthread_t *threads;
    pthread_attr_t attr;
    ThreadDataStruct *thread_data;
    int rc;
    uint32_t t;
    uint32_t component_count = 0;
    void *status;
    BucketStruct *buckets;
    double *Q;
    
    // fprintf(stdout, " s = %.8f, q = %.8f, N = %d, M = %d \n", options->s, options->q, options->N, options->M);

    SetSeed(options->seedval, options->ThreadCount); 
    
    Q = CreateQ(options->M, options->s);
    
    buckets = GenerateBuckets(options->M, options->N, x, y);
    
    // to do thread this function
    GenerateNodes(options->M, options->N, buckets);
    
    
    threads = calloc(options->ThreadCount, sizeof(pthread_t));
    thread_data = calloc(options->ThreadCount, sizeof(ThreadDataStruct));
    
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    
    openEdgeList(edges, options);
    
    for(t=0; t < options->ThreadCount; t++)
    {
        thread_data[t].thread_id = t;
        thread_data[t].thread_count = options->ThreadCount;
        
        thread_data[t].N = options->N;
        thread_data[t].M = options->M;
        thread_data[t].s = options->s;
        thread_data[t].q = options->q;
        thread_data[t].x = x;
        thread_data[t].y = y;
        thread_data[t].edges = edges;
        thread_data[t].options = options;
        thread_data[t].buckets = buckets;
        thread_data[t].Q = Q;
        thread_data[t].BufferSize = options->BufferSize;
        
        rc = pthread_create(&threads[t], &attr, BusyWork,(void *) &thread_data[t]);
        if (rc)
        {
            // TODO: some error handling here
            exit(-1);
        }
    }
    
    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    
    for(t = 0; t < options->ThreadCount; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            // TODO: some error handling here
            exit(-1);
        }
        
    }
    
    if(options->components_enabled)
    {
        // TODO memory allocation error handling
        // we want this memory zeroed hence the use of calloc
        edges->component = caller_calloc(options->N, sizeof(uint32_t));
        component_count = MarkConnectedComponents(options->N, edges);
    }
   
    // if labeling of connected components is enabled then
    // we will return the labels of the components as they where
    // before we added links to connect them up as overwise use of this
    // option while ensuring the graph is connected would
    // be redundant
    if(options->connected)
    {
        
        // if the labeling of connected components is not enabled
        // then label them so that we can work out what to connnect
        if (!options->components_enabled)
        {
            edges->component = caller_calloc(options->N, sizeof(uint32_t));
            component_count = MarkConnectedComponents(options->N, edges);
        }
        MakeConnected(options->N, edges, options->BufferSize,
                      x, y, component_count);
    }
    
    
    /* cleanup */
    free(buckets);
    free(Q);
    free(threads);
    free(thread_data);
    closeEdgeList();
    
    // free up the random number generator buffer;
    FreeUint();
    
    // fprintf(stdout, " N = %d, M = %d, edges = %d \n", options->N, options->M, edges->count);  

    return 0;
}

