/*
 *  conSERN.c
 *  conSERN
 *
 *  Created by Eric Parsonage on 8/2/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "conSERN.h"
#include "FastSERN.h"
#include "mex.h"
#include "edgeprobfuncs.h"
#include <time.h>

#define NUM_ELEMS(arr)                                                 \
(sizeof (struct {int not_an_array:((void*)&(arr) == &(arr)[0]);}) * 0 \
+ sizeof (arr) / sizeof (*(arr)))

#define DISTANCE_FUN_rhs    0
#define S_rhs               1
#define Q_rhs               2
#define N_rhs               3
#define M_rhs               4
#define THREADS_rhs         5
#define BUFFER_SIZE_rhs     6
#define ALGORITH_rhs        7
#define REGION_SHAPE_rhs    8
#define REGION_GEOMETRY_rhs 9
#define CONNECTED_rhs       10
#define SEED_rhs            11

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *data;
    mwSize node_array_sz[2];
    double *pEdgeCount;
    double *s;
    NodeList nodes;
    EdgeList edges;
    Options options;
    GeometryStruct *geometry;
    PolygonStruct *polygon;
    VectorStruct *vector;
    int i;
    int arg;

    int nonZero[] = {0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
    int dimM[]    = {1,  1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1};
    int dimN[]    = {1, -1, 1, 1, 1, 1, 1, 1, 1, -2, 1, 1};
    
    
    char *ar[] = { "distance function", "s", "q", "N", "M",
                    "number of threads", "buffer size",
                    "algorithm", "region shape", "region geometry",
                    "connected graph", "seed", NULL};
    char **c = ar;

    
    
    /* Check for proper number of input and output arguments. */
    if ((nrhs < 10 ) || (nrhs > 12))
        mexErrMsgTxt("Incorrect number of input arguments: "
                     "[x,y,n_edges,edge_i,edge_j,edge_weights, "
                     "node_components] = conSERN( Distancefunction, "
                     "s, q, N, M, #Threads, "
                     "BufferSize , algorithm, Distancefunction, shape, "
                     " geometry,[connected][, seed])");
    
    if (nlhs > 7)
        mexErrMsgTxt("Too many output arguments: "
                     "[x,y,n_edges,edge_i,edge_j,edge_weights, "
                     "node_components] = conSERN( Distancefunction, "
                     "s, q, N, M, #Threads, "
                     "BufferSize , algorithm, Distancefunction, shape, "
                     " geometry,[connected][, seed])");
  

    // check arguments have the correct dimension and
    // are positive scalars if necessary
    for (arg = 0; arg < nrhs; arg++)
    {
        
        // check dimensions match
        if (dimN[arg] > 0)
        {
            if (mxGetN(prhs[arg]) != dimN[arg])
            {
                mexPrintf("\nArgument %i (%s) must have dimension N = %i\n",
                          arg + 1, *c, dimN[arg]);
                mexErrMsgTxt("Error in arguments");
            }
        }
        else // dimN represents a minimum dimension
        {
            if (mxGetN(prhs[arg]) < -dimN[arg])
            {
                mexPrintf("\nArgument %i (%s) must have dimension N >= %i\n",
                          arg + 1, *c, -dimN[arg]);
                mexErrMsgTxt("Error in arguments");
            }
        }
        
        
        // check dimensions match
        if (dimM[arg] > 0)
        {
            if (mxGetM(prhs[arg]) != dimM[arg])
            {
                mexPrintf("\nArgument %i (%s) must have dimension M = %i\n",
                          arg + 1, *c, dimM[arg]);
                mexErrMsgTxt("Error in arguments");
            }
        }
        else // dimM represents a minimum dimension
        {
            if (mxGetM(prhs[arg]) < -dimM[arg])
            {
                mexPrintf("\nArgument %i (%s) must have dimension M >= %i\n",
                          arg + 1, *c, -dimM[arg]);
                mexErrMsgTxt("Error in arguments");
            }
            
        }
        
        // all the scalars are positive and some of them must be non zero
        if ((dimN[arg] == 1) && (dimN[arg] == 1))
        {
            
            if (*(mxGetPr(prhs[arg])) <= 0.0  && nonZero[arg] == 1)
            {
                mexPrintf("\nArgument %i (%s) must be positive\n",
                          arg + 1, *c);
                mexErrMsgTxt("Error in arguments");
            }
            
            if (*(mxGetPr(prhs[arg])) < 0.0  && nonZero[arg] == 0)
            {
                mexPrintf("\nArgument %i (%s) must be non negative\n",
                          arg + 1, *c);
                mexErrMsgTxt("Error in arguments");
            }
        }
        c++;
    }
    
    // before we start clear all our structures
    memset(&nodes, 0, sizeof(NodeList));
    memset(&edges, 0, sizeof(EdgeList));
    memset(&options, 0, sizeof(Options));
    
    
    /* check what optional outputs are required */
    /* only allocate space for the "distances" if required */
    if (nlhs >= 6) options.weights_enabled = 1;
    
    /* only allocate space for node component identifiers if required */
    if(nlhs == 7) options.components_enabled = 1;
    
    /* Get the input arguments */
    options.distanceFunction =  (uint32_t) mxGetScalar(prhs[DISTANCE_FUN_rhs]);
    
    if (options.distanceFunction >= NUM_ELEMS(distfunc))
    {
        mexErrMsgTxt("unimplemented distance function selected");
    }
    
    /* probability (i,j) connected various functions */
    options.distFunction = distfunc[options.distanceFunction];
    
    data = mxDuplicateArray(prhs[S_rhs]);
    
   // we already checked this
    
    if (mxGetN(data) == 1)
        options.s1 = mxGetScalar(prhs[S_rhs]);
    else // we have two S parameters
    {
        if (mxGetN(data) == 2)
        {
            s = (double *) mxGetData(data);
            options.s1 = s[1];
            options.s2 = s[2];
        }
        else
        {
            mexErrMsgTxt("Only two values for s allowed");
        }
        
    }
   
    options.q = mxGetScalar(prhs[Q_rhs]);
    
    /* number of nodes in the  graph */
    options.N = (uint32_t) mxGetScalar(prhs[N_rhs]);
    
    /* dimension of array of buckets */
    options.M = (uint32_t) mxGetScalar(prhs[M_rhs]);
    
    /* number of threads to use */
    options.ThreadCount = (uint32_t) mxGetScalar(prhs[THREADS_rhs]);
    
    /* number of edges in buffer for each thread */
    options.BufferSize = (uint32_t) mxGetScalar(prhs[BUFFER_SIZE_rhs]);
    
    /* fast algorithm = 0 N^2 algorithm = 1 so default is fast */
    options.algorithm =  (uint32_t) mxGetScalar(prhs[ALGORITH_rhs]);

   

    
    /* default is to leave the graph disconnected */
    options.connected =
        (nrhs >= 10) ? (uint32_t) mxGetScalar(prhs[CONNECTED_rhs]) : 0;
    
    /* initialize random number generator */
    options.seedval = (nrhs == 12) ? (uint32_t) mxGetScalar(prhs[SEED_rhs]) :
    (uint32_t) time(0);
    

    
    /* handle the shape and the geometry*/
    data = mxDuplicateArray(prhs[REGION_GEOMETRY_rhs]);
    vector  = (VectorStruct *) mxGetData(data);
    assert(mxGetM(data) == 2); // we already checked this
    
    polygon = polygonNew();
    
    for (i = 0;i <  mxGetN(data); i++) polygonAppend(polygon, vector + i);
    
    geometry = geometryGenerate( &options,
                                (GeometryType)
                                mxGetScalar(prhs[REGION_SHAPE_rhs]),
                                polygon);

    
    /* create output matrices */
    node_array_sz[0] = options.N;
    node_array_sz[1] = 1;
    plhs[0] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    nodes.x = (float *) mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    nodes.y = (float *) mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    pEdgeCount = mxGetPr(plhs[2]);
    
    /* create the graph */
    // TODO function to set initial size close to that needed
    edges.growth = 4 * options.N;
    
    GenSERN(&nodes, &edges, &options, geometry);
    polygonFree(polygon);
    pEdgeCount[0] = edges.count;
    
    plhs[3] = mxCreateNumericArray(0, 0, mxUINT32_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(0, 0, mxUINT32_CLASS, mxREAL);
    mxSetM(plhs[3], edges.count);
    mxSetN(plhs[3], 1);
    mxSetData(plhs[3], edges.from);
    mxSetM(plhs[4], edges.count);
    mxSetN(plhs[4], 1);
    mxSetData(plhs[4], edges.to);
    /* only allocate space for the "distances" if required */
    if (edges.weights_enabled)
    {
        plhs[5] = mxCreateNumericArray(0, 0, mxSINGLE_CLASS, mxREAL);
        mxSetM(plhs[5], edges.count);
        mxSetN(plhs[5], 1);
        mxSetData(plhs[5], edges.weight);
    }
    
    if(options.components_enabled)
    {
        plhs[6] = mxCreateNumericArray(0, 0, mxUINT32_CLASS, mxREAL);
        mxSetM(plhs[6], options.N);
        mxSetN(plhs[6], 1);
        mxSetData(plhs[6], edges.component);
    }
    
}
