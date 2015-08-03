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
#include <time.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // we nned the edgelist to be initialised empty
    EdgeList edges = {NULL, NULL, NULL, NULL, 0, 0, 0, 0};
    Options options = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int positive_arg[] = {1, 1, 1, 1, 1, 1, 0, 0, 1};
    
    int arg;
    //long int seedval;
    float *x,*y;
    double *n_edges;
    mwSize node_array_sz[2];
    char *ar[] = {"s", "q", "N", "M", "# of threads",
        "buffer size", "algorithm", "connected", "seed", NULL};
    char **i = ar;
    //uint32_t ThreadCount, BufferSize;
    
    
    /* Check for proper number of input and output arguments. */
    if ((nrhs < 6 ) || (nrhs > 9))
        mexErrMsgTxt("Incorrect number of input arguments: "
                     "[x,y,n_edges,edge_i,edge_j,edge_weights, node_components] "
                     "= waxman_gen3(s, q, N, M, #Threads, "
                     "BufferSize [, algorithm][, connected][, seed])");
    
    if (nlhs > 7)
        mexErrMsgTxt("Too many output arguments: "
                     "[x,y,n_edges,edge_i,edge_j,edge_weights, node_components] "
                     "= waxman_gen3(s, q, N, M, #Threads, "
                     "BufferSize [, algorithm][, connected][, seed]).");
    
    /* only allocate space for the "distances" if required */
    if (nlhs >= 6) options.weights_enabled = 1;
    
    /* only allocate space for node component identifiers if required */
    if(nlhs == 7) options.components_enabled = 1;
    
    // check arguments are positive scalars
    for (arg = 0; arg < nrhs; arg++)
    {
        
        
        if (mxGetM(prhs[arg]) != 1 ||
            mxGetN(prhs[arg]) != 1)
            
        {
            mexPrintf("\nArgument %i (%s) must be a scalar\n",
                      arg + 1, *i);
            mexErrMsgTxt("Error in arguments");
        }
        
        if    ((*(mxGetPr(prhs[arg]))) <= 0.0  && positive_arg[arg] == 1)
        {
            mexPrintf("\nArgument %i (%s) must be positive\n",
                      arg + 1, *i);
            mexErrMsgTxt("Error in arguments");
        }
        
        if    ((*(mxGetPr(prhs[arg]))) < 0.0  && positive_arg[arg] == 0)
        {
            mexPrintf("\nArgument %i (%s) must be non neagtive\n",
                      arg + 1, *i);
            mexErrMsgTxt("Error in arguments");
        }
        i++;
    }
    
    /* Get the input arguments */
    /* probability (i,j) connected = beta * exp(-s * distance(i,j) */
    options.s = mxGetScalar(prhs[0]);
    options.q = mxGetScalar(prhs[1]);
    
    /* number of nodes in the Waxman graph */
    options.N = (uint32_t) mxGetScalar(prhs[2]);
    
    /* dimension of array of buckets */
    options.M = (uint32_t) mxGetScalar(prhs[3]);
    
    /* number of threads to use */
    options.ThreadCount = (uint32_t) mxGetScalar(prhs[4]);
    
    /* number of edges in buffer for each thread */
    options.BufferSize = (uint32_t) mxGetScalar(prhs[5]);
    
    /* fast algorithm = 0 N^2 algorithm = 1 so default is fast */
    options.algorithm = (nrhs >= 7) ? (uint32_t) mxGetScalar(prhs[6]) : 0;
    
    /* default is to leave the graph disconnected */
    options.connected =  (nrhs >= 8) ? (uint32_t) mxGetScalar(prhs[7]) : 0;
    
    
    /* initialize random number generator */
    // TODO connect this to the random number generator in WaxmanGen
    options.seedval = (nrhs == 9) ? (uint32_t) mxGetScalar(prhs[8]) :
    (uint32_t) time(0);
    
    /* create output matrices */
    node_array_sz[0] = options.N;
    node_array_sz[1] = 1;
    plhs[0] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    x = (float *) mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    y = (float *) mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    n_edges = mxGetPr(plhs[2]);
    
    /* create the graph */
    edges.growth = 1.5 * options.N;
    
    WaxmanGen(x, y, &edges, &options);
    
    n_edges[0] = edges.count;
    
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
