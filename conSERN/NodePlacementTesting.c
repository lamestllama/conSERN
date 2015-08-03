/*
 *  NodePlacementTesting.c
 *  NodePlacementTesting
 *
 *  Created by Eric Parsonage on 7/24/15.
 *  Copyright 2015 __MyCompanyName__. All rights resulterved.
 *
 */

#include "NodePlacementTesting.h"
#include "binomial.h"
#include "mex.h"
#include <math.h>
#include <time.h>
#include <pthread.h>


static uint32_t *eric_m_z;
static uint32_t *eric_m_w;

#define INTP_TO_DOUBLEP 2.328306435454494e-10

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct
{
    double x;
    double y;
} VectorStruct2d; // maybe one day we will do 3d?


typedef struct
{
    uint32_t count;
    uint32_t allocated;
    VectorStruct2d *vertices;
} PolygonStruct;

typedef enum
{
    empty,
    partial,
    full
} BucketStatus;

typedef struct
{
    float    *x;
    float    *y;
    uint64_t count;
    uint64_t start;
    int32_t i;
    int32_t j;
    double  p;
    BucketStatus status;
    PolygonStruct *polygon;
} BucketStruct;


typedef enum
{
    rectangle,
    ellipse,
    polygon
} GeometryType;


typedef struct
{
    uint32_t N;
    uint32_t M;
    uint32_t Mx;
    uint32_t My;
    
    GeometryType type;
    VectorStruct2d origin;
    double xSize;
    double ySize;
    double totalArea;
    double bucketSize;
    PolygonStruct *polygon;
    
} GeometryStruct;


typedef struct
{
    uint32_t     thread_id;
    uint32_t     thread_count;
    BucketStruct *buckets;
    GeometryStruct *geometry;
} ThreadDataStruct;


typedef struct
{
    double a;
    double b;
    double c;
    double d;
} ellipseHelperStruct;


typedef struct
{
    VectorStruct2d v1;
    VectorStruct2d v2;
    VectorStruct2d v3;
    VectorStruct2d v4;
} vHelperStruct;

#define F(u, v) (asin(sqrt(1 - u * u) * sqrt(1 - v * v) -u * v) - u * sqrt(1 - u * u) - v * sqrt(1 - v * v) + 2 * u * v)
#define ellipse(v) ((((v.x) / (A)) * ((v.x) / (A))) + (((v.y) / (B)) * ((v.y) / (B))))

double ellipseRectangleIntersectArea(const GeometryStruct *g, const PolygonStruct *p)
{
    ellipseHelperStruct Ai[4];
    vHelperStruct V[4];
    VectorStruct2d pOffset[4];
    double S[4];
    double s = 0;
    int i;
    int j;
    double xOffset, yOffset;
    double A, B, Wdiv2, Hdiv2, xBar, yBar;
    
    
    A = g->xSize / 2.0;
    B = g->ySize / 2.0;
    Wdiv2 = g->bucketSize / 2.0;
    Hdiv2 = g->bucketSize / 2.0;
    
    // offsets to centre around the ellipse
    xOffset = g->origin.x + A;
    yOffset = g->origin.y + B;
    
    xBar = p->vertices[0].x - xOffset + Wdiv2;
    yBar = p->vertices[0].y - yOffset + Hdiv2;

    for(i = 0, j = 1; i < 4; i++, j++)
    {
        Ai[i].a= max(0,pow(-1, (j * j - j) / 2.0) * xBar - Wdiv2);
        Ai[i].b= max(0,pow(-1, (j * j + j - 2) / 2.0) * yBar - Hdiv2);
        Ai[i].c= max(0,pow(-1, (j * j - j)/ 2.0) * xBar + Wdiv2 - Ai[i].a);
        Ai[i].d= max(0,pow(-1, (j * j + j - 2) / 2.0) * yBar + Hdiv2 -Ai[i].b);
    }
    
    // set of 4 vertices for each quadrant
    for(i = 0; i < 4; i++)
    {
        V[i].v1.x = Ai[i].a;
        V[i].v1.y = Ai[i].b;
        V[i].v2.x = Ai[i].a;
        V[i].v2.y = Ai[i].b + Ai[i].d;
        V[i].v3.x = Ai[i].a + Ai[i].c;
        V[i].v3.y = Ai[i].b + Ai[i].d;
        V[i].v4.x = Ai[i].a + Ai[i].c;
        V[i].v4.y = Ai[i].b;
    }
    
    // look at each of the 6 cases in each of the 4 quadrants
    for(i = 0 ; i < 4; i++)
    {
        
        if (Ai[i].c < DBL_EPSILON || Ai[i].d <  DBL_EPSILON)
        {
            S[i] = 0;
            continue;
        }
        
        // V1 outside the ellipse meaning all vertices are outside the ellipse
        if (ellipse(V[i].v1) >= 1)
        {
            S[i] = 0;
            continue;
        }
        
        // only V1 inside the ellipse
        if ((ellipse(V[i].v1) < 1) &&
            (ellipse(V[i].v2) >= 1) &&
            (ellipse(V[i].v4) >= 1))
        {
            S[i] = (A * B / 2) * F((Ai[i].a / A), (Ai[i].b / B));
            continue;
        }
        
        // V4 inside V2 outside
        if ((ellipse(V[i].v4) < 1) &&
            (ellipse(V[i].v2) >= 1))
        {
            S[i] =  (A * B / 2) * (F((Ai[i].a / A), (Ai[i].b / B)) -
                                   F(((Ai[i].a + Ai[i].c )/ A), (Ai[i].b / B)));
            continue;
        }
        
        // V2 inside V4 outside
        if ((ellipse(V[i].v2) < 1) &&
            (ellipse(V[i].v4) >= 1))
        {
            S[i] =  (A * B / 2) * (F((Ai[i].a / A), (Ai[i].b / B)) -
                                   F((Ai[i].a / A), ((Ai[i].b + Ai[i].d) / B)));
            continue;
        }
        
        // only V2 and V4 inside V3 outside
        if ((ellipse(V[i].v2) < 1) &&
            (ellipse(V[i].v4) < 1) &&
            (ellipse(V[i].v3) >= 1))
        {
            S[i] =  (A * B / 2) * (F((Ai[i].a / A), (Ai[i].b / B)) -
                                   F(((Ai[i].a + Ai[i].c) / A), (Ai[i].b / B)) -
                                   F((Ai[i].a / A), ((Ai[i].b + Ai[i].d) / B)));
            continue;
        }
        
        // V3 inside V4 outside
        if (ellipse(V[i].v2) < 1)
        {
            S[i] = Ai[i].c * Ai[i].d;
            continue;
        }
        
        
    }
    
    for (int i = 0; i < 4; i++)
        s += S[i];
    
    return s;
}



static inline void SetSeed(uint32_t u, uint32_t ThreadCount)
{
    
    int i;
    
    srand(u);
    eric_m_z = malloc(sizeof(uint32_t) * ThreadCount);
    eric_m_w = malloc(sizeof(uint32_t) * ThreadCount);
    
    for(i = 0; i < ThreadCount; i++)
    {
        eric_m_z[i] = rand();
        eric_m_w[i] = rand();
    }
    
}

/* inline to avoid overhead of function call */
static inline uint32_t GetUint(uint32_t thread)
{
    uint32_t t;
    eric_m_z[thread] = 36969 * (eric_m_z[thread] & 65535) + (eric_m_z[thread] >> 16);
    eric_m_w[thread] = 18000 * (eric_m_w[thread] & 65535) + (eric_m_w[thread] >> 16);
    t = (eric_m_z[thread] << 16) + eric_m_w[thread];
    return t;
}


/* inline to avoid overhead of function call */
static inline double GetUniform(uint32_t thread)
{
    /* // 0 <= u < 2^32 */
    uint32_t u = GetUint(thread);
    /* // The magic number below is 1/(2^32 + 2). */
    /* // The resultult is strictly between 0 and 1. */
    return (u + 1.0) * INTP_TO_DOUBLEP;
}


// cross product of vectors
static inline double vecCross(const VectorStruct2d *a, const VectorStruct2d *b)
{
    return a->x * b->y - a->y * b->x;
}

// subtract vector b from vector a
static inline VectorStruct2d *vecSub(const VectorStruct2d *a,
                                     const VectorStruct2d *b,
                                    VectorStruct2d *result)
{
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    return result;
}


// dot product of vectors
static inline double vecDot(const VectorStruct2d *a, const VectorStruct2d *b)
{
    return a->x * b->x + a->y * b->y;
}

// Calculates the area of a parallelogram with two sides that have length
// and direction of a->b and b->c the area can be positive or negative
// depending on the anfle between the edges and this is what helps us
// determine what side c is on of the edge a->b.
// If c lies on the left side of directed edge a->b then return 1
// if c lies on the right side of directed edge a->b return -1
// otherwise c lies on the directed edge a->b return 0.

int edgeSide(const VectorStruct2d *a,
             const VectorStruct2d *b,
             const VectorStruct2d *c)
{
    double t;
    VectorStruct2d t1;
    VectorStruct2d t2;
    vecSub(b, a, &t1);
    vecSub(c, b, &t2);
    t = vecCross(&t1, &t2);
    return t > 0.0 ? 1 : t < 0.0 ? -1 : 0;
}


// we want to solve a0 + s * (a1 - a0) = b0 + t * (b1 - a0) for s and t
// cross product both sides by (a1 - a0) cancels out (a1 - a0) from lhs
// gives a0 X (a1 - a0)  = b0 X (a1 - a0) + t * (a1 - a0) X (b1 - a0)
// re-arranging for t gives:
// t = (a0 - b0) X (a1 - a0) / ((b1 - a0) X (a1 - a0))

int edgeIntersect(const VectorStruct2d *a0, const VectorStruct2d *a1,
                  const VectorStruct2d *b0, const VectorStruct2d *b1,
                  VectorStruct2d *result)
{
    VectorStruct2d aDirection;
    VectorStruct2d bDirection;
    VectorStruct2d baDirection;
    double t;
    double baCross;
    
    vecSub(a1, a0, &aDirection);
    vecSub(b1, b0, &bDirection);
    
    baCross = vecCross(&bDirection, &aDirection);
    
    // the edges are parallel so cannot intersect
    if (0 == baCross)return 0;
    
    vecSub(a0, b0, &baDirection);
    t = vecCross(&baDirection, &aDirection) / baCross;
    
    // intersection must lie between the points b0 and b1
    // ie on the line segment rather than anywhere on the line
    if (t >= 1.0 || t <= 0.0) return 0;
    
    result->x = b0->x + t * bDirection.x;
    result->y = b0->y + t * bDirection.y;
    return 1;
}


PolygonStruct *polygonNew()
{
    // calloc fills in allocated structure with zeros
    // TODO memory allocation error handling.
    return (PolygonStruct *) calloc(1, sizeof(PolygonStruct));
}

void polygonFree(PolygonStruct *p)
{
    // first free up the array of VectorStruct2D
    // that is used by the PolygonStruct this
    // POlygonstruct pointer points too
    free(p->vertices);
    // then free the memory that the PolygonStruct points too
    free(p);
}

void polygonAppend(PolygonStruct *p, const VectorStruct2d *v)
{
    // count should never get bigger than allocated
    assert(p->allocated >= p->count);
    
    if (p->count >= p->allocated)
    {
        // double the buffer
        p->allocated *= 2;
        // inital buffer is 4 enougth for a square
        // i.e. one of our buckets
        if (p->allocated == 0) p->allocated = 4;
        p->vertices = (VectorStruct2d *)
                    realloc(p->vertices, sizeof(VectorStruct2d) * p->allocated);
    }
    p->vertices[p->count++] = *v;
}


// ray casting to the right to see how many edges we cross
int32_t polygonInside(const PolygonStruct *polygon,
                      const float x, const float y)
{
    uint32_t  i;
    uint32_t  j;
    int32_t  c = 0;
    
    for (i = 0, j = polygon->count - 1; i < polygon->count; j = i++)
    {
        if (((polygon->vertices[i].y > y) != (polygon->vertices[j].y > y)) &&
            (x <
             (polygon->vertices[j].x - polygon->vertices[i].x) *
             (y-polygon->vertices[i].y) /
             (polygon->vertices[j].y-polygon->vertices[i].y) +
             polygon->vertices[i].x))
            c = !c;
    }
    return c;
}


// dependent on direction of polygon you can get negative areas
// this is useful for determining the direction of the polygon
double polygonAreaImp(const PolygonStruct *p)
{
    uint32_t i;
    uint32_t j;
    double area = 0;
    
    for (i = 0;i < p->count; i++)
    {
        j = (i + 1) % p->count;
        area += p->vertices[i].x * p->vertices[j].y;
        area -= p->vertices[i].y * p->vertices[j].x;
    }
    
    return area /= 2;
}

// always returns a positive area
double polygonArea(const PolygonStruct *p)
{
    return fabs(polygonAreaImp(p));
}


// this  only works if the polygon has at least 3 vertices
// that needs to be error checked on input regardless
int polygonDirection(const PolygonStruct *p)
{
    return (polygonAreaImp(p) > 0) ? 1 : -1;
}


void polygonClipEdge(PolygonStruct *subPoly,
                     VectorStruct2d *edgeOrigin, VectorStruct2d *edgeEnd,
                    int polygonDirection, PolygonStruct *result)
{
    uint32_t i;
    int32_t prevSide;
    int32_t currentSide;
    VectorStruct2d tmp;
    VectorStruct2d *prevVertex;
    VectorStruct2d *currentVertex;
    
    result->count = 0;
    
    
    // the last vertice defined in the polygon
    prevVertex  = subPoly->vertices + subPoly->count - 1;
    
    // which side of the edge under consideration is the current vertex ?
    prevSide = edgeSide(edgeOrigin, edgeEnd, prevVertex);
    
    // if it is on the inside of the edge being considered
    //  considered we keep it for now
    if (prevSide != -polygonDirection) polygonAppend(result, prevVertex);
    
    
    for (i = 0; i < subPoly->count; i++)
    {
        currentVertex = subPoly->vertices + i;
        currentSide = edgeSide(edgeOrigin, edgeEnd, currentVertex);
        
        // if edge between the previous vertex and the current one
        // straddle the edge under consideration add the intersection
        // point to the polygon being returned
        if (prevSide + currentSide == 0 && prevSide)
            if (edgeIntersect(edgeOrigin, edgeEnd, prevVertex, currentVertex, &tmp))
                polygonAppend(result, &tmp);
        
        // we already considered the last vertex in the polygon
        if (i == subPoly->count - 1) break;
        
        if (currentSide != -polygonDirection) polygonAppend(result, currentVertex);
        prevVertex = currentVertex;
        prevSide = currentSide;
    }
}


PolygonStruct *polygonClip(PolygonStruct *poly, PolygonStruct *clip)
{
    uint32_t i;
    PolygonStruct *in;
    PolygonStruct *out;
    PolygonStruct *tmp;
    int32_t direction;
    
    in = polygonNew();
    out = polygonNew();
    
    // get the direction of the clipping polygon
    direction = polygonDirection(clip);
    
    // clip the boundary polygon along the edge that implicitly closes
    // the clipping polygon i.e. between the first point given and the
    // last given they must not be coincident
    polygonClipEdge(poly, clip->vertices + clip->count - 1,
                    clip->vertices, direction, out);
    
    
    // clip the boundary polygon along each of
    // the other edges in the clipping polygon
    for (i = 0; i < clip->count - 1; i++)
    {
        // take the output of the last step
        // and feed it into the next step as input
        tmp = out; out = in; in = tmp;
        
        // stop if there are no vertices left inside the clip polgon
        if(in->count == 0)
        {
            out->count = 0;
            break;
        }
        polygonClipEdge(in, clip->vertices + i, clip->vertices + i + 1, direction, out);
    }
    
    polygonFree(in);
    return out;
}

uint32_t ellipseInside(const GeometryStruct *g, double x, double y)
{
    double x0, y0, xRadius, yRadius;
    
    xRadius = g->xSize / 2.0;
    yRadius = g->ySize / 2.0;
    x0 = g->origin.x + g->xSize / 2.0;
    y0 = g->origin.y + g->ySize / 2.0;
    
    return  ((((x - x0) * (x - x0)) / (xRadius * xRadius)) +
             (((y - y0) * (y - y0)) / (yRadius * yRadius)) < 1.0) ? 1 : 0;
}
             

// This routine takes a N = number of nodes and M = the desired
// number of buckets in dimension of the region with the largest
// distance from min to max and P a pointer to a polygon defining
// the region and initialises a GeometryStruct .
GeometryStruct *geometryGenerate(uint32_t N, uint32_t M,
                                 GeometryType type, PolygonStruct *p)
{
    uint32_t i;
    GeometryStruct *g;
    
    g = calloc(1, sizeof(GeometryStruct));
    assert(NULL != g);
    
    g->N = N;
    g->M = M;
    g->type = type;
    
    g->polygon = p;
    
    g->xSize = g->origin.x = g->polygon->vertices[0].x;
    g->ySize = g->origin.y = g->polygon->vertices[0].y;
    
    for (i = 0; i < g->polygon->count; i++)
    {
        g->origin.x = min(g->origin.x, g->polygon->vertices[i].x);
        g->origin.y = min(g->origin.y, g->polygon->vertices[i].y);
        g->xSize = max(g->xSize, g->polygon->vertices[i].x);
        g->ySize = max(g->ySize, g->polygon->vertices[i].y);
    }
    
    g->xSize -= g->origin.x;
    g->ySize -= g->origin.y;
    
    switch (g->type)
    {
        case rectangle:
            g->totalArea = g->xSize * g->ySize;
            break;
            
        case ellipse:
            g->totalArea = g->xSize * g->ySize * M_PI / 4.0;
            break;
            
        case polygon:
            g->totalArea = polygonArea(g->polygon);
            break;
            
        default:
            // this should never happen
            assert(1 != 0);
            break;
    }
    
    
    if (g->xSize > g->ySize)
    {
        g->Mx = g->M;
        g->bucketSize = g->xSize / g->Mx;
        g->My = (uint32_t)ceil(g->ySize / g->bucketSize);
      
    }
    else
    {
        g->My = g->M;
        g->bucketSize = g->ySize / g->My;
        g->Mx = (uint32_t)ceil(g->xSize / g->bucketSize);
       
    }
    
    return g;
}





void * GenerateNodes(void *t)
{
    ThreadDataStruct *thread_data;
    GeometryStruct *g;
    BucketStruct *buckets;
    uint32_t i, j;
    BucketStruct *bucket;
    
    thread_data = (ThreadDataStruct *)t;
    
    g = thread_data->geometry;
    buckets = thread_data->buckets;
    
    for (i = thread_data->thread_id; i < g->Mx * g->My;
         i += thread_data->thread_count)
    {
        bucket = buckets + i;
        switch (bucket->status)
        {
                // no points needed in this bucket
            case empty:
                break;
                
                // the whole bucket is to be covered no need to check
                // if points lie within a polygon boundary.
            case full:
                for (j = 0; j < bucket->count; j++)
                {
                    bucket->x[j] =  g->origin.x +
                    (((double) bucket->i + GetUniform(thread_data->thread_id)) * g->bucketSize);
                    bucket->y[j] =  g->origin.y +
                    (((double) bucket->j + GetUniform(thread_data->thread_id)) * g->bucketSize);
                }
                break;
                
                // the polygon crosses this bucket we need to accept
                //  a generated point if it lies within the polygon
                //  boundary and reject it otherwise
            case partial:
                
            switch (g->type)
            {
                case rectangle:
                    for (j = 0; j < bucket->count; j++)
                    {
                        bucket->x[j] =  g->origin.x +
                        (((double) bucket->i + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        bucket->y[j] =  g->origin.y +
                        (((double) bucket->j + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        
                        // reject the point if not inside
                        if ((bucket->x[j] > g->xSize + g->origin.x) ||
                            (bucket->y[j] > g->ySize + g->origin.y))
                            j--;
                    }
                    break;
                    
                case ellipse:
                    for (j = 0; j < bucket->count; j++)
                    {
                        bucket->x[j] =  g->origin.x +
                        (((double) bucket->i + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        bucket->y[j] =  g->origin.y +
                        (((double) bucket->j + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        
                        // reject the point if not inside
                        if (0 == ellipseInside(g, bucket->x[j], bucket->y[j]))
                            j--;
                    }
                    break;
                    
                case polygon:
                    
                    for (j = 0; j < bucket->count; j++)
                    {
                        bucket->x[j] =  g->origin.x +
                        (((double) bucket->i + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        bucket->y[j] =  g->origin.y +
                        (((double) bucket->j + GetUniform(thread_data->thread_id)) * g->bucketSize);
                        
                        // reject the point if not inside
                        if (!polygonInside(bucket->polygon,
                                           bucket->x[j], bucket->y[j]))
                            j--;
                    }
                    break;
            }
            break;
                
        }
    }
    pthread_exit((void*) 7);
}


// Knuth's double relative error comparison algorithm
int32_t close(const double x1, const double x2, const double epsilon)
{
    int exponent;
    double max, delta, difference;
    
    // exponent of largest absolute value
    max = (fabs (x1) > fabs (x2)) ? x1 : x2;
    frexp (max, &exponent);
    
    
    // Form a neighborhood of size  2 * delta
    delta = ldexp(epsilon, exponent);
    difference = x1 - x2;
    
    return (difference > delta)  ? 1 : (difference < -delta) ? 1 : 0;
    
}


void SetBucketSizes(const GeometryStruct* g, BucketStruct * buckets)
{
    uint32_t k;
 
    double sum_p = 0.0;
    uint32_t sum_n = 0;
    double norm = 0.0;
    
    
    // with enough buckets there is bound to be some
    // rounding error in the sum of probabilites
    // normalise these so they sum to 1
    for (k = 0; k < g->Mx * g->My ; k++)
        norm += buckets[k].p;
    for (k = 0; k < g->Mx * g->My ; k++)
        buckets[k].p /= norm;
    
 
    /* multinomial.  we have N nodes to distribute amongst
     Mx x My boxes. Each of the boxes have possibly different
     probabilities of being chosen based on the area of the bucket.
     Implemented one box at a time using binomial and conditioning 
     on the distribution of the nodes in the previous boxes */
    for (k = 0; k < g->Mx * g->My ; k++)
    {
        buckets[k].count = GetBinomialRand(buckets[k].p /
                                (1.0 - sum_p), g->N - sum_n);
        sum_p += buckets[k].p;
        sum_n += buckets[k].count;
    }
    
}



BucketStruct *GenerateBuckets(const GeometryStruct *g, float *x, float *y)
{
    uint32_t i;
    uint32_t j;
    uint32_t k;
    uint32_t c;
    uint32_t offset;
    double area;
    double bucketArea;
    double xExtent;
    double yExtent;
    BucketStruct *bucket;
    
    
    VectorStruct2d vertices[4];
    PolygonStruct  clippingPoly = {4,4, vertices};
    PolygonStruct *result;
    
    /* calloc used to zero all fields */
    BucketStruct * buckets = calloc(g->Mx * g->My, sizeof(BucketStruct));
    assert(NULL != buckets);
    
    bucketArea = g->bucketSize * g->bucketSize;
    
    for (j = 0 ; j < g->My; j++ )
    {
        for (i = 0; i < g->Mx; i++)
        {
            bucket = buckets + i + g->Mx * j;
            bucket->i = i;
            bucket->j = j;
       
            
            
            switch (g->type)
            {
                   
                case rectangle:
                   
                    xExtent = ((double) (i + 1)) * g->bucketSize;
                    yExtent = ((double) (j + 1)) * g->bucketSize;
                    
                    if ((xExtent > g->xSize) || (yExtent > g->ySize))
                    {
                        bucket->p =
                            (min(xExtent, g->xSize) -
                             ((double) i) * g->bucketSize) *
                            (min(yExtent, g->ySize) -
                             ((double) j) * g->bucketSize) /
                             g->totalArea;
                        
                        bucket->status = partial;
                        break;
                    }
                    
                    bucket->p = bucketArea / g->totalArea;
                    bucket->status = full;
                    break;
                    
                    
                case ellipse:

                    clippingPoly.vertices[0].x = i * g->bucketSize + g->origin.x;
                    clippingPoly.vertices[0].y = j * g->bucketSize + g->origin.y;
                    clippingPoly.vertices[1].x = clippingPoly.vertices[0].x;
                    clippingPoly.vertices[1].y = clippingPoly.vertices[0].y +
                    g->bucketSize;
                    clippingPoly.vertices[2].x = clippingPoly.vertices[0].x +
                    g->bucketSize;
                    clippingPoly.vertices[2].y = clippingPoly.vertices[1].y;
                    clippingPoly.vertices[3].x = clippingPoly.vertices[2].x;
                    clippingPoly.vertices[3].y = clippingPoly.vertices[0].y;
                    
                    for (c = k = 0; k < clippingPoly.count; k++)
                        c += ellipseInside(g, clippingPoly.vertices[k].x,
                                              clippingPoly.vertices[k].y);
                    
                    // there are two special cases when the ellipse and the
                    // bucket might intersect without the ellipse containing
                    // at least one of the buckets corners.
                    // if we have M odd either xSize or ySize are less than the
                    // size of the region covered by the buckets.
                    // we can test this by checking if the extreme top or
                    // extreme right of the ellipse lie inside the polygon
                    // TODO make this fire to check it is doing something
                    // I have not managed to yet
                    if(c == 0)
                    {
                        c += polygonInside(&clippingPoly,
                                          g->origin.x + g->xSize / 2.0,
                                          g->origin.y + g->ySize);
                        c += polygonInside(&clippingPoly,
                                           g->origin.x + g->xSize,
                                           g->origin.y + g->ySize / 2.0);
                        //mexPrintf("\ni = %i j= %i c=%i", i, j, c);
                    }
                    
                    switch(c)
                    {
                        case 0:
                            bucket->status = empty;
                            break;
                            
                        case 4:
                            bucket->status = full;
                            bucket->p = bucketArea / g->totalArea;
                            break;
                        
                        default:
                            bucket->status = partial;
                            bucket->p =
                              ellipseRectangleIntersectArea(g, &clippingPoly) /
                                g->totalArea;
                            break;
                    }
                    break;
                    
                case polygon:
                    
                    clippingPoly.vertices[0].x = i * g->bucketSize + g->origin.x;
                    clippingPoly.vertices[0].y = j * g->bucketSize + g->origin.y;
                    clippingPoly.vertices[1].x = clippingPoly.vertices[0].x;
                    clippingPoly.vertices[1].y = clippingPoly.vertices[0].y +
                    g->bucketSize;
                    clippingPoly.vertices[2].x = clippingPoly.vertices[0].x +
                    g->bucketSize;
                    clippingPoly.vertices[2].y = clippingPoly.vertices[1].y;
                    clippingPoly.vertices[3].x = clippingPoly.vertices[2].x;
                    clippingPoly.vertices[3].y = clippingPoly.vertices[0].y;
                    
                    result = polygonClip(g->polygon, &clippingPoly);
                    
                    switch (result-> count)
                    {
                        // this bucket is completely outside the polygon boundary
                        // there will be no need to place any nodes in it
                    case 0: polygonFree(result);
                        bucket->status = empty;
                        break;
                        
                        // if this bucket is completely within the polygon boundary
                        // there will be no need for the accept reject algorithm
                        // in the placement of points inside it
                    case 4: area = polygonArea(result);
                        bucket->p = area / g->totalArea;
                        if (0 == close(area, bucketArea, 0.0001))
                        {
                            polygonFree(result);
                            bucket->status = full;
                        }
                        else
                        {
                            bucket->polygon = result;
                            bucket->status = partial;
                        }
                        break;
                        
                        // this bucket is partially within the polygon boundary
                        // we need to use the accept reject algorithm to place
                        // points inside it
                    default: area = polygonArea(result);
                        bucket->p = area / g->totalArea;
                        bucket->polygon = result;
                        bucket->status = partial;
                    }
                    break;
                default:
                    break;
            
            }
            
        }
    }
    
    SetBucketSizes (g, buckets);

    /* Set each bucket to point at the appropriate part of the node arrays */
    offset = 0;
    buckets[0].start = 1;
    buckets[0].x = x;
    buckets[0].y = y;
    buckets[0].i = 0;
    buckets[0].j = 0;
    
    for(i = 1; i < g->Mx * g->My; i++)
    {
        buckets[i].start = buckets[i -1 ].start + buckets[i - 1].count ;
        buckets[i].x = x + buckets[i].start - 1;
        buckets[i].y = y + buckets[i].start - 1;
    }
    
    return buckets;
}

// no longer relies on the surface being a unit square
double *CreateQ(const GeometryStruct *g, double s)
{
    double distance;
    double *Q;
    uint64_t i;
    uint64_t j;
 
    /* calloc initialises with zeros */
    double *t = calloc(g->M, sizeof(double));
    assert(NULL != t);
    
    for (i = 1; i < g->M; i++) t[i] = (i - 1) * g->bucketSize;
    
    Q = malloc((size_t) sizeof(double) * g->Mx * g->My);
    assert(NULL != Q);
    for(j = 0; j < g->My; j++)
    {
    for(i = 0; i < g->Mx; i++)
    {
       
            distance = sqrt(t[i] * t[i] + t[j] * t[j]);
            Q[i + g->Mx * j] =  exp(-s*distance);
        }
    }
    
    free(t);
    return Q;
}



// there is going to be no checking of parameters as this is only for testing

void mexFunction(int nlhs,
                mxArray *plhs[], 
				 int nrhs, 
				 const mxArray *prhs[])
{

    
    pthread_t *threads;
    pthread_attr_t attr;
    ThreadDataStruct *thread_data;
    void *status;;
    
    uint32_t i;
    uint32_t t;
    uint32_t threadCount;
    int rc;
    mxArray *data;
    float *x, *y;
    mwSize node_array_sz[2];
    
    PolygonStruct *polygon;
    VectorStruct2d *vector;
    
    GeometryStruct *geometry;
    BucketStruct *buckets;
    
    
    threadCount = mxGetScalar(prhs[2]);
    SetSeed(time(0), threadCount);
    
    data = mxDuplicateArray(prhs[4]);
    vector  = (VectorStruct2d *) mxGetData(data);
    assert(mxGetM(data) == 2);
    
    polygon = polygonNew();
    
    for (i = 0;i <  mxGetN(data); i++)
        polygonAppend(polygon, vector + i);
    
    geometry = geometryGenerate( mxGetScalar(prhs[0]),
                                 mxGetScalar(prhs[1]),
                                 (GeometryType) mxGetScalar(prhs[3]),
                                 polygon);
    
    
    
    
    
    node_array_sz[0] = geometry->N;
    node_array_sz[1] = 1;
    plhs[0] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    x = (float *) mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(1, node_array_sz, mxSINGLE_CLASS, mxREAL);
    y = (float *) mxGetPr(plhs[1]);
    
    buckets = GenerateBuckets(geometry, x, y);
    
    threads = calloc(threadCount, sizeof(pthread_t));
    thread_data = calloc(threadCount, sizeof(ThreadDataStruct));
    
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    
    
    for(t=0; t < threadCount; t++)
    {
        thread_data[t].thread_id = t;
        thread_data[t].thread_count = threadCount;
        thread_data[t].buckets = buckets;
        thread_data[t].geometry = geometry;

        rc = pthread_create(&threads[t], &attr, (void *) &GenerateNodes,(void *) &thread_data[t]);
        if (rc)
        {
            // TODO: some error handling here
            exit(-1);
        }
    }
    
    /* Free attribute and wait for the other threads */
    
    for(t = 0; t < threadCount; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            // TODO: some error handling here
            exit(-1);
        }
        
    }
    pthread_attr_destroy(&attr);
    
    
    free(threads);
    free(thread_data);
    
//    
//    data = mxDuplicateArray(prhs[3]);
//    vector  = (VectorStruct2d *) mxGetData(data);
//    assert(mxGetM(data) == 2);
//    
//
//    for (i = 0;i <  mxGetN(data); i++)
//        polygonAppend(clippingPolygon, vector + i);
//        
//    
//    
//    result = polygonClip(polygon, clippingPolygon);
//    
//    
//    plhs[0] = mxCreateNumericArray(0, 0, mxDOUBLE_CLASS, mxREAL);
//    mxSetM(plhs[0], 2);
//    mxSetN(plhs[0], result->count);
//    
//    
//    
//    vector = (VectorStruct2d *)mxMalloc(sizeof(VectorStruct2d) * result->count);
//    
//    for (i = 0;i <  result->count; i++)
//        vector[i] = result->vertices[i];
//    
//
//              
//    
//    mxSetData(plhs[0],vector);
    

};