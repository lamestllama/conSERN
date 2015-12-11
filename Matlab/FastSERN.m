function [x y e from to edge_weights node_components] = ...
    FastSERN(DistanceFunction, s, q, N, metric, shape, A, connected, M, threads, algorithm, buffersize, seed);
%
% FastSERN.m, (c) Matthew Roughan, 2015
%
% created: 	Thu Jul 16 2015 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% FASTSERN: shim file to make Matlab documentation work
%
% INPUTS:
% DistanceFunction =  probability deterrence functions
%                s = distance deterrence function parameter
%                q = thinning parameter
%                N = number of nodes
%           metric = distance metric
%            shape = region shape
%                A = parameters of region (e.g., length of sides)
%        connected = indicator of whether final graph should have arbitrary links added to
%                       make it connected (default = 0, leave as created)
%                M = grid size of array of buckets
%          threads = number of parallel threads to run
%        algorithm = 0 to use the fast algorith 1 to use the N^2 algorithm
%       buffersize = size of buffer to use before adding edges to main list
%             seed = seed in random number generator (optional)
%
% See defaults.m for lists of
%    1 - region shapes
%    2 - probability deterrence functions
%    3 - metrics
%         
% OUTPUTS:        
%           (x,y) = coordinates of node
%               e = number of edges
%       (from,to) = edge list
%    edge_weights = distances (optional)
% node_components = value corresponds to the component the node is part of (optional)
%
%
error('Don''t call this function -- its just a shim to create a help file.'); 
