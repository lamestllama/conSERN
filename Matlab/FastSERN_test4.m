%
% FastSERN_test3.m, (c) Matthew Roughan, 2015
%
% created: 	Wed Jul 15 2015 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% push the size limit
%     only try this on 'maths1'
%
clear;
common;
defaults;
linepicking;

s = [3];
% n = 10.^9;  
% k = 10; 
n = 4*10.^9;  
k = 1; 
Gs = laplace_trans(g, support, s);  
q = k ./ ((n-1) * Gs);
small_q = 0.00001  ./ ((n-1) * Gs); 
     
runs = 1; % number of simulations to run
threads = [20];  
%threads = [1];  
buffersizes = [10.^[5]]; 
M = [20]; % bucket dimension for large s
  
% call program once to make sure it is in memory
[x y e from to] = conSERN(DistanceFunction, s, 0.1, 10, metric, shape, A, connected, M, threads(1), algorithm, buffersize(1), seed);

for i=1:runs
  fprintf('run %d of %d\n', i, runs);
  seed = i;
      
  for j=1:length(threads)  
    fprintf('      threads = %d\n', threads(j));
    
    % get just the time to create the nodes (edge prob ~ 0)
    tic();
    [x y n_edges edge_i edge_j] = ...
	conSERN(DistanceFunction, s, small_q, n, metric, shape, A, connected, 1, threads(j), algorithm, buffersizes(1), seed);
    times_baseline(i,j) = toc(); 
    
    for k=1:length(buffersizes)
      buffersize = buffersizes(k);
      fprintf('       buffersize = %d\n', buffersize);
       
      tic();
      [x y n_edges edge_i edge_j] = ...
	  conSERN(DistanceFunction, s, q, n, metric, shape, A, connected, M, threads(j), algorithm, buffersize, seed);
      times(i,j,k) = toc();
      edges(i,j,k) = n_edges;
       
      legend_str(k,:) = sprintf('B = %7d   ', buffersize);
    end
  end 
end
 

times


