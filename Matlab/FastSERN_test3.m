%
% FastSERN_test3.m, (c) Matthew Roughan, 2015
%
% created: 	Wed Jul 15 2015 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% test how the algorithm scales with threads
%
%
clear; 
common;

s = [3];
n = 10.^6;  
k = 10; 
Gs = laplace_trans(g, support, s);  
q = k ./ ((n-1) * Gs);
small_q = 0.00001  ./ ((n-1) * Gs); 
    
runs = 100; % number of simulations to run
threads = [1:1:12];  
buffersizes = [10.^[4:6]]; 
% buffersizes = [10.^[5]]; 
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
 
T = squeeze(min(times));
E = round(squeeze(mean(edges)));

T_baseline = squeeze(min(times_baseline));

figure(1)
hold off
plot(0,0)
hold on
plots = plot(threads, T);
plots(length(buffersizes)+1) = plot(threads, max(max(T))./threads, 'r--');
plots(length(buffersizes)+2) = plot(threads, T_baseline, 'g--');
legend_str(length(buffersizes)+1,1:5) = 'ideal';
legend_str(length(buffersizes)+2,1:13) = 'node creation';
legend(plots, legend_str, 'location', 'eastoutside');
set(gca, 'fontsize', 18);
% set(gca, 'xtick',  10.^[-1:1]);
set(gca, 'xlim',  [1  max(threads)]);
% set(gca, 'ytick',  10.^[-2:1:1]);
% set(gca, 'ylim',  [0 0.12]);
xlabel('threads');
ylabel('time (seconds)');

filename = sprintf('Plots/FastSERN_test3a.eps')
print('-depsc', filename);


figure(2)
hold off
plot(0,0)
hold on
plots2 = plot(threads, T ./ max(max(T)));
plots2(length(buffersizes)+1) = plot(threads, 1./threads, 'r--');
legend_str(length(buffersizes)+1,1:5) = 'ideal';
legend(plots2, legend_str, 'location', 'eastoutside');
set(gca, 'fontsize', 18);
% set(gca, 'xtick',  10.^[-1:1]);
set(gca, 'xlim',  [1  max(threads)]);
% set(gca, 'ytick',  10.^[-2:1:1]);
% set(gca, 'ylim',  [0 0.12]);
xlabel('threads');
ylabel('time (relative)');

filename = sprintf('Plots/FastSERN_test3b.eps')
print('-depsc', filename);


