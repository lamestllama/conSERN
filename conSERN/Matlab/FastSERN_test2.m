%
% FastSERN_test2.m, (c) Matthew Roughan, 2015
%
% created: 	Wed Jul 15 2015 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% compare performance for a rank of values of s, for different M
%
%
clear;
common;

s_step = 0.25;
s = 10.^[-1:s_step:1];
n = 10^6;
k = 1; 
Gs = laplace_trans(g, support, s);  
q = k ./ ((n-1) * Gs)
 
runs = 100; % number of simulations to run
Ms = [1, 5, 10, 20, 40]; % bucket dimension for large s
 
% call program once to make sure it is in memory
[x y e from to] = conSERN(DistanceFunction, s(1), 0.1, 10, metric, shape, A, connected, Ms(1), threads, algorithm, buffersize, seed);

for i=1:runs
  fprintf('run %d of %d\n', i, runs);
  seed = i;
  for j=1:length(s)
    for k=1:length(Ms)
      M = Ms(k);
      
      tic();
      [x y n_edges edge_i edge_j] = ...
	  conSERN(DistanceFunction, s(j), q(j), n, metric, shape, A, connected, M, threads, algorithm, buffersize, seed);
      times(i,j,k) = toc();
      edges(i,j,k) = n_edges;
      
      legend_str(k,:) = sprintf('M = %2d', M);
    end
  end
end
 
T = squeeze(min(times));
E = round(squeeze(mean(edges)));

figure(1)
hold off
semilogx(100,0)
hold on
plots = plot(s, T);
legend(plots, legend_str, 'location', 'northwest');
set(gca, 'fontsize', 18);
% set(gca, 'xtick',  10.^[2:ceil(max_n)]);
xlabel('n');
ylabel('time (seconds)');

figure(2)
hold off
loglog(100,1)
hold on
plots = plot(s, T);
legend(plots, legend_str, 'location', 'northwest');
set(gca, 'fontsize', 18);
set(gca, 'xtick',  10.^[-1:1]);
set(gca, 'xlim',  [min(s) max(s)]);
set(gca, 'ytick',  10.^[-2:1:1]);
set(gca, 'ylim',  10.^[-1.5 1]);
xlabel('s');
ylabel('time (seconds)');
grid on

filename = sprintf('Plots/FastSERN_test2_loglog.eps', s)
print('-depsc', filename);
