%
% FastSERN_test1.m, (c) Matthew Roughan, 2015
%
% created: 	Wed Jul 15 2015 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% test for large s, how the algorithms scale with n
%
%
clear;
common;

s = 10;
n_step = 0.25;
max_n = 8;
% max_n = 6; 
ns = ceil( 10.^[2:n_step:max_n] );
max_naive = 10^5;
max_matlab = 10^4.25;
% max_naive = 10^4;
% max_matlab = 10^4;
k = 1;
Gs = laplace_trans(g, support, s);  
q = k ./ ((ns-1) * Gs)
 
runs = 100; % number of simulations to run 
Ms = [1, 5, 10, 15]; % bucket dimension for large s
 
% call programs once to make sure they are in memory
[N, E, d] = waxman_gen(ns(1), s, q(1), i, problem, parameters);
[x y e from to] = conSERN(DistanceFunction, s, q(1), ns(1), metric, shape, A, connected, Ms(1), threads, algorithm, buffersize, seed);

for i=1:runs
  fprintf('run %d of %d\n', i, runs);
  seed = i;
  for j=1:length(ns)
    n = ns(j);
          
    if n <= max_matlab
      tic();
      [N, E, d] = waxman_gen(n, s, q(j), seed, problem, parameters);
      times_matlab(i,j) = toc();
      edges_matlab(i,j) = length(d);
    end

    if n <= max_naive
      tic();
      [x y n_edges edge_i edge_j] = ...
	  conSERN(DistanceFunction, s, q(j), n, metric, shape, A, connected, 1, threads, naive, buffersize, seed);
      times_naive(i,j) = toc();
      edges_naive(i,j) = n_edges;
    end
    
    for k=1:length(Ms)
      M = Ms(k);
      
      tic();
      [x y n_edges edge_i edge_j] = ...
	  conSERN(DistanceFunction, s, q(j), n, metric, shape, A, connected, M, threads, algorithm, buffersize, seed);
      times(i,j,k) = toc();
      edges(i,j,k) = n_edges;
      
      legend_str(k,:) = sprintf('M = %2d', M);
    end
  end
end

T = squeeze(min(times));
E = squeeze(mean(edges));

T_naive = min(times_naive);
E_naive = mean(edges_naive);

T_matlab = min(times_matlab);
E_matlab = mean(edges_matlab);

figure(1)
hold off
semilogx(100,0)
hold on
plots = plot(ns, T);
plots(length(Ms)+1) = plot(ns(find(ns <= max_naive)), T_naive, 'r--');
plots(length(Ms)+2) = plot(ns(find(ns <= max_matlab)), T_matlab, 'b--');
legend_str(length(Ms)+1,:) = sprintf('naive ');
legend_str(length(Ms)+2,:) = sprintf('matlab');
legend(plots, legend_str, 'location', 'northwest');
set(gca, 'fontsize', 18);
set(gca, 'xtick',  10.^[2:ceil(max_n)]);
xlabel('n');
ylabel('time (seconds)');

filename = sprintf('Plots/FastSERN_test1_semilogx.eps', s)
print('-depsc', filename);


figure(2)
hold off
loglog(100,1)
hold on
plots = plot(ns, T);
plots(length(Ms)+1) = plot(ns(find(ns <= max_naive)), T_naive, 'r--');
plots(length(Ms)+2) = plot(ns(find(ns <= max_matlab)), T_matlab, 'b--');
legend_str(length(Ms)+1,:) = sprintf('naive ');
legend_str(length(Ms)+2,:) = sprintf('matlab');
legend(plots, legend_str, 'location', 'southeast');
set(gca, 'fontsize', 18);
set(gca, 'fontsize', 18);
set(gca, 'xlim',  10.^[2 ceil(max_n)]);
set(gca, 'xtick',  10.^[2:2:ceil(max_n)]);
set(gca, 'ylim',  10.^[-4 3]);
set(gca, 'ytick',  10.^[-4:2:3]);
xlabel('n');
ylabel('time (seconds)');
grid on

filename = sprintf('Plots/FastSERN_test1_loglog.eps', s)
print('-depsc', filename);

figure(3)
plot(ns, 2*E./repmat(ns'-1, 1, 2) );
