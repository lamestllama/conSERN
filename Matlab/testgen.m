%
% basic test that the code functions
%
clear;
defaults;

s = 3;
q = 0.01; 
N = 100;
[x y e from to] = conSERN(DistanceFunction, s, q, N, metric, shape, A, connected, M, threads, algorithm, buffersize, seed);

figure(2)
clf
hold on
plot(x, y, '*');
plot([x(from), x(to)]', [y(from), y(to)]');
hold off

