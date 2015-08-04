% geometry shapes
rectangle = 0;
ellipse = 1;
polygon = 2;

% probability functions
waxman=0
clipped_waxman=1
constant =2
threshold=3
powerlaw=4
cauchy=5
exponential=6
maxentropy=7
clipped_waxman2= 8

%distance functions
euclidean=0
manhattan=1
max=2
min=3


%some sample regions
A1= [[9;11], [6;8], [4;5], [10;3]]
A2 = [[3;6],[4;8], [7;11],[10;7]]
A3 = [[50;150],[200;50],[350;150], [350;300],[250;300],[200;250],[150;350],[100;250],[100;200]]
A4 = [[1;1],[4;2]];



DistanceFunction = waxman
s=0.01
q=0.00075
N = 45000
Metric = euclidean
shape = rectangle
geometry = A4
connected = 0

M=1
Threads=1
Algorithm=0
Buffersize=5000
seed=1962


[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(1)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A4
shape = ellipse
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(2)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A1
shape = polygon
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(3)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A1
shape = polygon
Metric = manhattan
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(4)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A2
shape = polygon
Metric = euclidean
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(5)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A3
shape = polygon
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(6)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A3
shape = polygon
DistanceFunction = powerlaw
s = [0.01 0.02]
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A, connected, M, Threads, Algorithm, Buffersize, seed);
figure(7)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

% TEST if defaults work too
A = A3
shape = polygon
DistanceFunction = powerlaw
s = [0.01 0.02]
[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A);
figure(8)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

% Show how to handle maxentropy
% the model defines the thining parameter as being used in the shape function
% I have not allowed the shape function to take q as a parameter yet
% so for the moment the maxentropy shape function takes two parameters
% s1 and q = s2  so you have to manually couple the q inside the shape function
% with the q used a thinning parameter thus
A = A3
shape = polygon
DistanceFunction = maxentropy
s = [0.01 0.02]
q = 0.02


[x y n from to] = conSERN(DistanceFunction, s, q, N, Metric, shape, A);
figure(9)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off


