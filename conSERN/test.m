
rectangle = 0;
ellipse = 1;
polygon = 2;

A1= [[9;11], [6;8], [4;5], [10;3]]
A2 = [[3;6],[4;8], [7;11],[10;7]]
A3 = [[50;150],[200;50],[350;150], [350;300],[250;300],[200;250],[150;350],[100;250],[100;200]]
A4 = [[1;1],[4;2]];


s=0.01;
q = 0.00075;
N = 45000;
M= 20;
Threads = 12;
Buffersize = 10000;
Algorithm = 0;
DistanceFunction = 0; % only q*exp(-s*dij) implemented
Regionshape = 0; %0 is a rectangle 1 is an ellipse and 2 is a polygon
RegionGeometry = A4; %bounding box coordinates for a Rectangle or ellipse
                     % 2 x n array for a poilygon
connected = 0; % cant do this at any scale yet
seed = 37;


A = A4
shape = rectangle
[x y n from to] = conSERN(s, q, N, M , Threads, Buffersize, Algorithm, DistanceFunction, shape, A, connected, seed);
figure(1)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A4
shape = ellipse
[x y n from to] = conSERN(s, q, N, M , Threads, Buffersize, Algorithm, DistanceFunction, shape, A, connected, seed);
figure(2)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A1
shape = polygon
[x y n from to] = conSERN(s, q, N, M , Threads, Buffersize, Algorithm, DistanceFunction, shape, A, connected, seed);
figure(3)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A2
shape = polygon
[x y n from to] = conSERN(s, q, N, M , Threads, Buffersize, Algorithm, DistanceFunction, shape, A, connected, seed);
figure(4)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off

A = A3
shape = polygon
[x y n from to] = conSERN(s, q, N, M , Threads, Buffersize, Algorithm, DistanceFunction, shape, A, connected, seed);
figure(5)
hold on
plot([A(1,:), A(1,1)] , [A(2,:), A(2,1)], '-')
plot(x, y, '*');
hold off


