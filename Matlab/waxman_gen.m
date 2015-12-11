function [points, E, d] = waxman_gen(n, s, q, seed, problem, parameters)
%
% waxman_gen.m, (c) Matthew Roughan, 2015
% 
% Copyright 2012-15 Matthew Roughan <matthew.roughan@adelaide.edu.au>
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%   
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% WAXMAN_GEN: generate Waxman random graph where
%       n nodes are distributed uniformly in the unit square (for now) and 
%       where Waxman links are assigned with probability
%                p(d) = q * exp(-s*d)
%       for a potential link of distance d, where parameters
%                s > 0
%                q in (0,1]
%
%
% INPUTS:
%         n = number of nodes
%         (s, q) = Waxman parameters
%         problem says something about space
%         
% OUTPUTS:        
%         points = nx2 array (x_i, y_i) coordinates of points
%         E = ex2 array (i,j) of edges
%         d = ex1 array d_k of distance of k'th edge
%

if nargin < 5
  % default problem is the square
  problem = 0
end

if nargin < 6
  % default parameters
  parameters = parset(problem);
end

% random points from the region in question
[ points ] =  LinePickingSimPoints( n, problem, parameters, seed );

switch problem
 case {0, 1, 2, 3, 4, 5, 6, 11, 13, 15} % Euclidean distance metrics
  tmp = zeros(n, n);
  for j=1:size(points,1)
    x = points(j,:);
    [X] = meshgrid(x);
    tmp = tmp + (X - X').^2;
  end
  D = sqrt(tmp);
 case {9} % Manhattan distance
  D = zeros(n, n);
  for j=1:size(points,1)
    x = points(j,:);
    [X] = meshgrid(x);
    D = D + abs(X - X');
  end
 case {10,16,17} % Max distance
  D = zeros(n, n);
  for j=1:size(points,1)
    x = points(j,:);
    [X] = meshgrid(x);
    D = max(D, abs(X - X'));
  end
 otherwise
  % case {7,12} % Geodesic distance on sphere and hypersphere
  % case {8}    % Geodesic distance on prism
  % case {14}   % Geodesic distance on cylinder
  error('can''t currently generate this problem');
end

% Waxman probabilities
P = q * exp(-s*D);
P = triu(P, 1);
R = rand(size(P));
[k] = find(R <= P);
E = [ceil(k/n), mod(k,n)];

% length of selected edges
d = D(k);

