function [F, errbnd] = laplace_trans(f, range, s, n, varargin)
%
% laplace_trans.m
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
% LAPLACE_TRANS: compute the nth derivative of the Laplace transform of a function:
%
%               F^{(n)}(s) = \int_0^{\infty}  f(t) (-t)^n e^{-st} dt
% 
%    see http://en.wikipedia.org/wiki/Laplace_transform
%
% INPUTS:
%         f     = a pointer to the function
%         range = the range over which to integrate
%         s     = a (row) vector of values for which to calculate the Laplace transform
%         n (default value = 0) includes a power of (-t)^n in the Laplace transform
%                        n = 0, standard Laplace transform F(s)
%                        n = 1, negative of first derivative of the LT, F'(s)
%                              in general
%                                 L[ (-t)^n f(t) ] =  F^{(n)}(s)
%         varagin = is a variable set of arguments that will passed to the 
%                   quadrature function to be used (we use quadv so that we can 
%                   compute the integrals for a range of values of s together).
%
% OUTPUTS:        
%         F     = a row vector giving the Laplace transform at points s
%         
%
%
if (nargin < 4)
  n = 0;
end
s = s(:)'; % turn it into a row vector

tol = 1.0e-12;
h = @(t) exp(-s*t) .* (-t).^n .* f(t);
[F] = integral(h, range(1), range(2), 'RelTol', tol, 'ArrayValued', true);
errbnd = Inf;

