%
% laplace_transform_test.m
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
% do a set of tests of Laplace transform calculations         
%
% standard cases from http://en.wikipedia.org/wiki/Laplace_transform
%   
%         f(t)                   F(s)
%       t^n / n!               1/s^(n+1)             Re{s}>0, n>-1
%       exp(-alpha*t)          1/(s+alpha),          Re{s}>-alpha       
%       sin(omega*t)           omega/(s^2+omega^2)   Re{s}>0
%       sinh(omega*t)          omega/(s^2-omega^2)   Re{s}>|omega|
%       cos(omega*t)           s/(s^2+omega^2)       Re{s}>0
%       cosh(omega*t)          s/(s^2-omega^2)       Re{s}>|omega|
%
% NB: there is a symbolic Laplace transform in the sym-manip toolkit we could use for
% verification
% 

clear;
ds = 0.1;
s_max = 10;
colors = [[1 0 0];
	  [1 0.5 0];
	  [1 1 0];
	  [0 1 0];
	  [0 1 1];
	  [0 0 1];
	  [0 0.5 0.5];
	  [0.5 0 0.5];
	 ];
range = [0, Inf];
range = [0, -log(eps)];  % t_max = - ln(epsilon)/min_s

%%% 
%%% Test the cases above
%%% 

%       t^n / n!               1/s^(n+1)             Re{s}>0, n>-1
figure(1)
hold off
plot(0,0)
hold on
s = 0.5:ds:s_max;
Ns = [1:5];
for i=1:length(Ns)
  n = Ns(i);
  f = @(t) t.^n / gamma(n+1);
  F = @(s) s.^(-n-1);
  Fcalc = laplace_trans(f, range, s);
  p1(i) = plot(s,F(s),'color', colors(i,:), 'linewidth', 2);
  plot(s,Fcalc, 'k--', 'linewidth', 2);
  legend_str(i,:) = sprintf('n = %2d', n);
end
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p1, legend_str);
xlabel('s');
ylabel('F(s)');
title('t^n/n!');
print('-depsc', 'Plots/laplace_trans_test_1.eps');

%       exp(-alpha*t)          1/(s+alpha),          Re{s}>-alpha       
figure(2)
hold off
plot(0,0)
hold on
alpha = 0.3;
s = 0:ds:s_max;
f = @(t) exp(-alpha*t);
F = @(s) 1./(s+alpha);
range1 = [0, Inf];
Fcalc1 = laplace_trans(f, range1, s);
range2 = range;
Fcalc2 = laplace_trans(f, range2, s);
p2(1) = plot(s,F(s),'b', 'linewidth', 2);
p2(2) = plot(s,Fcalc1, 'g:', 'linewidth', 4);
p2(3) = plot(s,Fcalc2, 'r--', 'linewidth', 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p2, 'F(s)', 'infinite interval', 'finite interval');
xlabel('s');
ylabel('F(s)');
title('L[exp(-\alpha t)]');
print('-depsc', 'Plots/laplace_trans_test_2.eps');

%       sin(omega*t)           omega/(s^2+omega^2)   Re{s}>0
figure(3)
hold off
plot(0,0)
hold on
omega = 4;
s = ds:ds:s_max;
f = @(t) sin(omega*t);
F = @(s) omega./(s.^2+omega.^2);
range1 = [0, Inf];
Fcalc1 = laplace_trans(f, range1, s);
range2 = range;
Fcalc2 = laplace_trans(f, range2, s);
p3(1) = plot(s,F(s),'b', 'linewidth', 2);
p3(2) = plot(s,Fcalc1, 'g:', 'linewidth', 4);
p3(3) = plot(s,Fcalc2, 'r--', 'linewidth', 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p3, 'F(s)', 'infinite interval', 'finite interval');
xlabel('s');
ylabel('F(s)');
title('L[sin]');
print('-depsc', 'Plots/laplace_trans_test_3.eps');


%       sinh(omega*t)          omega/(s^2-omega^2)   Re{s}>|omega|
figure(4)
hold off
plot(0,0)
hold on
omega = 4;
s = omega+ds:ds:s_max;
f = @(t) sinh(omega*t);
F = @(s) omega./(s.^2-omega.^2);
range1 = [0, Inf];
Fcalc1 = laplace_trans(f, range1, s);
range2 = range;
Fcalc2 = laplace_trans(f, range2, s);
p4(1) = plot(s,F(s),'b',  'linewidth', 2);
p4(2) = plot(s,Fcalc1, 'g:', 'linewidth', 4);
p4(3) = plot(s,Fcalc2, 'r--', 'linewidth', 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p4, 'F(s)', 'infinite interval', 'finite interval');
xlabel('s');
ylabel('F(s)');
title('L[sinh]');
print('-depsc', 'Plots/laplace_trans_test_4.eps');

%       cos(omega*t)           s/(s^2+omega^2)       Re{s}>0
figure(5)
hold off
plot(0,0)
hold on
omega = 4;
s = ds:ds:s_max;
f = @(t) cos(omega*t);
F = @(s) s./(s.^2+omega.^2);
range1 = [0, Inf];
Fcalc1 = laplace_trans(f, range1, s);
range2 = range;
Fcalc2 = laplace_trans(f, range2, s);
p5(1) = plot(s,F(s),'b',  'linewidth', 2);
p5(2) = plot(s,Fcalc1, 'g:', 'linewidth', 4);
p5(3) = plot(s,Fcalc2, 'r--', 'linewidth', 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p5, 'F(s)', 'infinite interval', 'finite interval');
xlabel('s');
ylabel('F(s)');
title('L[cos]');
print('-depsc', 'Plots/laplace_trans_test_5.eps');

%       cosh(omega*t)          s/(s^2-omega^2)       Re{s}>|omega|
figure(6)
hold off
plot(0,0)
hold on
omega = 4;
s = omega+ds:ds:s_max;
f = @(t) cosh(omega*t);
F = @(s) s./(s.^2-omega.^2);
range1 = [0, Inf];
Fcalc1 = laplace_trans(f, range1, s);
range2 = range;
Fcalc2 = laplace_trans(f, range2, s);
p6(1) = plot(s,F(s),'b',  'linewidth', 2);
p6(2) = plot(s,Fcalc1, 'g:', 'linewidth', 4);
p6(3) = plot(s,Fcalc2, 'r--', 'linewidth', 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p6, 'F(s)', 'infinite interval', 'finite interval');
xlabel('s');
ylabel('F(s)');
title('L[cosh]');
print('-depsc', 'Plots/laplace_trans_test_6.eps');


%%% 
%%% Test the first derivative of a couple of cases
%%% 

%       t^n / n!               1/s^(n+1)             Re{s}>0, n>-1
figure(10)
hold off
plot(0,0)
hold on
s = 0.5:ds:s_max;
Ns = [1:5];
for i=1:length(Ns)
  n = Ns(i);
  f = @(t) t.^n / gamma(n+1);
  F = @(s) (-1)*(n+1)*s.^(-n-1-1);
  Fcalc = laplace_trans(f, range, s, 1);
  p10(i) = plot(s,F(s),'color', colors(i,:), 'linewidth', 2);
  plot(s,Fcalc, 'k--', 'linewidth', 2);
  legend_str(i,:) = sprintf('n = %2d', n);
end
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p10, legend_str);
xlabel('s');
ylabel('F(s)');
print('-depsc', 'Plots/laplace_trans_test_11.eps');




%%% Test the second derivative of the cases above
%       t^n / n!               1/s^(n+1)             Re{s}>0, n>-1
figure(20)
hold off
plot(0,0)
hold on
s = 0.5:ds:s_max;
Ns = [1:5];
for i=1:length(Ns)
  n = Ns(i);
  f = @(t) t.^n / gamma(n+1);
  F = @(s) (n+1)*(n+2)*s.^(-n-1-2);
  Fcalc = laplace_trans(f, range, s, 2);
  p10(i) = plot(s,F(s),'color', colors(i,:), 'linewidth', 2);
  plot(s,Fcalc, 'k--', 'linewidth', 2);
  legend_str(i,:) = sprintf('n = %2d', n);
end
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'xlim', [0 10]);
legend(p10, legend_str);
xlabel('s');
ylabel('F(s)');
print('-depsc', 'Plots/laplace_trans_test_21.eps');







