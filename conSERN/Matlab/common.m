
format long;

colors = [[0 0.9 0];
	  [0 0.9 0.9];
	  [0 0 0.9];
	  [1.0 0.8 0];
	  [0.8 0 0];
	  [0.9 0.4 0];
	  [0 0.5 0.5];
	  [0.5 0 0.5];
	 ];
symbols = ['p', 'o', 'd', 's', 'p', '^', 'v', '*'];
linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':'];
device = '-depsc';
suffix = 'eps';
seed = 1;
plotdir = 'Plots';
set(0,'DefaultTextFontsize', 16); % not working
set(0,'DefaultAxesFontsize', 16); % not working
set(0,'DefaultLineLinewidth', 1);
set(0,'DefaultAxesLineWidth', 1); 
version = '0.03';

% get code for doing line distances in a random region
path(path, '~/Reports/Networking/Topology/Waxman/Matlab/LinePicking/Matlab/');

% and Laplace transform code
path(path,'../../../Waxman/Matlab/Estimation/');

%  and MEX code for the Fast generation
path(path,'../');

% set default plot size
left = 0;
bottom = 0;
width = 16;
height = 10;
position = [ left bottom width height ];
set(0, 'DefaultFigurePaperPosition', position);

% problem settings
problem = 0; % square 
parameters = [1]; % side length of the square
[problem_name, description] = LinePickingProblemLookup(problem);
support = LinePickingSupport(problem, parameters);
g = @(t) LinePickingPDF(t, problem, parameters);

% algorithms
skipping = 0;
naive = 1;

% geometry shapes
region_rectangle = 0;
region_ellipse = 1;
region_polygon = 2;

% probability functions
pf_waxman=0;
pf_clipped_waxman=1;
pf_waxman_transition_threshold=2;
pf_threshold=3;
pf_constant=4;
pf_powerlaw=5;
pf_cauchy=6;
pf_exponential=7;
pf_maxentropy=8;

%distance functions
metric_euclidean=0;
metric_manhattan=1;
metric_discrete=2; 
metric_max=3;

%some sample regions
unit_square = [[0;0], [1;1]]; 
A1= [[9;11], [6;8], [4;5], [10;3]];
A2 = [[3;6],[4;8], [7;11],[10;7]];
A3 = [[50;150],[200;50],[350;150], [350;300],[250;300],[200;250],[150;350],[100;250],[100;200]];
A4 = [[1;1], [4;2]];

% default parameters
M=1;
threads=1;
algorithm=0;
buffersize=5000;
seed=1962;
connected=0;
shape = region_rectangle;
A = unit_square;
metric = metric_euclidean;
DistanceFunction = pf_waxman;

