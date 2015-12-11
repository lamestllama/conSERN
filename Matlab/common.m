%%% Various plotting options

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

% set default plot size
left = 0;
bottom = 0;
width = 16;
height = 10;
position = [ left bottom width height ];
set(0, 'DefaultFigurePaperPosition', position);
