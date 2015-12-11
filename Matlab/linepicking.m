% use this to bring in code from LinePicking, and 
% to set up basic test scenario used in tests


% get code for doing line distances in a random region
path(path, '~/Reports/Networking/Topology/Waxman/Matlab/LinePicking/Matlab/');

% problem settings
problem = 0; % square 
parameters = [1]; % side length of the square
[problem_name, description] = LinePickingProblemLookup(problem);
support = LinePickingSupport(problem, parameters);
g = @(t) LinePickingPDF(t, problem, parameters);
