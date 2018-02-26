clc; clear;

startPoint = [3.0;1.0];

tree_point = [0.0, 0.0;
              0.2, -0.05;
              -0.3, 0.08]';

finishPoint = tree_point(:, 1);

finishPoint_upp = [max(tree_point(1, :)); max(tree_point(2, :))];
finishPoint_low = [min(tree_point(1, :)); min(tree_point(2, :))];

max_force = 3.0;

problem.func.dynamics = @(t,x,u)( dynamics(x,u) );
problem.func.pathObj = @(t,x,u)( objective(x,u) );
problem.func.bndCst = @(t0,x0,tF,xF)( finalconstraint(x0,xF) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Set up bounds on state and control                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 5;
problem.bounds.finalTime.upp = 10;

problem.bounds.state.low = [-10; -5]; % bound along trajectory
problem.bounds.state.upp = [ 10; 5];

problem.bounds.initialState.low = startPoint - [0.1;0.1];
problem.bounds.initialState.upp = startPoint + [0.1;0.1];

problem.bounds.finalState.low = finishPoint_low;
problem.bounds.finalState.upp = finishPoint_upp;

problem.bounds.control.low = -max_force;
problem.bounds.control.upp = max_force;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Initialize trajectory with guess                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Car travels at a speed of one, and drives in a straight line from start
% to finish point.

del = finishPoint - startPoint;  % vector from start to finish

problem.guess.time = [0, norm(del)];   % time = distance/speed
problem.guess.state = [ startPoint, finishPoint];
problem.guess.control = [0, 0];  

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Options for Transcription                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'display','iter',...
    'MaxFunEval',1e5,...
    'tolFun',1e-6);

% problem.options.method = 'hermiteSimpson';
% problem.options.hermiteSimpson.nSegment = 25;

% problem.options.method = 'gpops';

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = 15;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);

t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
z = soln.interp.state(t);
x = z(1,:);
y = z(2,:);
u = soln.interp.control(t);
