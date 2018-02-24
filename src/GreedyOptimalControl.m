% This code simulates n magnetically steered spheres with catalytic Janus caps.
% The goal is to use a linear program to find the shortest path.
%  The paths consists of alternately revolving the magnetic field around
%  the current x and y axes, then allowing the Janus particles to self
%  propel a distance, then repeating.
function GreedyOptimalControl(n)
clc
if nargin < 1
% Set default number of Janus spheres 
   n = 3;
end

format compact
%  rng(74)

%  rng(15)
%  rng(18)
rng(78)

%% Initialization
%<<<<<<<<<<<< Variables Init>>>>>>>>>>>>>
% Initialize Janus particle positions  
% x_init= [x1 x2 x3 ... xn;
%          y1 y2 y3 ... yn;
%          z1 z2 z3 ... zn];
x_init = randn(3,n)*5;
%Initialize Janus particle goal positions
% x_goal= [x1 x2 x3 ... xn;
%          y1 y2 y3 ... yn;
%          z1 z2 z3 ... zn];
x_goal = randn(3,n)*2;
%Generate random thrust vectors  (unit magnitude, in R^3)
% thrustV= [u1 u2 u3 ... un;
%           v1 v2 v3 ... vn;
%           w1 w2 w3 ... wn];
thrustV =rand(3,n);

% Thrust vector normalization
for i = 1:n
   thrustV(:,i) = thrustV(:,i)./norm(thrustV(:,i));
end


%<<<<<<<<<<<< Graph Init>>>>>>>>>>>>>
% Draw  the Janus spheres
figure(1); clf;
sf = 10;
[sx,sy,sz] = sphere;
sx = sx/sf;sy=sy/sf;sz=sz/sf;

% Draw the sphere magnetic orientation (blue) and thrust orientation
% (green)
sphereHandler = ones(n,1);      % handles to spheres
trajHandler = ones(n,1);        % handles to paths
colors  = hsv(n);               % unique color for each sphere

for i = 1:n
    sphereHandler(i) =  surf(sx+x_init(1,i),sy+x_init(2,i),sz+x_init(3,i),'FaceColor',colors(i,:),'EdgeColor',.8*colors(i,:),'LineWidth',4);
    hold on
%   Draw the trajectory of each sphere
    trajHandler(i) = line(x_init(1,i),x_init(2,i),x_init(3,i),'Color',colors(i,:));
end
axis equal

% Draw the goal locations
th = 0:pi/12:2*pi;
flat = zeros(1,numel(th));
rad = 1/4;
for i = 1:n
    line( x_goal(1,i)+flat,x_goal(2,i)+rad*cos(th),x_goal(3,i)+rad*sin(th),'color','green');
    line( x_goal(1,i)+rad*cos(th),x_goal(2,i)+flat,x_goal(3,i)+rad*sin(th),'color','green');
    line( x_goal(1,i)+rad*cos(th),x_goal(2,i)+rad*sin(th),x_goal(3,i)+flat,'color','green');
end


%% Closed-loop control
% Greedy optimal control mode switch
% mode = 1: Rt(k) = R(k-1)*Rx(a) or Rt(k) = R(k-1)*Ry(a)
% mode = 2: Rt(k) = R(k-1)*Rx(a)*Ry(b)
% mode = 3: Rt(k) = R(k-1)*Rx(a)*Ry(b)*Rx(c)
mode = 3;       
numSteps = 5000;                    % Set the maximum number of time steps
errors = ones(numSteps,1)*NaN;      % Sum of squared distance to goals

% Plot error
figure(2)
errHandler = plot(errors);
xlabel('Time')
ylabel('Sum Squared Distance')
set(gca,'FontSize',20);

x = x_init;                         % Current states (Janus sphere position)
Vt = (x_goal-x).^2;                 % Position error (Lyapunov fcn)
iop = 1;                            % Operation counter
t_step = 0.0.*ones(numSteps,1);     % time step
Rt = eye(3);        % rotation matrix for the next move

figure(1)
set(gca,'FontSize',20);
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%  While sum squared distance error is larger than a treshold
while(sum(Vt(:))>0.1)

    Vt = (x_goal-x).^2;
    errors(iop) = sum(Vt(:));
 
    % Rotation matrix index for actuation
    % Consider a Lyapunov fcn candidate V(t)=(x_goal-x)'(x_goal-x) and its derivative
    % Vdot(t)= (x_goal-x)'(-x_dot) = (x_goal-x)'(-Rv)
    % V(t+t_opt) = V(t)+Vdot(t)*t_opt.
    % Find the optimal rotation matrix s.t. the magnitude of V_dot is
    % minized
    
    switch mode
        case 1
            if mod(iop,2)==1
                Vdot = @(theta) trace((x_goal-x)'*Rt*Rx(theta)*(-thrustV));
                theta = fminsearch(Vdot,rand*2*pi-pi);
                Rt = Rt*Rx(theta);
            else
                Vdot = @(theta) trace((x_goal-x)'*Rt*Ry(theta)*(-thrustV));
                theta = fminsearch(Vdot,rand*2*pi-pi);
                Rt = Rt*Ry(theta);
            end                
        case 2
            Vdot = @(theta) trace((x_goal-x)'*Rt*Rx(theta(1))*Ry(theta(2))*(-thrustV));
            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2));
        case 3
            Vdot = @(theta) trace((x_goal-x)'*Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(-thrustV));
            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3));
    end

    
    % Compute an optimal actuation time (>=0) after each rotation 
    topt = max(0.0, trace((x_goal-x)'*Rt*thrustV)/n);

    iop = iop + 1;          % Record time step    
    t_step(iop) = t_step(iop-1)+topt;
    
    % All Janus spheres move one step
    x = x + topt.*Rt*thrustV;
   
    for j = 1:n
        % Update the position of each sphere
        set( sphereHandler(j), 'XData', sx+x(1,j),'YData',sy+x(2,j),'ZData',sz+x(3,j));

        %update the path of each sphere
        xp= get(trajHandler(j),'XData');
        yp= get(trajHandler(j),'YData');
        zp= get(trajHandler(j),'ZData');
        set(trajHandler(j),'XData',[xp,x(1,j)],'YData',[yp,x(2,j)],'ZData',[zp,x(3,j)],'LineWidth',2);
    end
        
        set(errHandler,'Ydata',errors, 'Xdata',t_step, 'LineWidth', 2);
        drawnow
        
end


function R = Rx(a)
R = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];

function R = Ry(b)
R = [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];




