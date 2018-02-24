% This code simulates n magnetically steered spheres with catalytic Janus caps.
% The goal is to use a linear program to find the shortest path.
%  The paths consists of alternately revolving the magnetic field around
%  the current x and y axes, then allowing the Janus particles to self
%  propel a distance, then repeating.
function ClosedloopLP(n)
clc
if nargin < 1
% Set default number of Janus spheres 
   n = 3;
end


%% Select linear programming mode
 LP_mode = 0;    % Low computation: do linear programming once till the gradient is too small
% LP_mode = 1;  % High compuation: do linear programming every time step


format compact
% rng(15)
% rng(18)
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

% Generate random rotation matrices
N = 201; 
f = ones(1,N);
% rotation matrices
revs = zeros(3,3,N);
revs(:,:,1) = eye(3);

for i = 2:2:N
    % for even number i, rotate about current magnetic x-axis 
    revs(:,:,i)=revs(:,:,i-1)*Rx(2*pi*rand);
    % for odd number i, rotate about current magnetic y-axis 
    revs(:,:,i+1)=revs(:,:,i)*Ry(2*pi*rand);
end


%% Linear Programming
% % Construct a linear programming problem 
% https://www.mathworks.com/help/optim/ug/linprog.html 
% [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options) 
% minimizes with the optimization options specified by options.

%<<<<<<<<<<<< LP Equality >>>>>>>>>>>>>
% [R R1v1 R2v1 R3v1 ... RNv1;
%  R R1v2 R2v2 R3v2 ... RNv2;
%  R R1v3 R2v3 R3v3 ... RNv3;
%  :
%  R R1vn R2vn R3vn ... RNvn;] t = [x_goal-x_init], where t is a 1 by N
% vector.
% => beq = [x_goal-x_init]
% => Aeq*t = beq.

Aeq = zeros(3*n,N);  % N rotation matrices candidates by 3 DOF of the robots

for i = 1:N
    newThrustOrients = revs(:,:,i)*thrustV;         % newThrustOrients = Ri*vi
    Aeq(:,i) = newThrustOrients(:);
end

beq = x_goal(:)-x_init(:);  % Location difference b/w the goals and the init positions

% [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options)
% Minimize the total time t from initial positions to goals given N rotation
% matrices, where t gives the actuation time of the thrust vector after each rotation,
% fval is the total time, and exitflag denotes if there exits a solution.
[t, fval, exitflag] = linprog(f,[],[],Aeq,beq,zeros(1,N),Inf*ones(1,N));
fprintf('Minimum time for the linear programming:%.02f\n',fval);
display(exitflag)

% Extract rotation matrices whose subsequent thrust vector has actuation time
% larger than a threhold
rotm = revs(:,:,t>1e-6);
% Extract the actuation time for the corresponding totation matrices rotm
tef = t(t>1e-6);


%<<<<<<<<<<<< Graph Init>>>>>>>>>>>>>
% Draw  the Janus spheres
figure(1); clf;
sf = 10;
[sx,sy,sz] = sphere;
sx = sx/sf;sy=sy/sf;sz=sz/sf;

% Draw the sphere magnetic orientation (blue) and thrust orientation
% (green)
sphereHandler = ones(n,1);      % handles to spheres
magneticHandler = ones(n,1);    % handles to magnet vectors
thrustHandler = ones(n,1);      % handles to thrust vectors
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
% Set the maximum number of time steps
numSteps = 5000;
errors = ones(numSteps,1)*NaN;

figure(2)
errHandler = plot(errors);
xlabel('Time')
ylabel('Sum Squared Distance')
set(gca,'FontSize',20);

x = x_init;                         % Current states (Janus sphere position)
Vt = (x_goal-x).^2;                 % Position error (Lyapunov fcn)
iop = 1;                            % Operation counter
t_step = 0.0.*ones(numSteps,1);     % time step
Rt = eye(3);                        % rotation matrix for the next move

figure(1)
set(gca,'FontSize',20);
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%  While sum squared distance error is larger than a treshold
while(sum(Vt(:))>0.1)
    topt = -1e-6;
    Vt = (x_goal-x).^2;
    errors(iop) = sum(Vt(:));
    
    % Rotation matrix index for actuation
    % Consider a Lyapunov fcn candidate V(t)=(x_goal-x)'(x_goal-x) and its derivative
    % Vdot(t)= (x_goal-x)'(-x_dot) = (x_goal-x)'(-Rv)
    % V(t+t_opt) = V(t)+Vdot(t)*t_opt.
    % Find the optimal rotation matrix s.t. the magnitude of V_dot*topt is
    % maximized
    % rot_axis = 0 -> perform rotation about current magnetic x-axis
    % rot_axis = 1 -> perform rotation about current magnetic y-axis

    for k = 1:length(tef)
        err_sum =trace((x_goal-x)'*(-rotm(:,:,k)*thrustV));            
        % Find the max topt by checking all rotation matrix candidates
        if topt < max(0,-err_sum/n) 
           Rt = rotm(:,:,k);
           topt = max(0.0,-err_sum/n); 
        end
    end          
  
    % Redo the linear programming
    if  topt<2e-2 || LP_mode    
        N = 51;
        f = ones(1,N);
        revs = zeros(3,3,N);
        revs(:,:,1) = eye(3);
        for i = 2:2:N
            % for even number i, rotate about current magnetic x-axis 
            revs(:,:,i)=revs(:,:,i-1)*Rx(2*pi*rand);
            % for odd number i, rotate about current magnetic y-axis 
            revs(:,:,i+1)=revs(:,:,i)*Ry(2*pi*rand);
        end
        Aeq = zeros(3*n,N);  % N rotation matrices candidates by 3 DOF of the robots
        for i = 1:N
            newThrustOrients = revs(:,:,i)*thrustV;             % newThrustOrients = Ri*vi
            Aeq(:,i) = newThrustOrients(:);
        end
        beq = x_goal(:)-x(:);         % Location difference b/w the goals and the init positions

        % [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options)
        % Minimize the total time t from initial positions to goals given N rotation
        % matrices, where t gives the actuation time of the thrust vector after each rotation,
        % fval is the total time, and exitflag denotes if there exits a solution.
        t = linprog(f,[],[],Aeq,beq,zeros(1,N),Inf*ones(1,N));
        clc
        % Extract rotation matrices whose subsequent thrust vector has actuation time
        % larger than a threhold
        rotm = revs(:,:,t>1e-6);
        % Extract the actuation time for the corresponding totation matrices rotm
        tef = t(t>1e-6);        
    end

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




