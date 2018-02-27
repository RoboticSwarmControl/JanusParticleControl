% Authors: Li Huang and Aaron T. Becker
% Email: lhuang21@uh.edu
% All rights reserved
%=======================================
% Steering 4 Janus spheres to collide with greedy optimal control
function FourSpheresControl(n)
clc
if nargin < 1
% Set default number of Janus spheres 
   n = 4;
end

%% Select the control problem: we can either let 4 spheres converge to their mean position 
cf = 0.25*ones(4,1);    % cf = 0.25*ones(4,1): mu = mean position of x;
%% Or let 3 spheres chasing the 4th sphere till they converge.
% cf = [0;0;0;1]; %mu = x4
                       
format compact
rng(229)

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
mu = x_init*cf;

for ii = 1:n    
    x_goal(:,ii) = mu;
end
%Generate random thrust vectors  (unit magnitude, in R^3)
% thrustV= [u1 u2 u3 ... un;
%           v1 v2 v3 ... vn;
%           w1 w2 w3 ... wn];
thrustV =rand(3,n);

% Thrust vector normalization
for i = 1:n
   thrustV(:,i) = thrustV(:,i)./norm(thrustV(:,i));
end

% Calculate the collision spot
Rv = [thrustV(:,2)'-thrustV(:,1)',zeros(1,6);...
zeros(1,3),thrustV(:,2)'-thrustV(:,1)',zeros(1,3);...
zeros(1,6),thrustV(:,2)'-thrustV(:,1)';...
thrustV(:,3)'-thrustV(:,2)',zeros(1,6);...
zeros(1,3),thrustV(:,3)'-thrustV(:,2)',zeros(1,3);...
zeros(1,6),thrustV(:,3)'-thrustV(:,2)';...
thrustV(:,4)'-thrustV(:,3)',zeros(1,6);...
zeros(1,3),thrustV(:,4)'-thrustV(:,3)',zeros(1,3);...
zeros(1,6),thrustV(:,4)'-thrustV(:,3)'];
R_ = Rv\[x_init(:,1)-x_init(:,2);x_init(:,2)-x_init(:,3);x_init(:,3)-x_init(:,4)];
Rt = [R_(1:3)';R_(4:6)';R_(7:9)'];
collision_spot = x_init+Rt*thrustV;
fprintf('The unique collision spot is \n [x,y,z] = [%.03f, %.03f, %.03f]\n', collision_spot(1,1),collision_spot(2,1),collision_spot(3,1))


%<<<<<<<<<<<< Graph Init>>>>>>>>>>>>>
% Draw  the Janus spheres
figure(1);clf; 
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


%% Closed-loop control
numSteps = 5000;                    % Set the maximum number of time steps
Rt = eye(3);                        % rotation matrix for the next move
errors = ones(numSteps,1)*NaN;      % Sum of squared distance to goals

% Greedy optimal control mode switch
% mode = 2: Rt(k) = R(k-1)*Rx(a)*Ry(b)
% mode = 3: Rt(k) = R(k-1)*Rx(a)*Ry(b)*Rx(c)
mode = 3;                           

% Plot error
figure(2)
errHandler = plot(errors);
xlabel('Time')
ylabel('Sum Squared Distance')
set(gca,'FontSize',20);
x = x_init;                         % Current states (Janus sphere position)
Vt = 1e4;                           % Squared Euclidean distance error 
iop = 1;                            % Operation counter
t_step = 0.0.*ones(numSteps,1);     % time step

figure(1)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
set(gca,'FontSize',20);

%  While sum squared distance error is larger than a treshold
while((sum(Vt(:))>.01))

    mu = x*cf;
    for ii = 1:n    
        x_goal(:,ii) = mu;
    end
    Vt = (x(:)-x_goal(:))'*(x(:)-x_goal(:));
    errors(iop) = Vt;

    % Rotation matrix index for actuation
    % Consider a Lyapunov fcn candidate V(t)=(x_goal-x)'(x_goal-x) and its derivative
    % Vdot(t)= (x_goal-x)'(-x_dot) = (x_goal-x)'(-Rv)
    % V(t+t_opt) = V(t)+Vdot(t)*t_opt.
    % Find the optimal rotation matrix s.t. the magnitude of V_dot is
    % minimized

    
    switch mode     
        case 2            
            Vdot = @(theta)((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,4)-thrustV*cf )) );

            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2));

            % Compute an optimal actuation time (>=0) after each rotation 
            topt = max(0, -((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,4)-thrustV*cf )) )/n);
        case 3
            Vdot = @(theta)((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,4)-thrustV*cf )) );
            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3));

            % Compute an optimal actuation time (>=0) after each rotation 
            topt = max(0, -((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,4)-thrustV*cf )))/n);        
            
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

collision_spot = mean(x,2);
fprintf('Actual collision spot when error<0.01 where \n [x,y,z] = [%.03f, %.03f, %.03f]\n', collision_spot(1),collision_spot(2),collision_spot(3))


function R = Rx(a)
R = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];

function R = Ry(b)
R = [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];




