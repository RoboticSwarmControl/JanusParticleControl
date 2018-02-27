% Authors: Li Huang and Aaron T. Becker
% Email: lhuang21@uh.edu
% All rights reserved
%=======================================
% Steering 10 Janus spheres closer to each other with greedy optimal control
function TenSpheresControl(n)
clc
if nargin < 1
% Set default number of Janus spheres 
   n = 10;
end

cf = 0.1*ones(10,1);    
                       
format compact

%% Initialization
%<<<<<<<<<<<< Variables Init>>>>>>>>>>>>>
% Initialize Janus particle positions  
% x_init= [x1 x2 x3 ... xn;
%          y1 y2 y3 ... yn;
%          z1 z2 z3 ... zn];
x_init = randn(3,n)*5;
for ii = 1:n
    a = rand*2*pi;
    b = rand*2*pi;
    x_init(1,ii) = 5*cos(b)*sin(a);
    x_init(2,ii) = 5*sin(b)*sin(a);
    x_init(3,ii) = 5*cos(a);
end
    
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
    sphereHandler(i) =  surf(sx+x_init(1,i),sy+x_init(2,i),sz+x_init(3,i),'FaceColor',colors(i,:),'EdgeColor',.8*colors(i,:),'LineWidth',2);
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

x = x_init;                         % Current states (Janus sphere position)
Vt = 1e4;                           % Squared Euclidean distance error 
iop = 1;                            % Operation counter
t_step = 0.0.*ones(numSteps,1);     % time step

figure(1)
%  While sum squared distance error is larger than a treshold
while((sum(Vt(:))>1))
    % mu = x*cf = mean position of x;
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
    % minized

    
    switch mode     
        case 2            
            Vdot = @(theta)((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,4)-thrustV*cf ))...
                           +(x(:,5)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,5)-thrustV*cf ))...
                           +(x(:,6)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,6)-thrustV*cf ))...
                           +(x(:,7)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,7)-thrustV*cf ))...
                           +(x(:,8)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,8)-thrustV*cf ))...
                           +(x(:,9)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,9)-thrustV*cf ))...
                           +(x(:,10)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,10)-thrustV*cf )) );

            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2));

            % Compute an optimal actuation time (>=0) after each rotation 
            topt = max(0, -((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,4)-thrustV*cf ))...
                           +(x(:,5)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,5)-thrustV*cf ))...
                           +(x(:,6)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,6)-thrustV*cf ))...
                           +(x(:,7)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,7)-thrustV*cf ))...
                           +(x(:,8)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,8)-thrustV*cf ))...
                           +(x(:,9)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,9)-thrustV*cf ))...
                           +(x(:,10)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*(thrustV(:,10)-thrustV*cf ))  )/n);
        case 3
            Vdot = @(theta)((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,4)-thrustV*cf ))...
                           +(x(:,5)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,5)-thrustV*cf ))...
                           +(x(:,6)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,6)-thrustV*cf ))...
                           +(x(:,7)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,7)-thrustV*cf ))...
                           +(x(:,8)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,8)-thrustV*cf ))...
                           +(x(:,9)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,9)-thrustV*cf ))...
                           +(x(:,10)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,10)-thrustV*cf )) );
            theta = fminsearch(Vdot,[rand*2*pi-pi,rand*2*pi-pi,rand*2*pi-pi]);
            Rt = Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3));

            % Compute an optimal actuation time (>=0) after each rotation 
            topt = max(0, -((x(:,1)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,1)-thrustV*cf ))...
                           +(x(:,2)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,2)-thrustV*cf ))...
                           +(x(:,3)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,3)-thrustV*cf ))...
                           +(x(:,4)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,4)-thrustV*cf ))...
                           +(x(:,5)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,5)-thrustV*cf ))...
                           +(x(:,6)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,6)-thrustV*cf ))...
                           +(x(:,7)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,7)-thrustV*cf ))...
                           +(x(:,8)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,8)-thrustV*cf ))...
                           +(x(:,9)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,9)-thrustV*cf ))...
                           +(x(:,10)-mu)'*(Rt*Rx(theta(1))*Ry(theta(2))*Rx(theta(3))*(thrustV(:,10)-thrustV*cf )) )/n);

         
            
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




