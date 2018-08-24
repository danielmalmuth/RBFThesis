clc

% Functions
% q = @(x) exp(-(60.*(x-0.1)).^2); % initial distribution
% q_sol = @(u,x,t) exp(-(60.*((x-u*t)-0.1)).^2); % true solution
% q = @(x) cos(2*pi*x); % initial distribution
% q_sol = @(u,x,t) cos(2*pi*(x-u*t)); % true solution
q = @(x) exp(cos(2*pi*x)); % initial distribution
q_sol = @(u,x,t) q(x-u*t); % true solution

% Variables
u = -1; % constant in advection equation
target_time = 1;
numtimes = 6; % number of delta t to test
error = zeros(numtimes,1);
delta_t = zeros(numtimes,1); % vector of delta t

a = 0; % endpoints
b = 1;

delta_t = 0.125./2.^(1:numtimes);

h = (abs(u).*delta_t); % spatial difference based on delta t

for j = 1:numtimes
    n = (b-a)/h(j);
    x = linspace(a,b,n+1); x = x(1:n);
    q_vec = q(x)';

% % Centered Difference 2nd order
%     pos_center = (1/(2*h(j)))*ones((n),1);
%     neg_center = (-1/(2*h(j)))*ones((n),1);
%     D = spdiags([pos_center neg_center],[1 -1],n,n);
%     D(n,1) = (1/(2*h(j)));
%     D(1,n) = -(1/(2*h(j)));
%     q_approx = rk4_adv(delta_t(j),u,D,q_vec,target_time);
    
% % Centered Difference 4th order
%     pos2_center = (-1/(12*h(j)))*ones((n),1);
%     pos1_center = (8/(12*h(j)))*ones((n),1);
%     neg1_center = (-8/(12*h(j)))*ones((n),1);
%     neg2_center = (1/(12*h(j)))*ones((n),1);
%     D = spdiags([pos2_center pos1_center neg1_center neg2_center ...
%         pos2_center neg2_center],[2 1 -1 -2 -(n-2) (n-2)],n,n);
%     D(n,1) = (8/(12*h(j)));
%     D(1,n) = (-8/(12*h(j)));
%     q_approx = rk4_adv(delta_t(j),u,D,q_vec,target_time);

% Upwinding Scheme
    q_approx = fe_adv_up(delta_t(j),u,q_vec,target_time,h(j));
    
    
    q_true = q_sol(u,x,target_time)';
%     error(j) = max(abs(q_true - q_approx))/max(abs(q_true));
    error(j) = max(abs(q_vec - q_approx))/max(abs(q_vec));
end

loglog(delta_t,error,'x-')