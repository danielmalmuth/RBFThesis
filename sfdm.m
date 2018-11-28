% clc
close all

warning('off','MATLAB:nearlySingularMatrix')

% Functions
% q = @(x) exp(-(60.*(x-0.1)).^2); % initial distribution
% q = @(x) cos(2*pi*x); % initial distribution
q = @(x) exp(cos(2.*pi.*x)); % initial distribution
q_sol = @(u,x,t) q(x-u.*t); % true solution

% Variables
a = 0; % endpoints
b = 1;

u = @(t,x) -sin(2.*pi.*t).*(1+sin(2.*pi.*x)); % u in advection equation
target_time = 1;
numtimes = 4; % number of delta_t to test
max_u = 2;%u(fminbnd(@(t,x) -u(t,x),0,target_time)); % Max value of u function
cfl = 3; % CFL number

powers = [3 5 7 9];
degrees = [4 4 4 4];
numtests = length(powers);

error = zeros(numtimes,numtests);

delta_t = 1/16*(1./2.^(1:numtimes))'; % vector of delta_t values

h = delta_t*max_u/cfl; % hs based on delta_t values
for i = 1:numtests
    for j = 1:numtimes
        n = floor((b-a)/h(j));
        x = linspace(a,b,n+1); x = x(1:n)';
        q_vec = q(x);

    % % Centered Difference 2nd order with RK4
    %     pos_center = (1/(2*h(j)))*ones((n),1);
    %     neg_center = (-1/(2*h(j)))*ones((n),1);
    %     D = spdiags([pos_center neg_center],[1 -1],n,n);
    %     D(n,1) = (1/(2*h(j)));
    %     D(1,n) = -(1/(2*h(j)));
    %     q_approx = rk4_adv(delta_t(j),u,D,q_vec,target_time);

    % % Centered Difference 4th order with RK4
    %     pos2_center = (-1/(12*h(j)))*ones((n),1);
    %     pos1_center = (8/(12*h(j)))*ones((n),1);
    %     neg1_center = (-8/(12*h(j)))*ones((n),1);
    %     neg2_center = (1/(12*h(j)))*ones((n),1);
    %     D = spdiags([pos2_center pos1_center neg1_center neg2_center ...
    %         pos2_center neg2_center],[2 1 -1 -2 -(n-2) (n-2)],n,n);
    %     D(n,1) = (8/(12*h(j)));
    %     D(1,n) = (-8/(12*h(j)));
    %     q_approx = rk4_adv(delta_t(j),u,D,q_vec,target_time);

    % % Upwinding Scheme
    %     q_approx = fe_adv_up(delta_t(j),u,q_vec,target_time,h(j));

    % Semi-Legrangian Advection
        q_approx = rk4_sla(delta_t(j),u,x,(b-a),q_vec,target_time,...
            powers(i),degrees(i));



    %     q_true = q_sol(u,x,target_time)';
    %     error(j) = max(abs(q_true - q_approx))/max(abs(q_true));
        error(j,i) = max(abs(q_vec - q_approx))/max(abs(q_vec));
    end
end
%%
for i = 1:numtests
    loglog(h,error(:,i),'x-')
    hold on
end
xlabel('$h$','Interpreter','Latex')
ylabel('Relative $L_2$ Error','Interpreter','Latex')
title('Relative $L_2$ Error, Degree = 4','Interpreter','Latex')
legend('r^3 - O(5.0083)','r^5 - O(4.0957)','r^7 - O(4.2178)',...
    'r^9 - O(5.7320)')

orders = zeros(numtests,1);

for i = 1:numtests
    x = polyfit(log(h(2:end)),log(error(2:end,i)),1);
    orders(i) = x(1);
end