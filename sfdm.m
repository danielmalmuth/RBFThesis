clc

% Functions
q = @(x) exp(-(20.*(x-0.1)).^2); % initial distribution
q_sol = @(u,x,t) exp(-(20.*((x-u*t)-0.1)).^2); % true solution

% Variables
u = 1; % constant in advection equation
target_time = 1;
numtimes = 20; % number of delta t to test
error = zeros(numtimes,1);
delta_t = zeros(numtimes,1); % vector of delta t
a = 0; % endpoints
b = 1;
for i = 1:numtimes
    delta_t(i) = (b-a)/(4*i*u);
end

h = (u.*delta_t)/2; % spatial difference based on delta t

for j = 1:numtimes
    n = (b-a)/h(j);
    q_vec = q(linspace(a,b,n))';
    pos_center = (1/(2*h(j)))*ones((n-1),1);
    neg_center = (-1/(2*h(j)))*ones((n-1),1);

    % Differentiation operator using centered differences
    D = diag(pos_center,1) + diag(neg_center,-1) + ...
        diag((-1/(2*h(j))),(n-1)) + diag((1/(2*h(j))),-(n-1));
    
    q_approx = be_adv(delta_t(j),u,D,q_vec,target_time);
    q_true = q_sol(u,linspace(a,b,n),target_time)';
    error(j) = abs(max(q_true - q_approx));
end

plot(delta_t,error)