% Testing RBF as saddle point system becomes linearly dependent
% close all
clc

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

% Variables
n = 50; % number of nodes 
m = 10*n;
ngrid = 200; % number of evaluation locations for each coordinate
nloops = 250; % number of loops to test different alpha values
alpha = 0.01; % variable for projection
alphastart = 3e-0;
deg = 2;

solvers = {};
schurs = [];
truncs = [];

% Schur Versions

solvers = [solvers, "Backslash"];
schurs = [schurs, true];
truncs = [truncs, 0];
% 
% solvers = [solvers, "Pseudoinverse"];
% schurs = [schurs, true];
% truncs = [truncs, 0];
% 
solvers = [solvers, "SVD"];
schurs = [schurs, true];
truncs = [truncs, 0];
% 
% solvers = [solvers, "QR"];
% schurs = [schurs, true];
% truncs = [truncs, 0];
% 
% solvers = [solvers, "GECP"];
% schurs = [schurs, true];
% truncs = [truncs, 0];

% Non-Schur Versions

solvers = [solvers, "Backslash"];
schurs = [schurs, false];
truncs = [truncs, 0];
% 
% solvers = [solvers, "Pseudoinverse"];
% schurs = [schurs, false];
% truncs = [truncs, 1e-10];
% 
% solvers = [solvers, "SVD"];
% schurs = [schurs, false];
% truncs = [truncs, 0];
% 
% solvers = [solvers, "QR"];
% schurs = [schurs, false];
% truncs = [truncs, 1e-15];
% 
solvers = [solvers, "GECP"];
schurs = [schurs, false];
truncs = [truncs, 0];


alphavec = zeros(nloops+1,1);
errors = zeros(nloops+1,1);
% errors2 = zeros(nloops+1,1);
% errors3 = zeros(nloops+1,1);
% errors4 = zeros(nloops+1,1);
% errors5 = zeros(nloops+1,1);
% gamma = 10e-6; % parameter for pinv method
% ranks = zeros(nloops+1,1);

% Function to approximate
f = @(x,y) (x+y).*exp((1-x.^2-y.^2));
% f = @(x,y) x+y + exp(y-x);

% Generate alpha value
for i = 1:length(alphavec)
    alphavec(i) = alphastart*0.90^(i-1);
end

% Generate 2D pseudo-random point set
% [ux,uy] = meshgrid(linspace(0,1,ngrid),linspace(0,1,ngrid));
ux = linspace(0,1,ngrid);
uy = linspace(0,1,ngrid);

% p = haltonset(1,'Skip',1);
% % vec = net(p,n);
% rng(1234);
% vec = randn(n,1)/100;
% x = linspace(0,1,n);
% y = linspace(0,1,n);
% v = linspace(0,(n-1),n);
% xadd = (((-1).^(v+1))'.*vec)';
% yadd = (((-1).^(v))'.*vec)';
% x = x + xadd;
% y = y + yadd;
% 
% scatter(x,y,'.')

x0 = 0; % x0 and y0 are center coordinates
y0 = 0;  
r = 1;  % radius
angle = linspace(0,pi,n);   %(2*pi-(2*pi)/n),n);
rng(2345);
vec = randn(n,1)/10;
x=(r+vec').*cos(angle) + x0;
y=(r+vec').*sin(angle) + y0;
% scatter(x,y,'.')
% axis square

ux=r.*cos(linspace(0,pi,m)) + x0;  %(2*pi-(2*pi)/m),m)) + x0;
uy=r.*sin(linspace(0,pi,m)) + y0;   %(2*pi-(2*pi)/m),m)) + y0;

for test = 1:length(solvers)
    for i = 1:length(alphavec)
        % Scale the normal component of points by alpha
        % alpha = 0 projects point onto line y = x
        % alpha = 1 does not affect points
        
        alpha = alphavec(i);
        
%         [xp,yp] = diagproj(x,y,alpha);

        xp = (r+(alpha.*vec')).*cos(angle) + x0;
        yp = (r+(alpha.*vec')).*sin(angle) + y0;
        xp = xp';
        yp = yp';
        
       

        lam = rbffit_test(xp,yp,f(xp,yp),solvers(test),schurs(test),truncs(test),deg);
        g = rbfval2(lam, xp, yp, ux, uy,deg);

        errors(i) = max(max(abs(g-f(ux,uy))));
        % alpha = alpha - (alphastart/nloops);
    end
    loglog(alphavec,errors)
    if test == 1
        hold on
    end
end
xlabel('Alpha')
ylabel('Maximum error')
title('Max Error')

legendarray = {};
for i = 1:length(solvers)
    if schurs(i)
        legendarray = [legendarray, solvers(i) + " " + "(Schur)"];
    else
        legendarray = [legendarray, solvers(i) + " " + "(Non-Schur)"];
    end
end
legend(legendarray)
set(gcf, 'Position', get(0, 'Screensize'));
hold off

% figure
% semilogx(alphavec,ranks)
% xlabel('Alpha')
% ylabel('Rank')
% title('Rank over Tested Alpha Values')

% surf(ux,uy,f(ux,uy))
% title('Surface Plot')
% shading interp
% colormap summer
% hold on
% scatter3(xp,yp,f(xp,yp),'k','.')