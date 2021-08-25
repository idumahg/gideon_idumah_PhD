clear;clc
tic

% parameter values
n = 10; % number of circles
kappa = [1 0;0 1]; % diffusion tensor
C = 32; % number of unknown
lamda = 1;
tspan = 20; % time span

% boundary values for each of the metabolites in ECS and Astrocytic compartments 
Ub = [1.1289; 0.0165; 83.75; 1.247; 0.0103; 0.0100; 0.7536; 0.0098; 85.46; 1.2472; ...
    0.0103; 0.0100; 0.1202; 2.189; 0.0227; 0.0299; 0.0031];

Ubar = [];
for j = 1:length(Ub)
    [V,M,K,B,U] = matrices(n,kappa,Ub(j));
    Ubar = [Ubar;U];
end
N = size(V,2); % number of nodes

% boundary vector
BigU0 = sparse(C*N,1);

BigU0(4*N+1:10*N) = Ubar(1:6*N);
BigU0(21*N+1:end) = Ubar(6*N+1:end);

% other matrices
BigK = kron(eye(C), K); % coupled stiffness matrix
BigM = kron(eye(C), M); % coupled mass matrix
BigB = kron(eye(C), B); % coupled boundary matrix


eta1 = 0.04; eta2 = 0.26; eta3 = 0.4; eta4 = 0.3; % volume fraction of each compartment

Psi = spdiags([eta1*ones(4*N,1); eta2*ones(6*N,1); eta3*ones(11*N,1); eta4*ones(11*N,1)], 0, C*N, C*N); % coupled volume fraction

%xi = @(t) 20*gampdf(t,22,0.8); % time function
mu = 2; sigma = 1;
xi = @(t) 20*gampdf(t,10,1); % proxy for electrophysiology input into metabolic model 

[u_sol,t_sol] = diffusion(N,BigK,BigM,BigB,BigU0,lamda,tspan,Psi,xi);

timeElapsed = toc;

% [U1,U2,U3,U4] = unpack(u_sol', N); % struct of solution {'Glc'; 'oxy'; 'CO2'; 'Lac'; 'Glu'; 'Gln'; 'Pyr'; 'ATP'; 'ADP'; 'NADplus'; 'NADH'};

T = delaunay(V(1,:)',V(2,:)');

% for k = 1:size(u_sol, 1)
%     h = figure;
%     u = u_sol(k,:);
%     %trisurf(T,V(1,:),V(2,:),u,'EdgeColor','none')
%     %shading interp
%  
%     trisurf(T,V(1,:),V(2,:),u,'EdgeColor','none')
%     caxis([0.199 0.201])
%     hold on
%     view(2)
%     hold off
%     axis('square')
%     title (sprintf('%dth solution',k))
%     set(gca,'FontSize',15)
%     colorbar
%     pause(0.01)
% end


