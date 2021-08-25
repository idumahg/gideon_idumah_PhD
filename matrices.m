function [V,M,K,B,U] = matrices(n,kappa,U_b)
% calculates the stiffness matrix K, the mass matrix M, the boundary vector
% U, the boundary mass matrix B, the vector of nodes V.

% It takes as input:
% the number of circles n, diffusion tensor kappa, the concentration at the bdy U_b, 

j = 0:n; r = j/n;
nj = round((2*pi*r(2:end))./(r(2:end) - r(1:end-1))); % # of intervals on each of the circles
r = r(2:end);

V = zeros(2,1);

for j = 1:n
    k = 0:nj(j)-1;
    theta_rand = ((2*pi*k)/nj(j)) + rand;
    x_rand = r(j)*cos(theta_rand);
    y_rand = r(j)*sin(theta_rand);
    V = [V, [x_rand;y_rand]];
end

% % to visualize
% figure(1)
% plot(V(1,:),V(2,:),'k.','MarkerSize',1)
% hold on
% axis([-1,1,-1,1])
% axis('square')
% set(gca,'FontSize',15)

% Find the topology matrix
T = delaunay(V(1,:)',V(2,:)');
L = size(T,1);

% % check counterclockwise order
% for l = 1:L
%     I = T(l,:);
%     aux = [V(:, I(2)) - V(:, I(1)), V(:,I(3)) - V(:,I(1))];
%     if det(aux) < 0
%         I([2,3]) = I([3,2]);
%         T(l,:) = I; % swap the columns in T
%     end
%     
% end

% for ell = 1:L    
%     x = [V(1,T(ell,:)),V(1,T(ell,1))];    
%     y = [V(2,T(ell,:)),V(2,T(ell,1))];     
%     plot(x,y,'k-','LineWidth',0.5)    
% end
% hold off

[E,L] = GetEdgeMatrix(T); % to get the edge matrix

Nv = size(V,2); % original size of node vector

% Add new nodes
for l = 1:size(E,1)
    I = E(l,:);
    x_new = 0.5*sum(V(:,I),2);
    V = [V,x_new];
end


N = size(V,2);
Nt = size(T,1);

% add new columns to the topology matrix
for l = 1:Nt
    T_ext(l,:) = [T(l,:), Nv + L(l,:)];
end
T = T_ext;


% where V is a 2-by-N matrix where each column contains the indices of
% each nodal point and N is the number of nodal points.

% T is a Nt-by-3 matrix where each row contains the indices of each
% triangle and Nt is the number of triangles.


K = sparse(N, N);
M = sparse(N, N);
B = sparse(N, N);
U = sparse(N, 1);

    
% three point quadrature rule for reference triangle
w = [1/3 1/3 1/3];
p = [1/6 1/6 2/3; 1/6 2/3 1/6];

Phi = []; 
for j = 1:3
    Z = shape(p(:,j));
    Phi = [Phi, Z(:,1)];
    G{j} = Z(:,2:3);
end

% to account for round off errors
tol = 1e-8;
bdy0 = find(vecnorm(V) <= 1+tol); % finds the bdy nodes for circular domain
bdy1 = find(vecnorm(V) >= 1-tol); % finds the bdy nodes for circular domain
bdy = intersect(bdy0, bdy1);


% to find the boundary edges
J_bdy = [];
for j = 1:size(E,1)
    s = ismember(E(j,:),bdy); % checks if any of the nodes of Edge j is in bdy
    if nnz(s) == 2 % both nodes lie on bdy
        J_bdy = [J_bdy, j];% add such edge j.
    end
end

% to calculate boundary mass matrix and boundary vector.
for l = 1:length(J_bdy)
        
    % find end points of the edge
    i = E(J_bdy(l), 1);
    k = Nv + J_bdy(l); % midpoint node
    j = E(J_bdy(l), 2);
    
    el = norm(V(:,i) - V(:,j));

    U_loc = U_b*[el*(1/6);el*(2/3); el*(1/6)];
    B_loc = (el/15)*[2 1 -0.5; 1 8 1; -0.5 1 2];
    
    B([i,k,j], [i,k,j]) = B([i,k,j], [i,k,j]) + B_loc;
    U([i,k,j]) = U([i,k,j]) + U_loc; 
end

% to calculate mass matrix and stiffness matrix
for l = 1:Nt
     I = T(l, :); 
    
    % vertices and edges of the element
    z = V(:,I(1:3));
        
    % Jacobian matrix
    J = [z(:,2) - z(:,1) z(:,3) - z(:,1)];
    
    A1 = J'\G{1}';
    A2 = J'\G{2}';
    A3 = J'\G{3}';
    
    % Local stiffness matrix
    K_loc = (0.5)*abs(det(J))*(w(1)*A1'*kappa*A1 + w(2)*A2'*kappa*A2 + w(3)*A3'*kappa*A3);
    
    % Local mass matrix
    W = diag([w(1), w(2), w(3)]);
    M_loc = 0.5*abs(det(J))*(Phi*W*Phi');
    
    
    M(I,I) = M(I,I) + M_loc;
    K(I,I) = K(I,I) + K_loc;
end

% calculates the six shape functions for second order scheme
function Z = shape(z)
t = z(1); b= z(2);
aux = [1 -3 -3 2 4 2; 0 -1 0 2 0 0; 0 0 -1 0 0 2; 0 4 0 -4 -4 0; 0 0 0 0 4 0; 0 0 4 0 -4 -4];
aux2 = [ 1 0 0; t 1 0; b 0 1; t^2 2*t 0; t*b b t; b^2 0 2*b];
Z = aux*aux2;
end
end

