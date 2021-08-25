function [E,TE] = GetEdgeMatrix(T)
%
% Given the topology matrix, the progam returns the edges and the edge
% topology. The output is organized so that in the edge topology matrix TE,
% the first edge is between the first and second vertex, the second edge 
% between the second and third vertex, and the last one between the third
% and first.
%
% Input:  T  - (nt,3) topology matrix
% Output: E  - (ne,2) edges, each row giving the number of vertices defining 
%              the edge (orientation not important)
%         TE - (nt,3) edge topology, each row giving the number of edges of
%              the elements
% -------------------------------------------------------------------------
% E Somersalo
%--------------------------------------------------------------------------

nt = size(T,1);   % Number of elements in the mesh
nv = max(max(T)); % Number of nodes in the mesh
ne = nv+nt-1;     % Number of edges from Euler characteristics

E  = zeros(ne,2);
TE = zeros(nt,3);

% Initialize with the first element

Tj = T(1,:);
E(1,:) = [Tj(1),Tj(2)];
E(2,:) = [Tj(2),Tj(3)];
E(3,:) = [Tj(1),Tj(3)];
TE(1,:) = [1,2,3];

count   = 3;   % Number of edges found so far

for j = 2:nt
    % Scan through the elements
    Tj = T(j,:);
    % Candidate edges
    e1 = [Tj(1),Tj(2)];
    e2 = [Tj(2),Tj(3)];
    e3 = [Tj(1),Tj(3)];
    % Check if e1 are already in the edge list
    I1  = find(min(e1) == min(E(1:count,:)'));
    I2  = find(max(e1) == max(E(1:count,:)'));
    I12 = intersect(I1,I2);
    if length(I12) == 0
        % e1 is a new edge; add it to the edge list and update TE
        count = count + 1;
        E(count,:) = e1;
        TE(j,1) = count;
    else
        % e1 is already in the list. Just update TE
        TE(j,1) = I12;
    end
    % Check if e2 are already in the edge list
    I1  = find(min(e2) == min(E(1:count,:)'));
    I2  = find(max(e2) == max(E(1:count,:)'));
    I12 = intersect(I1,I2);
    if length(I12) == 0
        % e2 is a new edge; add it to the edge list and update TE
        count = count + 1;
        E(count,:) = e2;
        TE(j,2) = count;
    else
        % e1 is already in the list. Just update TE
        TE(j,2) = I12;
    end
    % Check if e3 are already in the edge list
    I1  = find(min(e3) == min(E(1:count,:)'));
    I2  = find(max(e3) == max(E(1:count,:)'));
    I12 = intersect(I1,I2);
    if length(I12) == 0
        % e3 is a new edge; add it to the edge list and update TE
        count = count + 1;
        E(count,:) = e3;
        TE(j,3) = count;
    else
        % e1 is already in the list. Just update TE
        TE(j,3) = I12;
    end

end

