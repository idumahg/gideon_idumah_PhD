function [u_sol,t_sol] = diffusion(N,BigK,BigM,BigB,BigU0,lamda,tspan,Psi,xi)

%% To DO LIST
% There was no u2_5 and u2_6 so i used value of u4_5 and u4_6.
% I used the value of Gln - Glu for GS reaction in neuron at constants.m
% I used 1's for the Glu-Gln transport constant in constants.m
% I used Gabriela's values for Phi_ATPase which is different from 2011
% I have taken away the blood in coupled system

aux = BigK + lamda*BigB;

u1_1 = 4.596*ones(N,1);
u1_2 = 7.009*ones(N,1);
u1_3 = 25.15*ones(N,1);
u1_4 = 1.194*ones(N,1);

u2_1 = 1.1289*ones(N,1);
u2_2 = 0.0165*ones(N,1);
u2_3 = 83.75*ones(N,1);
u2_4 = 1.247*ones(N,1);
u2_5 = 0.0103*ones(N,1);
u2_6 = 0.0100*ones(N,1);

u3_1 = 1.129*ones(N,1);
u3_2 = 0.0098*ones(N,1);
u3_3 = 85.46*ones(N,1);
u3_4 = 1.247*ones(N,1);
u3_5 = 14*ones(N,1);
u3_6 = 0.0010*ones(N,1);
u3_7 = 0.122*ones(N,1);
u3_8 = 2.189*ones(N,1);
u3_9 = 0.0226*ones(N,1);
u3_10 = 0.0299*ones(N,1);
u3_11 = 0.0031*ones(N,1);

u4_1 = 0.7536*ones(N,1);
u4_2 = 0.0098*ones(N,1);
u4_3 = 85.46*ones(N,1);
u4_4 = 1.2472*ones(N,1);
u4_5 = 0.0103*ones(N,1);
u4_6 = 0.0100*ones(N,1);
u4_7 = 0.1202*ones(N,1);
u4_8 = 2.189*ones(N,1);
u4_9 = 0.0227*ones(N,1);
u4_10 = 0.0299*ones(N,1);
u4_11 = 0.0031*ones(N,1);


%% initial value
U1_int = struct('Glc', u1_1, 'oxy', u1_2, 'CO2', u1_3, 'Lac', u1_4); % blood compartment

U2_int = struct('Glc', u2_1, 'oxy', u2_2, 'CO2', u2_3, 'Lac', u2_4, 'Glu', u2_5, 'Gln', u2_6); % ECS compartment

U3_int = struct('Glc', u3_1, 'oxy', u3_2, 'CO2', u3_3, 'Lac', u3_4, 'Glu', u3_5, 'Gln', u3_6, ...
    'Pyr', u3_7, 'ATP', u3_8, 'ADP', u3_9, 'NADplus', u3_10, 'NADH', u3_11); % Neuronal compartment

U4_int = struct('Glc', u4_1, 'oxy', u4_2, 'CO2', u4_3, 'Lac', u4_4, 'Glu', u4_5, 'Gln', u4_6, ...
    'Pyr', u4_7, 'ATP', u4_8, 'ADP', u4_9, 'NADplus', u4_10, 'NADH', u4_11); % Astrocyte compartment

BigU_int = [cell2mat(struct2cell(U1_int)); cell2mat(struct2cell(U2_int)); cell2mat(struct2cell(U3_int)); ...
    cell2mat(struct2cell(U4_int))]; %initial value of U.

% Transport flux coefficient matrix
BigF = sparse([-eye(4*N) eye(4*N) zeros(4*N,6*N) zeros(4*N,6*N) zeros(4*N,6*N) zeros(4*N,6*N);
    Eb(eye(4*N),N) -Eb(eye(4*N),N) -eye(6*N) eye(6*N) -eye(6*N) eye(6*N);
    zeros(11*N,4*N) zeros(11*N,4*N) EECS(eye(6*N),N) -EECS(eye(6*N),N) zeros(11*N,6*N) zeros(11*N,6*N) ;
    zeros(11*N,4*N) zeros(11*N,4*N) zeros(11*N,6*N) zeros(11*N,6*N) EECS(eye(6*N),N) -EECS(eye(6*N),N)]); 

% Reaction flux coefficient matrix
S3 = [-1 0 0 0 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 3 0 0 0; 0 1 -1 0 0 0 0;...
    0 0 0 0 0 1 0; 0 0 0 0 0 -1 0; 2 -1 1 -1 0 0 0; 2 0 0 1 5 0 -1; -2 0 0 -1 -5 0 1; ...
    -2 1 -1 -5 2 0 0; 2 -1 1 5 -2 0 0]; % Stoichiometry matrix for neuron

S4 = [-1 0 0 0 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 3 0 0 0; 0 1 -1 0 0 0 0;...
    0 0 0 0 0 -1 0; 0 0 0 0 0 1 0; 2 -1 1 -1 0 0 0; 2 0 0 1 5 -1 -1; -2 0 0 -1 -5 1 1; ...
    -2 1 -1 -5 2 0 0; 2 -1 1 5 -2 0 0]; % Stoichiometry matrix for astrocyte

Sigma3 = kron(S3, eye(N)); Sigma4 = kron(S4, eye(N));

BigV = sparse([zeros(4*N,7*N) zeros(4*N,7*N); zeros(6*N,7*N) zeros(6*N,7*N); ...
    Sigma3 zeros(11*N,7*N); zeros(11*N,7*N) Sigma4]);

%% ODE SOLVER
options = odeset('mass',Psi*BigM,'BDF','on','MaxOrder',2,'Jacobian',@(t,u) MyJacobian(t,u,N,aux,BigM,BigV,xi));
[t_sol,u_sol] = ode15s(@(t,u) rhs(t,u,aux,lamda,BigU0,BigM,BigV,BigF,N,xi), 0:0.4:tspan, BigU_int, options);


%% right hand side equation
function R = rhs(t,u,aux,lamda, BigU0, BigM, BigV, BigF, N, xi)
    
    [U1,U2,U3,U4] = unpack(u,N); % unpacks the array u into struct
    
%     [P] = blood(U1);
    
    [C_transport, C_reaction3, C_reaction4] = constants(); % gives the constant
    
    [Phi] = transport_flux(t,N,U1,U2,U3,U4,C_transport,xi); % computes the transport flux
   
    [gamma] = reaction_flux(t,N,U3,U4,C_reaction3,C_reaction4); % computes the reaction flux
    
    W = sparse(32*N,1); % Z = sparse(32*N,1);
    
%      Z(1:4*N) = P;

    W(4*N+1:10*N) = cell2mat(struct2cell(U2));
    W(21*N+1:end) = cell2mat(struct2cell(U4));
    
    R =  - aux*W + lamda*(BigU0) + BigM*BigV*gamma + BigM*BigF*Phi; 
end
% BigM*Z

%% Other functions
function out = Gb(A, N)
    out = [A, zeros(4*N, 2*N)];    
end

function out = Eb(A, N)
    out = [A; zeros(2*N, 4*N)];    
end

function out = GECS(A, N)
    out = [A, zeros(6*N, 5*N)];    
end

function out = EECS(A, N)
    out = [A; zeros(5*N, 6*N)];    
end

function out = EA(A, N)
    out = [A; zeros(5*N, 6*N)]; 
    out = [out, zeros(11*N, 5*N)];
end

function out = FA(A, N)
    out = [A; zeros(2*N, 4*N)]; 
    out = [out, zeros(6*N, 2*N)];
end


%% Jacobian matrix
function out = MyJacobian(t,u,N,aux,BigM,BigV,xi)
    
    C = 32;
    
    [U1,U2,U3,U4] = unpack(u,N); % unpacks the array u into struct
    
    [C_transport, C_reaction3, C_reaction4] = constants();
   
    [Dt,Dr] = derivatives(t,N,U1,U2,U3,U4,C_transport,C_reaction3,C_reaction4,xi); % computes the derivatives, Dt, Dr rep derivative of transport and reaction flux

    % for J1
    J1 = sparse(C*N, C*N);
    J1(:, 4*N+1:10*N) = aux(:, 4*N+1:10*N);
    J1(:, 21*N+1:end) = aux(:, 21*N+1:end);

    % for J2
    J2 = BigM * sparse([-diag(Dt.phi1) Gb(diag(Dt.phi2),N) zeros(4*N, 11*N) zeros(4*N, 11*N);...
        Eb(diag(Dt.phi1),N) -(FA(diag(Dt.phi2),N)+diag(Dt.phi3)+diag(Dt.phi5)) GECS(diag(Dt.phi4),N) GECS(diag(Dt.phi6),N); ...
        zeros(11*N, 4*N) EECS(diag(Dt.phi3),N) -EA(diag(Dt.phi4),N) zeros(11*N, 11*N); ...
        zeros(11*N, 4*N) EECS(diag(Dt.phi5),N) zeros(11*N, 11*N) -EA(diag(Dt.phi6),N)]);

    % for J3
    ON = zeros(N,N);
    J3_3 = [diag(Dr(1).gamma1.u1) ON ON ON ON ON ON diag(Dr(1).gamma1.u8) diag(Dr(1).gamma1.u9) diag(Dr(1).gamma1.u10) diag(Dr(1).gamma1.u11);...
        ON ON ON ON ON ON diag(Dr(1).gamma2.u7) ON ON diag(Dr(1).gamma2.u10) diag(Dr(1).gamma2.u11); ...
        ON ON ON diag(Dr(1).gamma3.u4) ON ON ON ON ON diag(Dr(1).gamma3.u10) diag(Dr(1).gamma3.u11); ...
        ON ON ON ON ON ON diag(Dr(1).gamma4.u7) diag(Dr(1).gamma4.u8) diag(Dr(1).gamma4.u9) diag(Dr(1).gamma4.u10) diag(Dr(1).gamma4.u11); ...
        ON diag(Dr(1).gamma5.u2) ON ON ON ON ON diag(Dr(1).gamma5.u8) diag(Dr(1).gamma5.u9) diag(Dr(1).gamma5.u10) diag(Dr(1).gamma5.u11); ...
        ON ON ON ON ON diag(Dr(1).gamma67.u6) ON ON ON ON ON; ...
        ON ON ON ON ON ON ON diag(Dr(1).gamma8.u8) ON ON ON];

    J3_4 = [diag(Dr(2).gamma1.u1) ON ON ON ON ON ON diag(Dr(2).gamma1.u8) diag(Dr(2).gamma1.u9) diag(Dr(2).gamma1.u10) diag(Dr(2).gamma1.u11);...
        ON ON ON ON ON ON diag(Dr(2).gamma2.u7) ON ON diag(Dr(2).gamma2.u10) diag(Dr(2).gamma2.u11); ...
        ON ON ON diag(Dr(2).gamma3.u4) ON ON ON ON ON diag(Dr(2).gamma3.u10) diag(Dr(2).gamma3.u11); ...
        ON ON ON ON ON ON diag(Dr(2).gamma4.u7) diag(Dr(2).gamma4.u8) diag(Dr(2).gamma4.u9) diag(Dr(2).gamma4.u10) diag(Dr(2).gamma4.u11); ...
        ON diag(Dr(2).gamma5.u2) ON ON ON ON ON diag(Dr(2).gamma5.u8) diag(Dr(2).gamma5.u9) diag(Dr(2).gamma5.u10) diag(Dr(2).gamma5.u11); ...
        ON ON ON ON diag(Dr(2).gamma67.u5) ON ON diag(Dr(2).gamma67.u8) diag(Dr(2).gamma67.u9) ON ON; ...
        ON ON ON ON ON ON ON diag(Dr(2).gamma8.u8) ON ON ON];

    J3 = BigM * BigV * sparse([zeros(7*N,4*N) zeros(7*N,6*N) J3_3 zeros(7*N,11*N);
        zeros(7*N,4*N) zeros(7*N,6*N) zeros(7*N,11*N) J3_4]);

    out = - J1 + J2 + J3;
end
        

end


