function [Phi] = transport_flux(t,N,U1,U2,U3,U4,C_transport,xi)
    % this function computes the transport flux (Phi). 
    
    % It takes as input:
    % time (t), # of nodes (N), the concentrations (U1-U4), 
    % the constants for transport flux (C_transport) and xi, which is the proxy for the electrophysiology input for now.
    
    %%
    T = C_transport(:,1); lamda = C_transport(:,2); M = C_transport(:,3); % constants for transport fluxes
    
    % for Phi1
    phi1_Glc = T(1).Glc*(U1.Glc./(U1.Glc + M(1).Glc));
    phi1_O2 = lamda(1).oxy.*(U1.oxy);
    phi1_CO2 = lamda(1).CO2*U1.CO2;
    phi1_Lac = T(1).Lac*(U1.Lac./(U1.Lac + M(1).Lac));

    % for Phi2
    phi2_Glc = T(2).Glc*(U2.Glc./(U2.Glc + M(2).Glc));
    phi2_O2 = lamda(2).oxy*U2.oxy;
    phi2_CO2 = lamda(2).CO2*U2.CO2;
    phi2_Lac = T(2).Lac*(U2.Lac./(U2.Lac + M(2).Lac));

    % for Phi3
    phi3_Glc = T(3).Glc*(U2.Glc./(U2.Glc + M(3).Glc));
    phi3_O2 = lamda(3).oxy*U2.oxy;
    phi3_CO2 = lamda(3).CO2*U2.CO2;
    phi3_Lac = T(3).Lac*(U2.Lac./(U2.Lac + M(3).Lac));
    phi3_Glu = zeros(N,1);
    phi3_Gln = T(3).Gln*(U2.Gln./(U2.Gln + M(3).Gln));

    % for Phi4
    phi4_Glc = T(4).Glc*(U3.Glc./(U3.Glc + M(4).Glc));
    phi4_O2 = lamda(4).oxy*U3.oxy;
    phi4_CO2 = lamda(4).CO2*U3.CO2;
    phi4_Lac = T(4).Lac*(U3.Lac./(U3.Lac + M(4).Lac));
    phi4_Glu = (T(4).Glu)*(1 + xi(t))*(U3.Glu./(U3.Glu + M(4).Glu));
    phi4_Gln = zeros(N,1);

    % for Phi5
    phi5_Glc = T(5).Glc*(U2.Glc./(U2.Glc + M(5).Glc));
    phi5_O2 = lamda(5).oxy*U2.oxy;
    phi5_CO2 = lamda(5).CO2*U2.CO2;
    phi5_Lac = T(5).Lac*(U2.Lac./(U2.Lac + M(5).Lac));
    phi5_Glu = T(5).Glu*(U2.Glu./(U2.Glu + M(5).Glu));
    phi5_Gln = zeros(N,1);

    % for Phi6
    phi6_Glc = T(6).Glc*(U4.Glc./(U4.Glc + M(6).Glc));
    phi6_O2 = lamda(6).oxy*U4.oxy;
    phi6_CO2 = lamda(6).CO2*U4.CO2;
    phi6_Lac = T(6).Lac*(U4.Lac./(U4.Lac + M(6).Lac));
    phi6_Glu = zeros(N,1);
    phi6_Gln = T(6).Gln*(1 + xi(t))*(U4.Gln./(U4.Gln + M(6).Gln));

    % transport fluxes
    phi1 = [phi1_Glc; phi1_O2; phi1_CO2; phi1_Lac];
    phi2 = [phi2_Glc; phi2_O2; phi2_CO2; phi2_Lac];
    phi3 = [phi3_Glc; phi3_O2; phi3_CO2; phi3_Lac; phi3_Glu; phi3_Gln];
    phi4 = [phi4_Glc; phi4_O2; phi4_CO2; phi4_Lac; phi4_Glu; phi4_Gln];
    phi5 = [phi5_Glc; phi5_O2; phi5_CO2; phi5_Lac; phi5_Glu; phi5_Gln];
    phi6 = [phi6_Glc; phi6_O2; phi6_CO2; phi6_Lac; phi6_Glu; phi6_Gln];

    Phi = [phi1; phi2; phi3; phi4; phi5; phi6];
    
end

