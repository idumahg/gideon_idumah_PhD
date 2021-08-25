function [Dt,Dr] = derivatives(t,N,U1,U2,U3,U4,C_transport,C_reaction3,C_reaction4,xi)
    % This function computes the derivatives of the transport flux and the
    % reaction flux used in the Jacobian.
    
    % It takes as input:
    % time (t), # of nodes (N), the concentrations (U1-U4), 
    % the constants for transport flux (C_transport) and for neuronal reaction (C_reaction3) and astrocytic reaction (C_reaction4). 
    % Lastly it also takes in xi, which is the proxy for the electrophysiology input for now.
    
    % It outputs the derivatives of transport fluxes (Dt) as a struct for
    % phi1-phi6. And also the derivatives of the reaction fluxes (Dr)
    
    %% derivative of transport flux
    T = C_transport(:,1); lamda = C_transport(:,2); M = C_transport(:,3); % constants for transport fluxes
    
    % for Phi1
    D_phi1_Glc = T(1).Glc*(M(1).Glc ./((U1.Glc + M(1).Glc).^2));
    D_phi1_O2 = lamda(1).oxy*ones(N,1);
    D_phi1_CO2 = lamda(1).CO2*ones(N,1);
    D_phi1_Lac = T(1).Lac*(M(1).Lac ./((U1.Lac + M(1).Lac).^2));

    % for Phi2
    D_phi2_Glc = T(2).Glc*(M(2).Glc ./((U2.Glc + M(2).Glc).^2));
    D_phi2_O2 = lamda(2).oxy*ones(N,1);
    D_phi2_CO2 = lamda(2).CO2*ones(N,1);
    D_phi2_Lac = T(2).Lac*(M(2).Lac ./((U2.Lac + M(2).Lac).^2));

    % for Phi3
    D_phi3_Glc = T(3).Glc*(M(3).Glc ./((U2.Glc + M(3).Glc).^2));
    D_phi3_O2 = lamda(3).oxy*ones(N,1);
    D_phi3_CO2 = lamda(3).CO2*ones(N,1);
    D_phi3_Lac = T(3).Lac*(M(3).Lac ./((U2.Lac + M(3).Lac).^2));
    D_phi3_Glu = zeros(N,1);
    D_phi3_Gln = T(3).Gln*(M(3).Gln ./((U2.Gln + M(3).Gln).^2));

    % for Phi4
    D_phi4_Glc = T(4).Glc*(M(4).Glc ./((U3.Glc + M(4).Glc).^2));
    D_phi4_O2 = lamda(4).oxy*ones(N,1);
    D_phi4_CO2 = lamda(4).CO2*ones(N,1);
    D_phi4_Lac = T(4).Lac*(M(4).Lac ./((U3.Lac + M(4).Lac).^2));
    D_phi4_Glu = T(4).Glu*(1 + xi(t))*( M(4).Glu ./((U3.Glu + M(4).Glu).^2));
    D_phi4_Gln = zeros(N,1);

    % for Phi5
    D_phi5_Glc = T(5).Glc*(M(5).Glc ./((U2.Glc + M(5).Glc).^2));
    D_phi5_O2 = lamda(5).oxy*ones(N,1);
    D_phi5_CO2 = lamda(5).CO2*ones(N,1);
    D_phi5_Lac = T(5).Lac*(M(5).Lac ./((U2.Lac + M(5).Lac).^2));
    D_phi5_Glu = T(5).Glu*(M(5).Glu ./((U2.Glu + M(5).Glu).^2));
    D_phi5_Gln = zeros(N,1);

    % for Phi6
    D_phi6_Glc = T(6).Glc*(M(6).Glc ./((U4.Glc + M(6).Glc).^2));
    D_phi6_O2 = lamda(6).oxy*ones(N,1);
    D_phi6_CO2 = lamda(6).CO2*ones(N,1);
    D_phi6_Lac = T(6).Lac*(M(6).Lac ./((U4.Lac + M(6).Lac).^2));
    D_phi6_Glu = zeros(N,1);
    D_phi6_Gln = T(6).Gln*(1 + xi(t))*(M(6).Gln ./((U4.Gln + M(6).Gln).^2));

    % struct of derivatives 
    Dt.phi1 = [D_phi1_Glc; D_phi1_O2; D_phi1_CO2; D_phi1_Lac];
    Dt.phi2 = [D_phi2_Glc; D_phi2_O2; D_phi2_CO2; D_phi2_Lac];
    Dt.phi3 = [D_phi3_Glc; D_phi3_O2; D_phi3_CO2; D_phi3_Lac; D_phi3_Glu; D_phi3_Gln];
    Dt.phi4 = [D_phi4_Glc; D_phi4_O2; D_phi4_CO2; D_phi4_Lac; D_phi4_Glu; D_phi4_Gln];
    Dt.phi5 = [D_phi5_Glc; D_phi5_O2; D_phi5_CO2; D_phi5_Lac; D_phi5_Glu; D_phi5_Gln];
    Dt.phi6 = [D_phi6_Glc; D_phi6_O2; D_phi6_CO2; D_phi6_Lac; D_phi6_Glu; D_phi6_Gln];
    
    %% derivative of reaction fluxes
    
    % for gamma3
    V3 = C_reaction3(1); mu3 = C_reaction3(2); K3 = C_reaction3(3); psi3 = C_reaction3(4); % constants for reaction flux (Neuron)
    
    p3 = U3.ATP./U3.ADP; % phosphorylation state for neuron
    p3_inv = U3.ADP./U3.ATP; % inverse of phosphorylation state for neuron
    r3 = U3.NADH./U3.NADplus; % redox state for neuron
    r3_inv = U3.NADplus./U3.NADH; % inverse of redox state for neuron
    
    g3_1.u1 = V3.Gcl * (r3_inv./(r3_inv + psi3.Gcl)).*(p3_inv./(p3_inv + mu3.Gcl)).*(K3.Gcl./((U3.Glc + K3.Gcl).^2));
    g3_1.u8 = V3.Gcl * (r3_inv./(r3_inv + psi3.Gcl)).*(-(mu3.Gcl./p3)./(U3.ATP.*((p3_inv + mu3.Gcl).^2))).*(U3.Glc./(U3.Glc + K3.Gcl));
    g3_1.u9 = V3.Gcl * (r3_inv./(r3_inv + psi3.Gcl)).*((mu3.Gcl)./(U3.ATP.*((p3_inv + mu3.Gcl).^2))).*(U3.Glc./(U3.Glc + K3.Gcl));
    g3_1.u10 = V3.Gcl * ((psi3.Gcl)./(U3.NADH.*((r3_inv + psi3.Gcl).^2))).*(p3_inv./(p3_inv + mu3.Gcl)).*(U3.Glc./(U3.Glc + K3.Gcl));
    g3_1.u11 = V3.Gcl * (-(psi3.Gcl./r3)./(U3.NADH.*((r3_inv + psi3.Gcl).^2))).*(p3_inv./(p3_inv + mu3.Gcl)).*(U3.Glc./(U3.Glc + K3.Gcl));
    
    g3_2.u7 = V3.LDH1 * (r3./(r3 + psi3.LDH1)).*(K3.LDH1./((U3.Pyr + K3.LDH1).^2));
    g3_2.u10 = V3.LDH1 * ((-psi3.LDH1*r3)./(U3.NADplus.*((r3 + psi3.LDH1).^2))).*(U3.Pyr./(U3.Pyr + K3.LDH1));
    g3_2.u11 = V3.LDH1 * (psi3.LDH1./(U3.NADplus.*((r3 + psi3.LDH1).^2))).*(U3.Pyr./(U3.Pyr + K3.LDH1));
    
    g3_3.u4 = V3.LDH2 * (r3_inv./(r3_inv + psi3.LDH2)).*(K3.LDH2./((U3.Lac + K3.LDH2).^2));
    g3_3.u10 = V3.LDH2 * (psi3.LDH2./(U3.NADH.*((r3_inv + psi3.LDH2).^2))).*(U3.Lac./(U3.Lac + K3.LDH2));
    g3_3.u11 = V3.LDH2 * ((-psi3.LDH2./r3)./(U3.NADH.*((r3_inv + psi3.LDH2).^2))).*(U3.Lac./(U3.Lac + K3.LDH2));
    
    g3_4.u7 = V3.TCA * (r3_inv./(r3_inv + psi3.TCA)).*(p3_inv./(p3_inv + mu3.TCA)).*(K3.TCA./((U3.Pyr + K3.TCA).^2));
    g3_4.u8 = V3.TCA * (r3_inv./(r3_inv + psi3.TCA)).*(-(mu3.TCA./p3)./(U3.ATP.*((p3_inv + mu3.TCA).^2))).*(U3.Pyr./(U3.Pyr + K3.TCA));
    g3_4.u9 = V3.TCA * (r3_inv./(r3_inv + psi3.TCA)).*((mu3.TCA)./(U3.ATP.*((p3_inv + mu3.TCA).^2))).*(U3.Pyr./(U3.Pyr + K3.TCA));
    g3_4.u10 = V3.TCA * ((psi3.TCA)./(U3.NADH.*((r3_inv + psi3.TCA).^2))).*(p3_inv./(p3_inv + mu3.TCA)).*(U3.Pyr./(U3.Pyr + K3.TCA));
    g3_4.u11 = V3.TCA * (-(psi3.TCA./r3)./(U3.NADH.*((r3_inv + psi3.TCA).^2))).*(p3_inv./(p3_inv + mu3.TCA)).*(U3.Pyr./(U3.Pyr + K3.TCA));
    
    g3_5.u2 = V3.OxPhos * (r3./(r3 + psi3.OxPhos)).*(p3_inv./(p3_inv + mu3.OxPhos)).*(K3.OxPhos./((U3.oxy + K3.OxPhos).^2));
    g3_5.u8 = V3.OxPhos * (r3./(r3 + psi3.OxPhos)).*(-(mu3.OxPhos./p3)./(U3.ATP.*((p3_inv + mu3.OxPhos).^2))).*(U3.oxy./(U3.oxy + K3.OxPhos));
    g3_5.u9 = V3.OxPhos * (r3./(r3 + psi3.OxPhos)).*((mu3.OxPhos)./(U3.ATP.*((p3_inv + mu3.OxPhos).^2))).*(U3.oxy./(U3.oxy + K3.OxPhos));
    g3_5.u10 = V3.OxPhos * ((-psi3.OxPhos*r3)./(U3.NADplus.*((r3 + psi3.OxPhos).^2))).*(p3_inv./(p3_inv + mu3.OxPhos)).*(U3.oxy./(U3.oxy + K3.OxPhos));
    g3_5.u11 = V3.OxPhos * (psi3.OxPhos./(U3.NADplus.*((r3 + psi3.OxPhos).^2))).*(p3_inv./(p3_inv + mu3.OxPhos)).*(U3.oxy./(U3.oxy + K3.OxPhos));
    
    g3_6.u6 = V3.PAG * (K3.PAG./((U3.Gln + K3.PAG).^2));
    
    g3_8.u8 = zeros(N,1);
    
    % for gamma4
    
    V4 = C_reaction4(1); mu4 = C_reaction4(2); K4 = C_reaction4(3); psi4 = C_reaction4(4); % constants for reaction flux (Astrocytes)
    
    p4 = U4.ATP./U4.ADP; % phosphorylation state for neuron
    p4_inv = U4.ADP./U4.ATP; % inverse of phosphorylation state for neuron
    r4 = U4.NADH./U4.NADplus; % redox state for neuron
    r4_inv = U4.NADplus./U4.NADH; % inverse of redox state for neuron

    g4_1.u1 = V4.Gcl * (r4_inv./(r4_inv + psi4.Gcl)).*(p4_inv./(p4_inv + mu4.Gcl)).*(K4.Gcl./((U4.Glc + K4.Gcl).^2));
    g4_1.u8 = V4.Gcl * (r4_inv./(r4_inv + psi4.Gcl)).*(-(mu4.Gcl./p4)./(U4.ATP.*((p4_inv + mu4.Gcl).^2))).*(U4.Glc./(U4.Glc + K4.Gcl));
    g4_1.u9 = V4.Gcl * (r4_inv./(r4_inv + psi4.Gcl)).*((mu4.Gcl)./(U4.ATP.*((p4_inv + mu4.Gcl).^2))).*(U4.Glc./(U4.Glc + K4.Gcl));
    g4_1.u10 = V4.Gcl * ((psi4.Gcl)./(U4.NADH.*((r4_inv + psi4.Gcl).^2))).*(p4_inv./(p4_inv + mu4.Gcl)).*(U4.Glc./(U4.Glc + K4.Gcl));
    g4_1.u11 = V4.Gcl * (-(psi4.Gcl./r4)./(U4.NADH.*((r4_inv + psi4.Gcl).^2))).*(p4_inv./(p4_inv + mu4.Gcl)).*(U4.Glc./(U4.Glc + K4.Gcl));

    g4_2.u7 = V4.LDH1 * (r4./(r4 + psi4.LDH1)).*(K4.LDH1./((U4.Pyr + K4.LDH1).^2));
    g4_2.u10 = V4.LDH1 * ((-psi4.LDH1*r4)./(U4.NADplus.*((r4 + psi4.LDH1).^2))).*(U4.Pyr./(U4.Pyr + K4.LDH1));
    g4_2.u11 = V4.LDH1 * (psi4.LDH1./(U4.NADplus.*((r4 + psi4.LDH1).^2))).*(U4.Pyr./(U4.Pyr + K4.LDH1));

    g4_3.u4 = V4.LDH2 * (r4_inv./(r4_inv + psi4.LDH2)).*(K4.LDH2./((U4.Lac + K4.LDH2).^2));
    g4_3.u10 = V4.LDH2 * (psi4.LDH2./(U4.NADH.*((r4_inv + psi4.LDH2).^2))).*(U4.Lac./(U4.Lac + K4.LDH2));
    g4_3.u11 = V4.LDH2 * ((-psi4.LDH2./r4)./(U4.NADH.*((r4_inv + psi4.LDH2).^2))).*(U4.Lac./(U4.Lac + K4.LDH2));

    g4_4.u7 = V4.TCA * (r4_inv./(r4_inv + psi4.TCA)).*(p4_inv./(p4_inv + mu4.TCA)).*(K4.TCA./((U4.Pyr + K4.TCA).^2));
    g4_4.u8 = V4.TCA * (r4_inv./(r4_inv + psi4.TCA)).*(-(mu4.TCA./p4)./(U4.ATP.*((p4_inv + mu4.TCA).^2))).*(U4.Pyr./(U4.Pyr + K4.TCA));
    g4_4.u9 = V4.TCA * (r4_inv./(r4_inv + psi4.TCA)).*((mu4.TCA)./(U4.ATP.*((p4_inv + mu4.TCA).^2))).*(U4.Pyr./(U4.Pyr + K4.TCA));
    g4_4.u10 = V4.TCA * ((psi4.TCA)./(U4.NADH.*((r4_inv + psi4.TCA).^2))).*(p4_inv./(p4_inv + mu4.TCA)).*(U4.Pyr./(U4.Pyr + K4.TCA));
    g4_4.u11 = V4.TCA * (-(psi4.TCA./r4)./(U4.NADH.*((r4_inv + psi4.TCA).^2))).*(p4_inv./(p4_inv + mu4.TCA)).*(U4.Pyr./(U4.Pyr + K4.TCA));

    g4_5.u2 = V4.OxPhos * (r4./(r4 + psi4.OxPhos)).*(p4_inv./(p4_inv + mu4.OxPhos)).*(K4.OxPhos./((U4.oxy + K4.OxPhos).^2));
    g4_5.u8 = V4.OxPhos * (r4./(r4 + psi4.OxPhos)).*(-(mu4.OxPhos./p4)./(U4.ATP.*((p4_inv + mu4.OxPhos).^2))).*(U4.oxy./(U4.oxy + K4.OxPhos));
    g4_5.u9 = V4.OxPhos * (r4./(r4 + psi4.OxPhos)).*((mu4.OxPhos)./(U4.ATP.*((p4_inv + mu4.OxPhos).^2))).*(U4.oxy./(U4.oxy + K4.OxPhos));
    g4_5.u10 = V4.OxPhos * ((-psi4.OxPhos*r4)./(U4.NADplus.*((r4 + psi4.OxPhos).^2))).*(p4_inv./(p4_inv + mu4.OxPhos)).*(U4.oxy./(U4.oxy + K4.OxPhos));
    g4_5.u11 = V4.OxPhos * (psi4.OxPhos./(U4.NADplus.*((r4 + psi4.OxPhos).^2))).*(p4_inv./(p4_inv + mu4.OxPhos)).*(U4.oxy./(U4.oxy + K4.OxPhos));

    g4_7.u5 = V4.GS * (p4./(p4 + mu4.GS)).*(K4.GS./((U4.Glu + K4.GS).^2));
    g4_7.u8 = V4.GS * (mu4.GS./(U4.ADP.*((p4 + mu4.GS).^2))).*(U4.Glu./(U4.Glu + K4.GS));
    g4_7.u9 = V4.GS *((-mu4.GS*p4)./(U4.ADP.*((p4 + mu4.GS).^2))).*(U4.Glu./(U4.Glu + K4.GS));


    g4_8.u8 = zeros(N,1);

    % struct of reaction derivatives
    Dr = [struct('gamma1', g3_1, 'gamma2', g3_2, 'gamma3', g3_3, 'gamma4', g3_4, 'gamma5', g3_5, 'gamma67', g3_6, 'gamma8', g3_8);
        struct('gamma1', g4_1, 'gamma2', g4_2, 'gamma3', g4_3, 'gamma4', g4_4, 'gamma5', g4_5, 'gamma67', g4_7, 'gamma8', g4_8)];
    
end

