function [gamma] = reaction_flux(t,N,U3,U4,C_reaction3,C_reaction4)

    % this function computes the reaction flux (gamma). 
    
    % It takes as input:
    % time (t), # of nodes (N), the concentrations (U3-U4), 
    % the constants for reaction flux (C_reaction3 and C_reaction4). 
    
    %% Using Gabriela's values for gamma_ATPase:
    H1 = 4.3; s = 0.23; H2 = 0.833*H1; Jpump_act = 0.4444; Jpump_base = 0.0811; Jglia_act = 0.1933; Jglia_base = 0.1897;
    eta2 = 0.26; eta3 = 0.4;
    
    %% Neuron (gamma3)

    V3 = C_reaction3(1); mu3 = C_reaction3(2); K3 = C_reaction3(3); psi3 = C_reaction3(4); % constants for transport fluxes
    
    p3 = U3.ATP./U3.ADP; % phosphorylation state for neuron
    p3_inv = U3.ADP./U3.ATP; % inverse of phosphorylation state for neuron
    r3 = U3.NADH./U3.NADplus; % redox state for neuron
    r3_inv = U3.NADplus./U3.NADH; % inverse of redox state for neuron

    gamma3_Gcl = V3.Gcl * (r3_inv./(r3_inv + psi3.Gcl)).*(p3_inv./(p3_inv + mu3.Gcl)).*(U3.Glc./(U3.Glc + K3.Gcl));
    gamma3_LDH1 = V3.LDH1 * (r3./(r3 + psi3.LDH1)).*(U3.Pyr./(U3.Pyr + K3.LDH1));
    gamma3_LDH2 = V3.LDH2 * (r3_inv./(r3_inv + psi3.LDH2)).*(U3.Lac./(U3.Lac + K3.LDH2));
    gamma3_TCA = V3.TCA * (r3_inv./(r3_inv + psi3.TCA)).*(p3_inv./(p3_inv + mu3.TCA)).*(U3.Pyr./(U3.Pyr + K3.TCA));
    gamma3_OxPhos = V3.OxPhos * (r3./(r3 + psi3.OxPhos)).*(p3_inv./(p3_inv + mu3.OxPhos)).*(U3.oxy./(U3.oxy + K3.OxPhos));
    gamma_PAG = V3.PAG * (U3.Gln./(U3.Gln + K3.PAG));
    
    if t >=2 && t<=5
        gamma3_ATPase = (H1 + s*eta3*Jpump_act)*ones(N,1);
    else
        gamma3_ATPase = (H1 + s*eta3*Jpump_base)*ones(N,1);
    end

    %% Astrocyte (gamma4)

    V4 = C_reaction4(1); mu4 = C_reaction4(2); K4 = C_reaction4(3); psi4 = C_reaction4(4);
    
    p4 = U4.ATP./U4.ADP; % phosphorylation state for astrocyte
    p4_inv = U4.ADP./U4.ATP; % inverse of phosphorylation state for neuron
    r4 = U4.NADH./U4.NADplus; % redox state for astrocyte
    r4_inv = U4.NADplus./U4.NADH; % inverse of redox state for neuron
    
    gamma4_Gcl = V4.Gcl * (r4_inv./(r4_inv + psi4.Gcl)).*(p4_inv./(p4_inv + mu4.Gcl)).*(U4.Glc./(U4.Glc + K4.Gcl));
    gamma4_LDH1 = V4.LDH1 * (r4./(r4 + psi4.LDH1)).*(U4.Pyr./(U4.Pyr + K4.LDH1));
    gamma4_LDH2 = V4.LDH2 * (r4_inv./(r4_inv + psi4.LDH2)).*(U4.Lac./(U4.Lac + K4.LDH2));
    gamma4_TCA = V4.TCA * (r4_inv./(r4_inv + psi4.TCA)).*(p4_inv./(p4_inv + mu4.TCA)).*(U4.Pyr./(U4.Pyr + K4.TCA));
    gamma4_OxPhos = V4.OxPhos * (r4./(r4 + psi4.OxPhos)).*(p4_inv./(p4_inv + mu4.OxPhos)).*(U4.oxy./(U4.oxy + K4.OxPhos));
    gamma_GS = V4.GS * (p4./(p4 + mu4.GS)).*(U4.Glu./(U4.Glu + K4.GS));
    
    if t >=2 && t<=5
        gamma4_ATPase = (H2 + s*0.5*eta2*Jglia_act)*ones(N,1);
    else
        gamma4_ATPase = (H2 + s*0.5*eta2*Jglia_base)*ones(N,1);
    end

    %% Phi reaction fluxes
    gamma3 = [gamma3_Gcl; gamma3_LDH1; gamma3_LDH2; gamma3_TCA; gamma3_OxPhos; gamma_PAG; gamma3_ATPase]; % for neuron
    gamma4 = [gamma4_Gcl; gamma4_LDH1; gamma4_LDH2; gamma4_TCA; gamma4_OxPhos; gamma_GS; gamma4_ATPase]; % for astrocytes
    
    gamma = [gamma3; gamma4];
end

