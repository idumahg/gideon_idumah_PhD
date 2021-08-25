function [C_transport, C_reaction3, C_reaction4] = constants()
% This function computes the struct for the constants for both transport
% flux, and reaction flux (neuronal and astrocytic).
% I first state the values of the constants and arrange them in a struct.

%% Data Source
% I used the values in Gabriela's thesis in Table 2.20 page 64 for reaction
% flux constants. I used the same values for lamda_O2 and lamda_CO2.

% I didnt find values for Glu and Gln, so I made that up.

% I used the values in Gabriela's thesis in Table 2.17 page 61 for
% transport flux constants. I used the values in Table 5 of 2011 paper for V3_PAG and V4_GS

%% transport flux constants
T1_Glc = 0.02; M1_Glc = 4.60; T1_Lac = 0.17; M1_Lac = 5.00; lamda1_O2 = 0.04; lamda1_CO2 = 0.04;
T2_Glc = T1_Glc; M2_Glc = M1_Glc; T2_Lac = T1_Lac; M2_Lac = M1_Lac; lamda2_O2 = lamda1_O2; lamda2_CO2 = lamda1_CO2;

T3_Glc = 83.33; M3_Glc = 5.00; T3_Lac = 66.67; M3_Lac = 0.40; lamda3_O2 = 0.94; lamda3_CO2 = 0.94;
T4_Glc = T3_Glc; M4_Glc = M3_Glc; T4_Lac = T3_Lac; M4_Lac = M3_Lac; lamda4_O2 = lamda3_O2; lamda4_CO2 = lamda3_CO2;
T3_Gln = 1; T4_Glu = 1; T5_Glu = 1; T6_Gln = 1; M3_Gln = 1; M4_Glu = 1; M5_Glu = 1; M6_Gln = 1; 

T5_Glc = 83.33; M5_Glc = 12500; T5_Lac = 66.67; M5_Lac = 0.40; lamda5_O2 = 0.68; lamda5_CO2 = 0.68;
T6_Glc = T5_Glc; M6_Glc = M5_Glc; T6_Lac = T5_Lac; M6_Lac = M5_Lac; lamda6_O2 = lamda5_O2; lamda6_CO2 = lamda5_CO2;


T = [struct('Glc', T1_Glc, 'oxy', '', 'CO2', '', 'Lac', T1_Lac, 'Glu', '', 'Gln', '');
    struct('Glc', T2_Glc, 'oxy', '', 'CO2', '', 'Lac', T2_Lac, 'Glu', '', 'Gln', '');
    struct('Glc', T3_Glc, 'oxy', '', 'CO2', '', 'Lac', T3_Lac, 'Glu', '', 'Gln', T3_Gln);
    struct('Glc', T4_Glc, 'oxy', '', 'CO2', '', 'Lac', T4_Lac, 'Glu', T4_Glu, 'Gln', '');
    struct('Glc', T5_Glc, 'oxy', '', 'CO2', '', 'Lac', T5_Lac, 'Glu', T5_Glu, 'Gln', '');
    struct('Glc', T6_Glc, 'oxy', '', 'CO2', '', 'Lac', T6_Lac, 'Glu', '', 'Gln', T6_Gln)];

lamda = [struct('Glc', '', 'oxy', lamda1_O2, 'CO2', lamda1_CO2, 'Lac', '', 'Glu', '', 'Gln', '');
    struct('Glc', '', 'oxy', lamda2_O2, 'CO2', lamda2_CO2, 'Lac', '', 'Glu', '', 'Gln', '');
    struct('Glc', '', 'oxy', lamda3_O2, 'CO2', lamda3_CO2, 'Lac', '', 'Glu', '', 'Gln', '');
    struct('Glc', '', 'oxy', lamda4_O2, 'CO2', lamda4_CO2, 'Lac', '', 'Glu', '', 'Gln', '');
    struct('Glc', '', 'oxy', lamda5_O2, 'CO2', lamda5_CO2, 'Lac', '', 'Glu', '', 'Gln', '');
    struct('Glc', '', 'oxy', lamda6_O2, 'CO2', lamda6_CO2, 'Lac', '', 'Glu', '', 'Gln', '')];

M = [struct('Glc', M1_Glc, 'oxy', '', 'CO2', '', 'Lac', M1_Lac, 'Glu', '', 'Gln', '');
    struct('Glc', M2_Glc, 'oxy', '', 'CO2', '', 'Lac', M2_Lac, 'Glu', '', 'Gln', '');
    struct('Glc', M3_Glc, 'oxy', '', 'CO2', '', 'Lac', M3_Lac, 'Glu', '', 'Gln', M3_Gln);
    struct('Glc', M4_Glc, 'oxy', '', 'CO2', '', 'Lac', M4_Lac, 'Glu', M4_Glu, 'Gln', '');
    struct('Glc', M5_Glc, 'oxy', '', 'CO2', '', 'Lac', M5_Lac, 'Glu', M5_Glu, 'Gln', '');
    struct('Glc', M6_Glc, 'oxy', '', 'CO2', '', 'Lac', M6_Lac, 'Glu', '', 'Gln', M6_Gln)];
    
C_transport = [T, lamda, M];

%% reaction flux constants

% Neuron 
V3_Gcl = 0.26; K3_Gcl = 4.60; mu3_Gcl = 0.09; psi3_Gcl = 10.00;
V3_LDH1 = 1436; K3_LDH1 = 2.15; psi3_LDH1 = 0.10;
V3_LDH2 = 1579.83; K3_LDH2 = 23.70; psi3_LDH2 = 10.00;
V3_TCA = 0.03; K3_TCA = 0.01; mu3_TCA = 0.01; psi3_TCA = 10.00;
V3_OxPhos = 8.18; K3_OxPhos = 1.00; mu3_OxPhos = 0.01; psi3_OxPhos = 0.10;
V3_PAG = 1.181; K3_PAG = 0.003;


V3 = struct('Gcl', V3_Gcl, 'LDH1', V3_LDH1, 'LDH2', V3_LDH2, 'TCA', V3_TCA, 'OxPhos', V3_OxPhos, 'PAG', V3_PAG);
mu3 = struct('Gcl', mu3_Gcl, 'LDH1', '', 'LDH2', '', 'TCA', mu3_TCA, 'OxPhos', mu3_OxPhos, 'PAG', '');
K3 = struct('Gcl', K3_Gcl, 'LDH1', K3_LDH1, 'LDH2', K3_LDH2, 'TCA', K3_TCA, 'OxPhos', K3_OxPhos, 'PAG', K3_PAG); 
psi3 = struct('Gcl', psi3_Gcl, 'LDH1', psi3_LDH1, 'LDH2', psi3_LDH2, 'TCA', psi3_TCA, 'OxPhos', psi3_OxPhos, 'PAG', '');

% Astrocyte
V4_Gcl = 0.25; K4_Gcl = 3.10; mu4_Gcl = 0.09; psi4_Gcl = 10.00;
V4_LDH1 = 4160; K4_LDH1 = 6.24; psi4_LDH1 = 0.10;
V4_LDH2 = 3245; K4_LDH2 = 48.66; psi4_LDH2 = 10.00;
V4_TCA = 0.01; K4_TCA = 0.01; mu4_TCA = 0.01; psi4_TCA = 10.00;
V4_OxPhos = 2.55; K4_OxPhos = 1.00; mu4_OxPhos = 0.01; psi4_OxPhos = 0.10;
V4_GS = 2.3614; mu4_GS = 100; K4_GS = 0.03;

V4 = struct('Gcl', V4_Gcl, 'LDH1', V4_LDH1, 'LDH2', V4_LDH2, 'TCA', V4_TCA, 'OxPhos', V4_OxPhos, 'GS', V4_GS);
mu4 = struct('Gcl', mu4_Gcl, 'LDH1', '', 'LDH2', '', 'TCA', mu4_TCA, 'OxPhos', mu4_OxPhos, 'GS', mu4_GS);
K4 = struct('Gcl', K4_Gcl, 'LDH1', K4_LDH1, 'LDH2', K4_LDH2, 'TCA', K4_TCA, 'OxPhos', K4_OxPhos, 'GS', K4_GS); 
psi4 = struct('Gcl', psi4_Gcl, 'LDH1', psi4_LDH1, 'LDH2', psi4_LDH2, 'TCA', psi4_TCA, 'OxPhos', psi4_OxPhos, 'GS', '');

C_reaction3 = [V3,mu3,K3,psi3];
C_reaction4 = [V4,mu4,K4,psi4];

end

