function [P] = blood(U1)
% this function computes the values for the blood compartment for Glucose,
% oxygen, CO2 and lactate.
C_art_Glc = 5;
C_art_Lac = 1.1;
C_art_O2 = 9.14;
C_art_CO2 = 23;


Q = 0.50;
F = 2/3;


b.Glc = (Q/F)*(C_art_Glc - U1.Glc);
b.O2 = (Q/F)*(XTotalToFree2(C_art_O2 - U1.O2, 2.5));
b.CO2 = (Q/F)*(XTotalToFree2(C_art_CO2 - U1.CO2, 1));
b.Lac = (Q/F)*(C_art_Lac - U1.Lac);

P = cell2mat(struct2cell(b));
end

