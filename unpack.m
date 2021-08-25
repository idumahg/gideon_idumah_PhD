function [U1,U2,U3,U4] = unpack(u,N)
% this function converts an array into a struct. It takes as input the
% array u and the number of nodes N, and outputs struct U1-U4

    C = 32;
    A = mat2cell(u, N*ones(1,C));
    
    E1 = A(1:4); fields1 = {'Glc'; 'oxy'; 'CO2'; 'Lac'};
    E2 = A(5:10); fields2 = {'Glc'; 'oxy'; 'CO2'; 'Lac'; 'Glu'; 'Gln'};
    E3 = A(11:21); fields3 = {'Glc'; 'oxy'; 'CO2'; 'Lac'; 'Glu'; 'Gln'; 'Pyr'; 'ATP'; 'ADP'; 'NADplus'; 'NADH'};
    E4 = A(22:32); fields4 = {'Glc'; 'oxy'; 'CO2'; 'Lac'; 'Glu'; 'Gln'; 'Pyr'; 'ATP'; 'ADP'; 'NADplus'; 'NADH'};
    
    U1 = cell2struct(E1,fields1,1); 
    U2 = cell2struct(E2,fields2,1); 
    U3 = cell2struct(E3,fields3,1); 
    U4 = cell2struct(E4,fields4,1); 
end

