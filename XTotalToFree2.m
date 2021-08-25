% Conversion from total oxygen concentration to free oxygen concentration,
% assuming that the total concentration follows the Hill model
%==========================================================================

function Cfree = XTotalToFree2(y,nHill)

Hct     = 0.45;     % Hematocrit
CHbrbc  = 5.18;     % Concentration of hemoglobin in RBC
sigma   = 1.4e-3;   % Oxygen solubility mM/mm Hg (Observe: Keener-Sneyd 
                    % gives sigma 1.4e-6 Molar/mm HG
KHill   = 26*sigma; %
%nHill   = 2.5;      %

beta = 4*Hct*CHbrbc;

Free2Tot = @(p) p + beta*(p^nHill)/(KHill^nHill + p^nHill);

Clow =  0;
Chigh = 0.25;

tol = 1e-10;
dy = 1;
while dy > tol
    Cmid = 0.5*(Clow + Chigh);
    ymid = Free2Tot(Cmid);
    if ymid<y
        Clow = Cmid;
    else
        Chigh = Cmid;
    end
    dy = abs(y - ymid);
end

Cfree = Cmid;



