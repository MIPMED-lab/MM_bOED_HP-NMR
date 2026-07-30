%--------------------------------------------------------------------------
% File that creates the Pitavastatin model. 
% It stores it in a mat-file named Pitavastatin.mat.
% The model is described in eqs. (103)-(105) in the paper:
%--------------------------------------------------------------------------
% Grandjean TR, Chappell MJ, Yates JW, Evans ND. (2014) Structural 
% identifiability analyses of candidate models for in vitro Pitavastatin
% hepatic uptake.
% Comput Methods Programs Biomed. 2014;114(3):e60-e69
%--------------------------------------------------------------------------
clear

syms Pout Php Xhp NADH NAD LDH LDHac LDHacE LDHna LDHP XhpHyper PhpHyper...                % states
     ...                    % known constants
     T1_X T1_P kin kpl kbn kun kbp kup kunE kbnE ki kr       % unknown parameters

% states:
x    = [Pout Php Xhp NADH NAD LDH LDHac LDHacE LDHna LDHP XhpHyper PhpHyper].';   

% output:
h    = [XhpHyper PhpHyper];   

% no input:
u    = [];

% parameters:
p    = [T1_X T1_P kin kpl kbn kun kbp kup kunE kbnE ki kr].';

% dynamic equations:
f    = [ - (Pout*kin);        
         (Pout*kin) + kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php;
         kpl*LDHP;
         kun*LDHac - kbn*LDH*NADH;
         kunE*LDHacE - kbnE*LDH*NAD;
         kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD;
         kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php;
         kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE;
         ki*LDH*Php - kr*LDHna;
         kbp*LDHac*Php - kup*LDHP - kpl*LDHP;
         LDH - (XhpHyper/T1_X);
         (Pout+Php+LDHna+LDHP) - PhpHyper/T1_P];  
         


% initial conditions: 
ics = [];

% which initial conditions are known:
known_ics = [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1]; 

save('PyruvateHP_CompetitiveRepression','x','h','u','p','f','ics','known_ics');




