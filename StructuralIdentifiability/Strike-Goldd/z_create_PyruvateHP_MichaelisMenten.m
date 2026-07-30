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

syms Pout Pin LDH LDHPyr Lac...                % states
     ...                    % known constants
     kin kf kr kpl T1_P T1_L       % unknown parameters

% states:
x    = [Pout Pin LDH LDHPyr Lac].';   

% output:
h    = [Pout+Pin+LDHPyr Lac];   

% no input:
u    = [];

% parameters:
p    = [kin kf kr kpl T1_P T1_L].';

% dynamic equations:
f    = [ -Pout*kin - Pout/T1_P;
         Pout*kin + kr*LDHPyr - kf*Pin*LDH - Pin/T1_P;
         kr*LDHPyr + kpl*LDHPyr - kf*Pin*LDH;
         kf*Pin*LDH - kr*LDHPyr - kpl*LDHPyr - LDHPyr/T1_P;
         kpl*LDHPyr - (Lac/T1_L)];  
         


% initial conditions: 
ics = [];

% which initial conditions are known:
known_ics = [1,1,0,1,1]; 

save('PyruvateHP_MichaelisMenten','x','h','u','p','f','ics','known_ics');




