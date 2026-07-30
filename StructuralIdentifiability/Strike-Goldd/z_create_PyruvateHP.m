

clear;

% 2 states
syms Php Xhp
x = [Php; Xhp;];

% 1 input
syms 
u = [];

% 5 parameters, 3 unknown
syms k_f1 k_r1 T1_P T1_X 
p =[k_f1; k_r1; T1_P; T1_X];

% 1 output
h = [Php; Xhp];

% initial conditions
syms Php Xhp
ics  = [Php; Xhp]; 
known_ics = [1, 1];

% dynamic equations
f = [k_r1*Xhp - ((1/T1_P)+k_f1)*Php;
     k_f1*Php - ((1/T1_X)+k_r1)*Xhp];

save('PyruvateHP','x','p','u','h','f','ics','known_ics');











