

clear;

% 2 states
syms Php_OUT Php_IN Xhp P X
x = [Php_OUT; Php_IN; Xhp; P; X];

% 1 input
syms 
u = [];

% 5 parameters, 3 unknown
syms k_in k_f1 k_r1 T1_P T1_X 
p =[k_in; k_f1; k_r1; T1_P; T1_X];

% 1 output
h = [Php_OUT+Php_IN; Xhp];

% initial conditions
syms Php_OUT Php_IN Xhp
ics  = [Php_OUT,Php_IN; Xhp]; 
known_ics = [1, 1, 1];

% dynamic equations
f = [-k_in*Php_OUT;
     k_r1*Xhp - ((1/T1_P)+k_f1)*Php_IN + k_in*Php_OUT;
     k_f1*Php_IN - ((1/T1_X)+k_r1)*Xhp];

save('PyruvateHP_Transporter','x','p','u','h','f','ics','known_ics');











