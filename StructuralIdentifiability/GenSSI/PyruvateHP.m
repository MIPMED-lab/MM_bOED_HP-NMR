function model = PyruvateHP()

    % Symbolic variables
    syms Php Xhp
    syms k_f1 k_r1 T1_P T1_X

    % Parameters
    model.sym.p = [k_f1; k_r1; T1_P; T1_X];

    % State variables
    model.sym.x = [Php; Xhp];

    % Control vectors (g)
    model.sym.g = [0
                   0
                   ];

    % Autonomous dynamics (f)
    model.sym.xdot = [k_r1*Xhp - ((1/T1_P)+k_f1)*Php
     k_f1*Php - ((1/T1_X)+k_r1)*Xhp];

    % Initial conditions
    model.sym.x0 = [240000;0;];
        
    % Observables
    model.sym.y = [Php
                   Xhp];
end