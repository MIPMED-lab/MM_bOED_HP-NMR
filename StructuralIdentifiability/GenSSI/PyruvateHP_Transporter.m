function model = PyruvateHP_Transporter()

    % Symbolic variables
    syms Php_OUT Php_IN Xhp
    syms k_in k_f1 k_r1 T1_P T1_X

    % Parameters
    model.sym.p = [k_in; k_f1; k_r1; T1_P; T1_X];

    % State variables
    model.sym.x = [Php_OUT; Php_IN; Xhp];

    % Control vectors (g)
    model.sym.g = [0
                   0
                   0
                   ];

    % Autonomous dynamics (f)
    model.sym.xdot = [-k_in*Php_OUT
     k_r1*Xhp - ((1/T1_P)+k_f1)*Php_IN + k_in*Php_OUT
     k_f1*Php_IN - ((1/T1_X)+k_r1)*Xhp];

    % Initial conditions
    model.sym.x0 = [240000;0;0;];
        
    % Observables
    model.sym.y = [Php_OUT+Php_IN
                   Xhp];
end