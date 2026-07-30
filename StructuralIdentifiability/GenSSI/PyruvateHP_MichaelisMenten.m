function model = PyruvateHP_MichaelisMenten()

    % Symbolic variables
    syms Pout Pin LDH LDHPyr Lac
    syms kin kf kr kpl T1_P T1_L

    % Parameters
    model.sym.p = [Pout; Pin; LDH; LDHPyr; Lac];

    % State variables
    model.sym.x = [kin; kf; kr; kpl; T1_P; T1_L];

    % Control vectors (g)
    model.sym.g = [0
                   0
                   0
                   0
                   0
                   ];

    % Autonomous dynamics (f)
    model.sym.xdot = [-Pout*kin - Pout/T1_P
         Pout*kin + kr*LDHPyr - kf*Pin*LDH - Pin/T1_P
         kr*LDHPyr + kpl*LDHPyr - kf*Pin*LDH
         kf*Pin*LDH - kr*LDHPyr - kpl*LDHPyr - LDHPyr/T1_P
         kpl*LDHPyr - (Lac/T1_L)];

    % Initial conditions
    model.sym.x0 = [240000;0;100;0;0];
        
    % Observables
    model.sym.y = [Pout+Pin+LDHPyr 
                    Lac];
end