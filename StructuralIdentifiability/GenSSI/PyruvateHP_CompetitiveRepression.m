function model = PyruvateHP_CompetitiveRepression()

    % Symbolic variables
    syms Pout Php Xhp NADH NAD LDH LDHac LDHacE LDHna LDHP XhpHyper PhpHyper
    syms T1_X T1_P kin kpl kbn kun kbp kup kunE kbnE ki kr

    % Parameters
    model.sym.p = [T1_X; T1_P; kin; kpl; kbn; kun; kbp; kup; kunE; kbnE; ki; kr];

    % State variables
    model.sym.x = [Pout; Php; Xhp; NADH; NAD; LDH; LDHac; LDHacE; LDHna; LDHP; XhpHyper; PhpHyper];

    % Control vectors (g)
    model.sym.g = [0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   ];

    % Autonomous dynamics (f)
    model.sym.xdot = [- (Pout*kin)       
         (Pout*kin) + kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php
         kpl*LDHP
         kun*LDHac - kbn*LDH*NADH
         kunE*LDHacE - kbnE*LDH*NAD
         kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD
         kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php
         kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE
         ki*LDH*Php - kr*LDHna
         kbp*LDHac*Php - kup*LDHP - kpl*LDHP
         LDH - (XhpHyper/T1_X)
         (Pout+Php+LDHna+LDHP) - PhpHyper/T1_P];

    % Initial conditions
    model.sym.x0 = [1; 1; 1; 0; 0; 0; 1; 1; 1; 1; 1; 1];
        
    % Observables
    model.sym.y = [XhpHyper 
                    PhpHyper];
end