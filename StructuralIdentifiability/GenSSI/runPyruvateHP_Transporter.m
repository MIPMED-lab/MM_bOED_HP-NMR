genssiAskForConfirmation(8);

% Symbolic parameters for identifiability analysis
syms k_in k_f1 k_r1 T1_P T1_X

% Structural identifiability analysis (for a subset of the parameters)
genssiMain('PyruvateHP_Transporter',12,[k_in; k_f1; k_r1; T1_P; T1_X]);