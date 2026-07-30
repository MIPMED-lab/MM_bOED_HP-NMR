genssiAskForConfirmation(8);

% Symbolic parameters for identifiability analysis
syms k_f1 k_r1 T1_P T1_X

% Structural identifiability analysis (for a subset of the parameters)
genssiMain('PyruvateHP',12,[k_f1; k_r1; T1_P; T1_X]);