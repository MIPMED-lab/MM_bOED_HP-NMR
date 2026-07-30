genssiAskForConfirmation(8);

% Symbolic parameters for identifiability analysis
syms kin kf kr kpl T1_P T1_L

% Structural identifiability analysis (for a subset of the parameters)
genssiMain('PyruvateHP_MichaelisMenten',12,[kin; kf; kr; kpl; T1_P; T1_L]);