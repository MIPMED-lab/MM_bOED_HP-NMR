genssiAskForConfirmation(8);

% Symbolic parameters for identifiability analysis
syms T1_X T1_P kin kpl kbn kun kbp kup kunE kbnE ki kr

% Structural identifiability analysis (for a subset of the parameters)
genssiMain('PyruvateHP_CompetitiveRepression',12,[T1_X; T1_P; kin; kpl; kbn; kun; kbp; kup; kunE; kbnE; ki; kr]);