%This is test script for the conjunctive and disjunctive distance based combination operators

% Some input mass functions 
m1 = [0 0 0 0.3 0 0 0.5 0.2];
m2 = [0 0 0.3 0 0 0 0.4 0.3];

% Conjunctive combination of m1 and m2 (commonality case)
conjQP_multi([m1 ; m2],'q')
% Conjunctive combination of m1 and m2 (plausibility case)
conjQP_multi([m1 ; m2],'pl')



% Disjunctive combination of m1 and m2 (commonality case)
disjQP_multi([m1 ; m2],'q')
% Disjunctive combination of m1 and m2 (plausibility case)
disjQP_multi([m1 ; m2],'pl')
