% Probabilistic facts:
0.5::heads1.
0.6::heads2.

% Rules:
twoHeads :- heads1, heads2.

% Queries:
query(heads1).
query(heads2).
query(twoHeads).

% P(twoHeads)=P(heads1∧heads2)
% =P(heads1)*P(heads2)=0.5*0.6=0.3