%  Multiple rules for the same head

% Probabilistic facts:
0.5::heads1.
0.6::heads2.

% Rules:
someHeads :- heads1.
someHeads :- heads2.

% If coin1 is heads → then “at least one coin is heads” is true.
% If coin2 is heads → then “at least one coin is heads” is true.
% If neither coin is heads, then someHeads won’t be true.

% Queries:
query(someHeads).

% what is the Probability of at least one of them landing on heads
% P(someHeads) = 1 - P(not heads1) * P(not heads2)
% = 1 - (1 - 0.5) * (1 - 0.6) = 1 - 0.5 * 0.4 = 1 - 0.2 = 0.8

