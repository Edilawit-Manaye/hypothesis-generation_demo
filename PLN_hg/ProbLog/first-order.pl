% Probabilistic facts:
0.6::lands_heads(_).

% Each coin has a 60% chance of landing heads.


% Background information:
coin(c1).
coin(c2).
coin(c3).
coin(c4).



% Rules:
heads(C) :- coin(C), lands_heads(C).
someHeads :- heads(_). 

% someHeads is true if there exists at least one coin C such that heads(C) is true.
% At least one of the four coins (c1, c2, c3, c4) lands heads.

% For any object C, if C is a coin and C lands heads, then we say "C lands heads".





% Queries:
query(someHeads).



% steps to evaluate the query
% 1. someHeads :- heads(_). 
% 2. heads(C) :- coin(C), lands_heads(C).
% 3. coin(c1). coin(c2). coin(c3). coin(c4) = True
% 4. heads(c1) :- lands_heads(c1).
%    heads(c2) :- lands_heads(c2).
%    heads(c3) :- lands_heads(c3).
% P(heads(ci))=P(landsh​eads(ci))=0.6  and 0.4 for tails
% 5. P(someHeads) = 1 - P(not heads(c1)) * P(not heads(c2)) * P(not heads(c3)) * P(not heads(c4))