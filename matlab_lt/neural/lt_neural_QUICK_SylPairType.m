%% lt 3/15/18 - given 2 syllables, determine if is same type
function [issame, syl1, syl2] = lt_neural_QUICK_SylPairType(motif1, motif2)
% % e.g. 
% motif1 = 'ab(c)';
% motif2 = 'dd(c)x';
% % the comparison syl must be in parantheses
% output:
% issame = 1 if same, 0 if diff


ind1 = strfind(motif1, '(');
syl1 = motif1(ind1+1);

ind2 = strfind(motif2, '(');
syl2 = motif2(ind2+1);

assert(length(ind1) ==1, 'sadfasd');
assert(length(ind2) ==1, 'sdafd');

if syl1 == syl2
    issame=1;
else
    issame=0;
end