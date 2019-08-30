%% lt 5/25/19 - to remove or73, since histology revealed that is not really confidently LMAN. 
% and neural activity is not good SNR.

bthis = find(strcmp({SummaryStruct.birds.birdname}, 'or74bk35'));

indsgood = OUTSTRUCT.All_birdnum~=bthis;

% === hold on to ptypes.
ptypes = OUTSTRUCT.PairTypesInOrder;
OUTSTRUCT = rmfield(OUTSTRUCT, 'PairTypesInOrder');

% === only keep good inds
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indsgood, 1);
OUTSTRUCT.PairTypesInOrder = ptypes;
disp('DONE!');