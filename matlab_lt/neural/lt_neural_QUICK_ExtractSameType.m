function [SameTypeSyls, DiffTypeSyls, SingleSyls] = ...
    lt_neural_QUICK_ExtractSameType(motif_regexpr_str, targsyl)

% % e.g.
% motif_regexpr_str = {'ab(b)', 'c(d)'}; 
% targsyl = 'c(d)';

%% lt - give it motif list and target syl. it tells you which are same type

        [SameTypeSyls, DiffTypeSyls, ~, SingleSyls] = ...
            lt_neural_v2_extractSameType(motif_regexpr_str, {targsyl});
