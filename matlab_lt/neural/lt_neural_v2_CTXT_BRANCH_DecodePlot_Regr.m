function lme = lt_neural_v2_CTXT_BRANCH_DecodePlot_Regr(DatAll, subtractshuffmean, ...
    UseDecodeShuff, LMANonly, combineContexts)
numctxts_categorical = 1; % if 1 then will not be linear relative to num contexts
% LMANonly = 1;
% maxctxts = 4;
%% lt 5/6/18

if ~exist('combineContexts', 'var');
    combineContexts=[];
end

%%
% ================= ORGANIZE DATA
% if numctxts_categorical ==1
%     Nctxts = categorical(DatAll.AllNumCtxts);
% else
%     Nctxts = DatAll.AllNumCtxts;
% end

if subtractshuffmean==1
    Decode = DatAll.AllDecode' - DatAll.AllDecode_neg_mean;
else
    Decode = DatAll.AllDecode';
end

%%
if UseDecodeShuff==1
    Decode = DatAll.AllDecode_neg_mean;
end

%%
tbl = table(Decode, DatAll.AllNumCtxts, DatAll.AllBirdNum, ...
    DatAll.AllBrainRegion, DatAll.AllBranchNum, DatAll.AllNeurNum, ...
    'VariableNames', {'Decode', 'NumCtxts', 'Bird', 'Bregion', 'Branch', 'Neuron'});

%% model
% mdl = 'Decode ~ NumCtxts + (NumCtxts|Bird) + (NumCtxts|Bird:Branch) + (NumCtxts|Bird:Neuron)';
if combineContexts==1
mdl = 'Decode ~ 1 + (1|Bird) + (1|Bird:Branch)';    
else
mdl = 'Decode ~ NumCtxts + (NumCtxts|Bird) + (NumCtxts|Bird:Branch)';
end
% mdl = 'Decode ~ NumCtxts';

%% run
% ================ LMAN ONLY
% if LMANonly==1
%     inds = strcmp(tbl.Bregion, 'LMAN') & tbl.NumCtxts<maxctxts;
% end
if LMANonly==1
    inds = strcmp(tbl.Bregion, 'LMAN');
end

% ===============
if numctxts_categorical==1
    tbl.NumCtxts = categorical(tbl.NumCtxts);
end

lme = fitlme(tbl(inds,:), mdl);
