function [OUTSTRUCT, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, ...
    SwitchStruct, dataset)
%% lt 1/4/19 - does dedault filtering of outstruct (i.e. good experements)
% GOOD (default) = first switch of day, and shows some learning.

OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT); % so followniog code works...

if ~exist('dataset', 'var')
    dataset = '';
end
% can save presets identified by dataset.
% dataset = 'xcov_spikes', same as defgault, but excludse experiemnts without good spiking data

%% === filter each dimension independenyl
% if ~isempty(birdstoplot)
%     indstokeep = ismember(OUTSTRUCT.bnum, birdstoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end
%
% if ~isempty(expttoplot)
%     indstokeep = ismember(OUTSTRUCT.enum, expttoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end
%
% if ~isempty(swtoplot)
%     indstokeep = ismember(OUTSTRUCT.switch, swtoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end

%% == filter by specific types of switches.
if isempty(dataset)
    % --- to get specific switch types. ... [is done in addition to above
    % fitlers]
    % swtoget = {}; % passes if matches ANY of these
    swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
    firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    % swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
    % firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    firstswitchofday=1;
    indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
        firstswitchfortarget_withinday, firstswitchofday);
elseif strcmp(dataset, 'xcov_spikes')
    birdtoget = [1 2 4];
    swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
    firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    % swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
    % firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    firstswitchofday=1;
    indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
        firstswitchfortarget_withinday, firstswitchofday, birdtoget);
end

%%
if isfield(OUTSTRUCT, 'switch')
    indstokeep = ismember([OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.switch], indtoget_b_e_s, 'rows');
elseif isfield(OUTSTRUCT, 'swnum')
    indstokeep = ismember([OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.swnum], indtoget_b_e_s, 'rows');
end
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);

