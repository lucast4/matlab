function [SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct, PARAMS, ...
    twind, fwind)
%% lt 11/11/18 - scalar correaltion between coherence and pitch


%% 

inds_t = PARAMS.tbins>=twind(1) & PARAMS.tbins<=twind(2);
inds_f = PARAMS.ffbins>=fwind(1) & PARAMS.ffbins<=fwind(2);


%% EXTRACT SCALAR VALUE OF COHERENCE

numbirds = length(SwitchCohStruct.bird);
for i=1:numbirds
    numexpt = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpt
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        for ss=1:numswitch
           nummotif = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
           for mm=1:nummotif
               disp([num2str(i) '-' num2str(ii) '-' num2str(ss)]);
               datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm);
               if isempty(datthis.chanpair)
                   continue
               end
               
               % ==================== extract coherene scalars for this dat
               fnamethis = [datthis.fileprefix '/Coh' datthis.filesuffix];
               pairstoget = datthis.chanpairstokeep;
               
               cohmat = lt_neural_LFP_QUICK_loadProcessDat(fnamethis, pairstoget);
                
               % ------------ COLLECT SCALAR
               cohscal = squeeze(nanmean(nanmean(cohmat(inds_t, inds_f, :,:),1),2));
               
               % ------------ SAVE OUTPUT
               SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).cohscalar = single(cohscal);
                              
           end
        end
    end
end


%%  save params
PARAMS.cohscal_twind = twind;
PARAMS.cohscal_fwind = fwind;

