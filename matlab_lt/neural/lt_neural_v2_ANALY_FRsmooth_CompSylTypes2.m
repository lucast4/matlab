function [All_RhoTargSame, All_RhoTargDiff] = lt_neural_v2_ANALY_FRsmooth_CompSylTypes2(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, timewind, usediffFromBase, syltypesneeded, plotOn, collectDiffType, epochtoplot)
if ~exist('plotOn', 'var')
    plotOn = 1;
end

if ~exist('collectDiffType', 'var')
    collectDiffType=0;
end
%% ###################### LIMIT TO NEURONS THAT CONTAIN ALL SYL TYPES?

% =-==== go thru all switches. if bad then throw ou
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);

% ======================== SHUFFLE SYL TYPE? % within each neuron
indstokeep = [];
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
                % ==== count how many of each syl type there exists
                numtarg = sum(OUTDAT.All_istarg(indsthis)==1);
                numsame = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==1);
                numdiff = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==0);
                
                if all([numtarg numsame numdiff] >= syltypesneeded)
                    % then keep
%                     disp([numtarg numsame numdiff]);
%                     disp(indsthis)
                    indstokeep = [indstokeep; indsthis];
                end
                
            end
        end
    end
end

disp(['Keeping ' num2str(length(indstokeep)) '/' num2str(length(OUTDAT.All_birdnum)) ' datapoints, passes syltypes required criterion']);

OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);

%% look at change in FR duirng learning. is it more similar to target for sametype?


%%
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);


%% COLLECT DAT

All_RhoTargSame = {};
All_RhoTargDiff = {};
All_birdnum = [];
All_enum = [];
All_sw = [];
All_neur = [];
All_LearnZ_TargSame = {};
All_LearnZ_TargDiff= {};
All_LearnZ_targ = [];

% 
% % ==========
% AllAll_IsSame = [];
% AllAll_IsTarg = [];
% AllAll_birdnum = [];
% AllAll_enum = [];
% AllAll_sw = [];
% AllAll_neur = [];
% AllAll_Rho = [];

for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ================== ONLY DO IF LEARNING IN SAME DIRECTION FOR
            % BOTH TARGS
            if SwitchStruct.bird(i).exptnum(ii).switchlist(ss).targsAreSameSyl==0 | ...
                    length(unique([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}]))>1
                disp('skip!!');
                disp('NOTE: currently only works with one target - if multiple, need to modify so that each target is compared to its matched same-types')
                continue
            end
            
            % ==== what is learndir?
            learndir = unique([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}]);
            assert(length(learndir)==1);
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
%                 numtargs = sum(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
%                     & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==1);
%                 if onlyifonetarget==1
%                     if numtargs>1
%                         disp('skip, since >1 targ');
%                         continue
%                     end
%                     targind = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
%                     & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==1);
%  
%                 else
%                     disp('HAVE NOT CODED FOR CASES WITH >1 TARG...')
%                     return
%                 end
%                        
                
                % ======== go thru all syls and compare FR diff to that of
                % target syl.
                % -- 1) get mean across all syls (of deviation from
                % baseline)
                frmatAll = cell2mat(OUTDAT.AllOnlyMinusBase_FRsmooth(indsthis));
                
                % -- 2) for each syl, subtract its own from the mean.
                if usediffFromBase==0
                allmean = mean(frmatAll);
                frmatAll = frmatAll - allmean;
                end
                
                % -- 3) restrict to temporal window
                t = OUTDAT.AllMinusBase_tbinAll{indsthis(1)};
                indstokeep_t = t>=timewind(1) & t<timewind(2);
                frmatAll = frmatAll(:, indstokeep_t);
                
                
                % ########### [SAME TYPE]FOR ALL PAIRS COLLECT PAIRWISE CORRELATIONS
                istarg = OUTDAT.All_istarg(indsthis);
                issame = OUTDAT.All_issame(indsthis);
                pitchz = OUTDAT.AllMinusBase_PitchZ(indsthis, epochtoplot);
                
                
                indssub_targ = find(istarg==1);
                indssub_same = find(istarg==0 & issame==1);
                indssub_diff = find(istarg==0 & issame==0);
                assert(length(indssub_same)>0);
                
%                 disp('---');
%                 disp(indssub_same);
%                 disp(indssub_targ);
                
                rhoall = [];
                learnZ_targsame = [];
                for j=indssub_targ'
                    for jj=indssub_same'
                        
                        % ==== get corr for this targ and same type
                        rho = corr(frmatAll(j, :)', frmatAll(jj,:)');
                            
                        rhoall = [rhoall; rho];
                        
                        
                        % =============== COLLECT LEARNING
                        learnZ_targsame = [learnZ_targsame; pitchz(jj)*learndir];
                    end
                    All_LearnZ_targ = [All_LearnZ_targ; learndir*pitchz(j)];
                end
                
                % ================== OUTPUT
                All_RhoTargSame = [All_RhoTargSame; rhoall];
                All_LearnZ_TargSame = [All_LearnZ_TargSame; learnZ_targsame];
                
                
                % ########### [DIFF TYPE]FOR ALL PAIRS COLLECT PAIRWISE CORRELATIONS
                if collectDiffType ==1
                    rhoall = [];
                    learnZ_targdiff = [];
                    for j=indssub_targ'
                        for jj=indssub_diff'
                            
                            % ==== get corr for this targ and same type
                            rho = corr(frmatAll(j, :)', frmatAll(jj,:)');
                            
                            rhoall = [rhoall; rho];
                            
                            
                            % =============== COLLECT LEARNING
                            learnZ_targdiff = [learnZ_targdiff; pitchz(jj)*learndir];

                        end
                    end
                    
                    % ================== OUTPUT
                    All_RhoTargDiff = [All_RhoTargDiff; rhoall];
                    All_LearnZ_TargDiff = [All_LearnZ_TargDiff; learnZ_targdiff];
                end
                
                
                
               
                
%                 
%                 
%                 % ======== COLLECT ALL PAIRWISE CORRELATIONS RELATIVE TO
%                 % TARGET
%                 assert(sum(istarg)==1, 'have not yet coded for multiple targs');
%                 
%                 frtarg = frmatAll(istarg==1, :);
%                 rhoall = corr(frtarg', frmatAll');
%                 
%                 % --- extract all corr values based on syllable type
%                 rho_same = rhoall(istarg==0 & issame==1);
%                 rho_diff = rhoall(istarg==0 & issame==0);
%                 
%                 rhoTargSameDiff{1} = []; % in future, might put for other targs...
%                 rhoTargSameDiff{2} = rho_same;
%                 rhoTargSameDiff{3} = rho_diff;
%                 
                All_birdnum = [All_birdnum; i];
                All_enum = [All_enum; ii];
                All_sw = [All_sw; ss];
                All_neur = [All_neur; nn];
% 
% 
%                 % =============== COLLECT ALL DATA INDIVIDUALLY
%                 AllAll_IsSame = [AllAll_IsSame; issame];
%                 AllAll_IsTarg = [AllAll_IsTarg; istarg];
%                 AllAll_birdnum = [AllAll_birdnum; ones(size(issame))*i];
%                 AllAll_enum = [AllAll_enum; ones(size(issame))*ii];
%                 AllAll_sw = [AllAll_sw; ones(size(issame))*ss];
%                 AllAll_neur = [AllAll_neur; ones(size(issame))*nn];
%                 AllAll_Rho = [AllAll_Rho; rhoall'];

            end
        end
    end
end


%%
if plotOn==1
lt_figure;

% ========= 1) distribtuion
xcenters = -0.95:0.1:0.95;
lt_subplot(3,2,1); hold on;
xlabel('corr')
title('targ-same(k), targ-diff(r) [indiv syls]');
lt_plot_histogram(cell2mat(All_RhoTargSame), xcenters, 1, 1, [], 1, 'k');
lt_plot_histogram(cell2mat(All_RhoTargDiff), xcenters, 1, 1, [], 1, 'r');
lt_plot_zeroline_vert;


% ========== one value for each switch
indsgrp = lt_tools_grp2idx({All_birdnum, All_enum, All_sw});
rho_same = grpstats(cellfun(@mean, All_RhoTargSame), indsgrp, {'mean'});
rho_diff = grpstats(cellfun(@mean, All_RhoTargDiff), indsgrp, {'mean'});
lt_subplot(3,2,2); hold on;
xlabel('rho(targ-same');
ylabel('rho(targ-diff');
title('one value each switch');
plot(rho_same, rho_diff, 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% ========== one value for each expt
indsgrp = lt_tools_grp2idx({All_birdnum, All_enum});
rho_same = grpstats(cellfun(@mean, All_RhoTargSame), indsgrp, {'mean'});
rho_diff = grpstats(cellfun(@mean, All_RhoTargDiff), indsgrp, {'mean'});
lt_subplot(3,2,3); hold on;
xlabel('rho(targ-same');
ylabel('rho(targ-diff');
title('one value each expt');
plot(rho_same, rho_diff, 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% ======== SPLIT INTO HIGH AND LOW CORR - RLEATED TO LEARNING?
% --- one datapoint for each switch
indsgrp = lt_tools_grp2idx({All_birdnum, All_enum, All_sw});
rho_same = grpstats(cellfun(@mean, All_RhoTargSame), indsgrp, {'mean'});
rho_diff = grpstats(cellfun(@mean, All_RhoTargDiff), indsgrp, {'mean'});
learntarg = grpstats(All_LearnZ_targ, indsgrp, {'mean'});

learn_same = grpstats(cellfun(@mean, All_LearnZ_TargSame), indsgrp, {'mean'});
learn_diff = grpstats(cellfun(@nanmean, All_LearnZ_TargDiff), indsgrp, {'mean'});

% ======================== SAME TYUPE
lt_subplot(3,2,5); hold on;
title('high corr(r); low corr (b)');
xlabel('targ learn (z)');
ylabel('sametype learn (z)');

% --- high corr
indstmp = rho_same>median(rho_same);
pcol = 'r';
x = learntarg(indstmp);
y = learn_same(indstmp);
plot(x,y, 'o', 'Color', pcol);

% --- low corr
indstmp = rho_same<=median(rho_same);
pcol = 'b';
x = learntarg(indstmp);
y = learn_same(indstmp);
plot(x,y, 'o', 'Color', pcol);

lt_plot_makesquare_plot45line(gca, 'k')


% ======================= DIFF TYPE
lt_subplot(3,2,6); hold on;
title('high corr(r); low corr (b)');
xlabel('targ learn (z)');
ylabel('diff learn (z)');

% --- high corr
indstmp = rho_diff>median(rho_diff);
pcol = 'r';
x = learntarg(indstmp);
y = learn_diff(indstmp);
plot(x,y, 'o', 'Color', pcol);

% --- low corr
indstmp = rho_diff<=median(rho_diff);
pcol = 'b';
x = learntarg(indstmp);
y = learn_diff(indstmp);
plot(x,y, 'o', 'Color', pcol);

lt_plot_makesquare_plot45line(gca, 'k')



end

%% ======== OUTPUT

