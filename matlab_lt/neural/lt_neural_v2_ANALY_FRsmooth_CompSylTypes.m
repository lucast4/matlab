function lt_neural_v2_ANALY_FRsmooth_CompSylTypes(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, timewind, onlyifonetarget, usediffFromBase, syltypesneeded)

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
                    disp([numtarg numsame numdiff]);
                    disp(indsthis)
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

All_RhoTargSameDiff = {};
All_birdnum = [];
All_enum = [];
All_sw = [];
All_neur = [];

% ==========
AllAll_IsSame = [];
AllAll_IsTarg = [];
AllAll_birdnum = [];
AllAll_enum = [];
AllAll_sw = [];
AllAll_neur = [];
AllAll_Rho = [];
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            rhoTargSameDiff = cell(1,3);
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
                numtargs = sum(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==1);
                if onlyifonetarget==1
                    if numtargs>1
                        disp('skip, since >1 targ');
                        continue
                    end
                    targind = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==1);
 
                else
                    disp('HAVE NOT CODED FOR CASES WITH >1 TARG...')
                    return
                end
                       
                
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
                
                
                % ======== COLLECT ALL PAIRWISE CORRELATIONS RELATIVE TO
                % TARGET
                istarg = OUTDAT.All_istarg(indsthis);
                issame = OUTDAT.All_issame(indsthis);
                assert(sum(istarg)==1, 'have not yet coded for multiple targs');
                
                frtarg = frmatAll(istarg==1, :);
                rhoall = corr(frtarg', frmatAll');
                
                % --- extract all corr values based on syllable type
                rho_same = rhoall(istarg==0 & issame==1);
                rho_diff = rhoall(istarg==0 & issame==0);
                
                rhoTargSameDiff{1} = []; % in future, might put for other targs...
                rhoTargSameDiff{2} = rho_same;
                rhoTargSameDiff{3} = rho_diff;
                
                % ================== OUTPUT
                All_RhoTargSameDiff = [All_RhoTargSameDiff; rhoTargSameDiff];
                All_birdnum = [All_birdnum; i];
                All_enum = [All_enum; ii];
                All_sw = [All_sw; ss];
                All_neur = [All_neur; nn];


                % =============== COLLECT ALL DATA INDIVIDUALLY
                AllAll_IsSame = [AllAll_IsSame; issame];
                AllAll_IsTarg = [AllAll_IsTarg; istarg];
                AllAll_birdnum = [AllAll_birdnum; ones(size(issame))*i];
                AllAll_enum = [AllAll_enum; ones(size(issame))*ii];
                AllAll_sw = [AllAll_sw; ones(size(issame))*ss];
                AllAll_neur = [AllAll_neur; ones(size(issame))*nn];
                AllAll_Rho = [AllAll_Rho; rhoall'];

            end
        end
    end
end

%% ============== PLOT

figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ============================ IGNORE TARGET FOR NOPW
if onlyifonetarget==1
    All_RhoTargSameDiff = All_RhoTargSameDiff(:,[2 3]);
    indsame = 1;
    inddiff = 2;
else
    indsame = 2;
    inddiff = 3;
    
end

% ==== 1) plot separetely for each expt/switch



% ==== 2) overall distributions (each expt take mean for each type)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('SAME -- DIFF');
title('all neur (only if have same and diff');
ylabel('corr with targ');

% lt_tools_grp2idx({All_birdnum, All_enum, All_sw, All_neur});
% 1) keep only those with both same and diff data.
indsgood = ~any(cellfun(@isempty, All_RhoTargSameDiff)'); % has both same and diff
Y = All_RhoTargSameDiff(indsgood, :);
% 2) get means
allrho_means = cell2mat(cellfun(@mean, Y, 'UniformOutput', 0));

for j=1:size(allrho_means,1)
   x = 1:size(allrho_means,2);
   y = allrho_means(j,:);
   plot(x, y, '-ok');
end

lt_plot(x+0.2, mean(allrho_means,1), {'Errors', lt_sem(allrho_means), 'Color', 'r'});


% --- significance
p = signrank(allrho_means(:,indsame), allrho_means(:,inddiff));
lt_plot_pvalue(p, 'srank', 1);

% ---
xlim([0 4]);
lt_plot_zeroline;


% ===== 3) overall distribution, for all motifs (not taking mean)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('SAME -- DIFF');
title('all neur and motifs (only if have same and diff');
ylabel('corr with targ');

% -- 1) just those with both smae and diff
indsgood = ~any(cellfun(@isempty, All_RhoTargSameDiff)'); % has both same and diff
Y = All_RhoTargSameDiff(indsgood, :);

ysame = cell2mat(cellfun(@transpose, Y(:,indsame), 'UniformOutput', 0));
ydiff = cell2mat(cellfun(@transpose, Y(:,inddiff), 'UniformOutput', 0));

lt_plot_MultDist({ysame, ydiff}, [1 2], 1, 'k');
lt_plot_zeroline;



% ========== LME
AllAll_idxcount = lt_tools_grp2idx({AllAll_birdnum, AllAll_enum, AllAll_sw, AllAll_neur});
% AllAll_idxcount = lt_tools_grp2idx({AllAll_birdnum, AllAll_enum, AllAll_sw});
AllAll_Rho = double(AllAll_Rho);
tbl = table(AllAll_Rho, AllAll_idxcount, AllAll_IsSame);

% --- remove all targs
tbl = tbl(AllAll_IsTarg==0, :);

mdl = 'AllAll_Rho ~ AllAll_IsSame + (AllAll_IsSame|AllAll_idxcount)';
lme = fitlme(tbl, mdl);















