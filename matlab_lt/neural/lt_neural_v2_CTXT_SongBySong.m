function OUTSTRUCT = lt_neural_v2_CTXT_SongBySong(analyfname, bregiontoplot, ...
    motifcorr_dirty, ignoreMelBirds)
%%

% === only same motif (dirty version)
% dirty: i.e. figures out if one same motif by making sure is exactly same
% number of renditions per song
% motifcorr_dirty =0; % keep at 0. don't have enough data to do 1. if 0, then takes average
% % over song.
%
% ignoreMelBirds = 1; % since I do not know song times for those birds

Nmin = 7; % min sample size (minimum of 2 contexts);

%%
apos =1;

%% load branch
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

% load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/ALLBRANCHv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
end

%% ======= Filter data (e.g. remove noise, poor labels, etc)

Params.LocationsToKeep = {};
Params.birdstoexclude = {};
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below

ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%% use first time bin only (can modify)

tt = 1;

%% ========= COMPUILE BY ACTUAL BRANCH IDS (based on regexp str)
DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, ...
    analyfname, tt, 0, [], prms);

disp('MANY EXTRACTION FAILURES BECUASE BRANCH POINTS DOESNT EXIST!! - IS OK')


%% ========= determine premotor windows

frinds = 1:1000*(prms.motifpredur + prms.motifpostdur)-1;

%%

AllClassnum = {};
AllPitch = {};
AllFRbinned = {};
AllFRmeans = {};
AllSongInds = {};
AllSongDnums = {};

AllBirdnum = [];
AllBranchnum = [];
AllNeurnum = [];
AllNeurnum_fakeInd = [];

numbirds = length(DATSTRUCT_BYBRANCH.bird);

numsongsall = [];
rendspersongall = [];

for i =1:numbirds
    
    if ignoreMelBirds==1
        % -------- is this mel bird?
        if isfield(ALLBRANCH.SummaryStruct.birds(i).neurons(1), 'isRAsobermel')
            % ------ if so, remove
            continue
        end
        
        
    end
    numbranch = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    for bb=1:numbranch
        
        DAT = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT;
        indsneur = find(strcmp(DAT.brainregion, bregiontoplot));
        
        
        for nn=indsneur'
            
            nclass = length(DAT.frdat(nn).classnum);
            
            
            % ================ pitch and class and frmat
            pitchall = [];
            classall = [];
            frmatall = [];
            tsongall = [];
            for cc=1:nclass
                ff = DAT.FFstruct(nn).classnum(cc).t_ff(:,2)';
                tdnum = DAT.FFstruct(nn).classnum(cc).t_ff(:,1)';
                frmat = DAT.frdat(nn).classnum(cc).frmat(frinds,:); % take frinds, so same timebins.
                
                tsongall = [tsongall tdnum];
                pitchall = [pitchall ff];
                classall = [classall cc*ones(size(ff))];
                frmatall = [frmatall frmat]; % fr over time bins
            end
            
            % ------ take premotor window for frmat (for frmat, determine which time inds to take (i.e.
            % premotor window)
            premwind = DAT.PREMOTORDECODE_struct(nn).window_relonset;
            xwind = prms.motifpredur + premwind; % window, in sec
            
            % ---- if window too short, lengthen to 10ms.
            if xwind(2)-xwind(1)<0.01
                xwind(2) = xwind(1)+0.01;
            end
            
            xtimes = lt_neural_QUICK_XfromFRmat(frmatall); % time of bins
            
            indsgood = xtimes>=xwind(1) & xtimes<xwind(2);
            frmatall = frmatall(indsgood, :);
            xtimes = xtimes(indsgood);
            
            % ------ convert frmatall to fr bins
            TrimDown = 1;
            frbinsize = 0.008;
            [frbinned, xtimes] = lt_neural_v2_QUICK_binFR(frmatall', xtimes, frbinsize, TrimDown);
            
            % ----- get one value for mean fr in window (i.e. spike count)
            frmeans = mean(frmatall,1);
            
            % --------- CONVERT FROM SONG TIMES TO JUST INDICES
            tsonginds = grp2idx(tsongall);
            
            % ======================= OUTPUT
            AllSongInds = [AllSongInds; tsonginds];
            AllSongDnums = [AllSongDnums; tsongall'];
            AllClassnum = [AllClassnum; classall];
            AllPitch = [AllPitch; pitchall];
            AllFRbinned = [AllFRbinned; frbinned];
            AllFRmeans = [AllFRmeans; frmeans];
            AllBirdnum = [AllBirdnum; i];
            AllBranchnum = [AllBranchnum; bb];
            
            AllNeurnum = [AllNeurnum; DAT.IndsOriginal_neuron(nn)];
            AllNeurnum_fakeInd = [AllNeurnum_fakeInd; nn];
            
            % ======================= sanity check
            % number of songs and renditions per song
            %             numsongsall = [numsongsall; length(unique(tsongall))];
            
            
        end
    end
end


%%

Maxbird = max(AllBirdnum);
Maxbranch = max(AllBranchnum);
Maxneur = max(AllNeurnum);

%% get correlations (song by song)

% ------------ for sanity cheeck
RendRhoAll = [];


% ============= for collect
AllPair_pitchrho = [];
AllPair_fraterho = [];

for i=1:Maxbird
    for ii=1:Maxbranch
        for nn=1:Maxneur
            
            indsthis = AllBirdnum==i & AllBranchnum==ii & AllNeurnum==nn;
            
            if ~any(indsthis)
                continue
            end
            
            assert(sum(indsthis)==1); % should only be one branch point
            
            % ====================== iterate thru all classes and collect
            classvals = AllClassnum{indsthis};
            tsongvals = AllSongDnums{indsthis};
            pitchvals = AllPitch{indsthis};
            fratevals = AllFRmeans{indsthis};
            
            if all(isnan(pitchvals))
                continue
            end
            
            songlist = unique(tsongvals);
            classlist = unique(AllClassnum{indsthis});
            assert(max(classlist) == length(classlist), 'skipped classes? need match to use class as ind in cell');
            
            pitch_songclass_cell = cell(length(songlist), length(classlist)); % song by class
            frate_songclass_cell = cell(length(songlist), length(classlist)); % song by class
            
            if any(isnan(songlist))
                keyboard
            end
            
            
            for cc = classlist
                
                for ss = 1:length(songlist)
                    
                    songthis = songlist(ss);
                    % ================ extract dat
                    indstmp = classvals==cc & tsongvals'==songthis;
                    
                    pitch_songclass_cell{ss, cc} = pitchvals(indstmp);
                    frate_songclass_cell{ss, cc} = fratevals(indstmp);
                    
                    
                    
                end
            end
            
            % ###################### 2 VERSIONS
            
            if motifcorr_dirty==1
                % %%%%%%%%%%%%%%%%%%%%%%%% 1) ONLY MOTIFS (DIRTY),
                % IS EASIER AS DOES NOT NEED TO
                % FIGURE OUT BEST WAY TO AVERAGE OVER SONG BEFORE
                % CORRELATIONS (FOR INSTANCE) - DIRTY MEANS THAT FIGURES
                % OUT IF MOTIFS BY ONLY LOOKING AT CASES WITH PERFECT
                % CORRELATION IN NUMBER OF CASES PER SONG
                
                nrends_songclass_mat = cellfun(@numel, pitch_songclass_cell);
                
                % ========= go thru all pairs of classes
                for c=1:length(classlist)
                    
                    for cc=c+1:length(classlist)
                        
                        class1 = classlist(c);
                        class2 = classlist(cc);
                        
                        rho = corr(nrends_songclass_mat(:,class1), nrends_songclass_mat(:,class2));
                        
                        RendRhoAll = [RendRhoAll; rho];
                        
                        % ============= NOTE: STOPPPED HERE. SINCE
                        % BASICALLY NO CASES FOR RA THAT HAVE ON SAME
                        % MOTIF, SINCE THIS IS CONSTRAINED TO XAA. SHOULD
                        % DO USING ALL MOTIFS. TO CONTINUE WOULD RESTRICT
                        % TO CASES WIOTH HIGH RHO, AND DO RENDITION BY
                        % RENDITION CORRELATION.
                        
                    end
                end
                
            else
                % %%%%%%%%%%%%%%%%%%%%%% 2) TAKE MEAN FOR AECH SONG, THEN
                % GET CORRELATION... (ASSUMES A SHARED LATENT VARIABLE FOR
                % PITCH FOR EACH SONG)
                pitchmat = cellfun(@mean, pitch_songclass_cell);
                fratemat = cellfun(@mean, frate_songclass_cell);
                
                % ====== go thru all pairs
                for c=1:length(classlist)
                    
                    for cc=c+1:length(classlist)
                        
                        class1 = classlist(c);
                        class2 = classlist(cc);
                        
                        indstmp = all(~isnan(pitchmat(:, [class1 class2]))');
                        if sum(indstmp)<Nmin
                            continue
                        end
                        
                        pitchcorr = corr(pitchmat(indstmp, class1), pitchmat(indstmp, class2));
                        fratecorr = corr(fratemat(indstmp, class1), fratemat(indstmp, class2));
                        
                        % ============ COOLLECT
                        AllPair_pitchrho = [AllPair_pitchrho; pitchcorr];
                        AllPair_fraterho = [AllPair_fraterho; fratecorr];
                        
                        
                    end
                end
                
            end
            
            
            
            
        end
    end
end

% =============== sanity check, plot all corr
if motifcorr_dirty==1
    lt_figure; hold on;
    lt_plot_histogram(RendRhoAll);
end





%% ================= OPUTPUT
OUTSTRUCT.AllPair_pitchrho = AllPair_pitchrho;
OUTSTRUCT.AllPair_fraterho = AllPair_fraterho;




























