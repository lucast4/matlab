function [OUTSTRUCT, birdnum] = lt_neural_POPLEARN_Summary(MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum, BregionWantedList)
%% lt 3/13/18 - population xcov, at each switch
% NOTE: currently only works with one pair in BregionWantedList

onlyPlotIfBothPrePostTrials = 0;
% birdnum = 1; % leave empty if include all birds
% exptnum = 1; % leave empty if include all birds
% BregionWantedList = {{'LMAN', 'RA'}};

%%
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    % === can choose to only plot specific bird
    if ~isempty(birdnum)
        if i ~= birdnum;
            continue
        end
    end
    
    
birdname = SwitchStruct.bird(i).birdname;

numexpts = length(SwitchStruct.bird(i).exptnum);
for ii=1:numexpts
    
    % === can choose to only plot specific expt
    if ~isempty(exptnum)
        if ii ~= exptnum;
            continue
        end
    end
    
    numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
    exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
    
    %% ======== collect xcov traces for PRE and POST rel to each switch.
    for iii=1:numswitches
        
        % ---- for this switch, figure out which populations have data
        % overlapping the onset (i.e. has data both pre and post swictch)
        swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
        
        % ============= Go thru all sets of neurons. ask whether part of this
        % switch.
        for ss = 1:numsets
            songfiles = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles{ss};
            songtimes = datenum(songfiles, 'yymmdd_HHMMSS');
            
            inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
            inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
            
            if isempty(inds_pre) & isempty(inds_post)
                continue
            end
            
            if onlyPlotIfBothPrePostTrials ==1
                if isempty(inds_pre) | isempty(inds_post)
                    continue
                end
            end
            
            disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
            
            
            % ############################################### ANALYSIS/PLOTS
            DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss);
            motiflist = {DAT.motif.regexpstr};
            
            %         assert(length(motiflist)==11, 'asdf');
            % ================= 2) COLLECT XCOV FOR ALL SYLS
            for bbbb = 1:length(BregionWantedList)
                %                     lt_figure; hold on;
                bregionwanted = BregionWantedList{bbbb};
                
                for m=1:length(motiflist)
                    
                    motif = motiflist{m};
                    if isempty(DAT.motif(m))
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).CCallpairs = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).xlags_sec = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).CCallpairs = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).xlags_sec = [];
                        
                        continue
                    end
                    
                    if isempty(DAT.motif(m).XCov_neurpair)
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).CCallpairs = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).xlags_sec = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).CCallpairs = [];
                        OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).xlags_sec = [];
                        
                        
                        continue
                    end
                    
                    % ====== find the neuron pairs that match desired brain
                    % regions - note if need to fliplr the xcov
                    bregionlist = {DAT.motif(m).XCov_neurpair(:).bregtmp};
                    indspairs = [];
                    flippair = [];
                    disp(['--------- DESIRED ' bregionwanted]);
                    for bb = 1:length(bregionlist)
                        if all(strcmp(bregionlist{bb}, bregionwanted))
                            % -- keep and don't flip
                            indspairs = [indspairs bb];
                            flippair = [flippair 0];
                        elseif all(strcmp(fliplr(bregionlist{bb}), bregionwanted))
                            % -- keep, and flip
                            indspairs = [indspairs bb];
                            flippair = [flippair 1];
                        end
                    end
                    flippair = logical(flippair);
                    
                    if isempty(indspairs)
                        % then no neuron pairs with these brain regions
                        continue
                    end
                    
                    %% =================================== COLLECT XCOV TRACES
                    ccRealAll = {DAT.motif(m).XCov_neurpair(indspairs).ccRealAll};
                    ccShiftAll = {DAT.motif(m).XCov_neurpair(indspairs).ccShiftAll};
                    xlags_sec = DAT.motif(m).XCov_neurpair(indspairs(1)).x;
                    
                    % -- flip any if needed
                    ccRealAll(flippair) = cellfun(@fliplr, ccRealAll(flippair), 'UniformOutput', 0);
                    ccShiftAll(flippair) = cellfun(@fliplr, ccShiftAll(flippair), 'UniformOutput', 0);
                    
                    
                    % ------------------ FIND BASE AND TRAINING INDS
                    indsundir = [DAT.motif(m).SegExtr_neurfakeID(1).SegmentsExtract.DirSong]==0;
                    songtimes = [DAT.motif(m).SegExtr_neurfakeID(1).SegmentsExtract(indsundir).song_datenum];
                    
                    preInds = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                    postInds = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                    
                    indtmp = length(postInds);
                    postInds_early = postInds(1:floor(indtmp/2));
                    postInds_late = postInds(floor(indtmp/2)+1:end);
                    
                    assert(length(songtimes) == size(ccRealAll{1},1), 'asasdfasd');
                    
                    % ======================================== PRE
                    indsthis = preInds;
                    CCallthis = [];
                    for j=1:length(ccRealAll)
                        
                        % -- for each pair of neurons
                        ccraw = mean(ccRealAll{j}(indsthis,:), 1);
                        ccshift = mean(ccShiftAll{j}(indsthis,:), 1);
                        ccfinal = ccraw - ccshift;
                        
                        CCallthis = [CCallthis; ccfinal];
                    end
                    
                    OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).CCallpairs = CCallthis';
                    OUTSTRUCT.bird(i).expt(ii).swnum(iii).PRE.motif(m).neurset(ss).xlags_sec = xlags_sec';
                    
                    % ======================================== POST
                    indsthis = postInds;
                    CCallthis = [];
                    for j=1:length(ccRealAll)
                        
                        % -- for each pair of neurons
                        ccraw = mean(ccRealAll{j}(indsthis,:), 1);
                        ccshift = mean(ccShiftAll{j}(indsthis,:), 1);
                        ccfinal = ccraw - ccshift;
                        
                        CCallthis = [CCallthis; ccfinal];
                    end
                    
                    OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).CCallpairs = CCallthis';
                    OUTSTRUCT.bird(i).expt(ii).swnum(iii).POST.motif(m).neurset(ss).xlags_sec = xlags_sec';
                end
            end
        end
    end
    
    
end

end

