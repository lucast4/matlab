%% lt 5/24/18 - does direction of AFP bias change during learning?
% ACROSS EXPERIMENTS FOR THE SAME BIRD/SYLLABLE?
%% LT 6/20/16
function lt_seq_dep_pitch_ACROSSBIRDS_Hamish2(SeqDepPitch_AcrossBirds, Params, plotExptRawDat)

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%%
useLearningRelLastBlineDay=0; % otherwise will use rel entire baseline
useWNday2=0; % otherwise will use 1. using 2 allows to get some expt without songs on day 1
takeDay2forExptWithNoDay1=1; % otherwise will throw out those experiments.


%% ========================= PLOT EACH EXPERIMENTS
count = 0;
for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    %     figcount=1;
    %     subplotrows=4;
    %     subplotcols=2;
    %     fignums_alreadyused=[];
    %     hfigs=[];

    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;

        % ====== PLOT RAW DAT FOR THIS DAY TO COMPARE TO EXTRACTED STATS
        plotLarge=1;
        BirdToPlot=birdname;
        ExptToPlot=exptname;
        SylsToPlot={targsyl};
        overlayMeans=1;
        UseSylColors=0; % 0 is black;
        flipsign=1; % plot pos, if neg
        use_std=0; % applies to mean plots only!! (std, instead of sem)
        plotRawFF=1; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
        OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
        OverlayMUSC_days=[];
        plotRunningCV = 0;
        lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds,...
            Params, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, ...
            plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, ...
            OverlayMUSC_days, plotLarge, plotRunningCV)
       
        count = count+1;
    end
end
disp(['Num expts plotted: ' num2str(count)]);


