function TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Outli(TrialStruct, ...
    plotRawFF)
%% 6/13/18 - LT, removes outliers from all syllables
% uses something like Tukey criterion, here default is quartile +/- 2*IQR
Numbirds = length(TrialStruct.birds);


%% ===== go thru each syl in each expt, remove data that is outlier
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

AllNumRends = [];
AllNumRemoved =[];
for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:Numexpt
        
        % ---------- SKIP IF NO DATA
        if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
            continue
        end
        
        Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        birdname =  TrialStruct.birds(i).birdname;
        exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        
        
        for ss =1:Numsyls
            
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
            
            T = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
            FF = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
            
            % ====== running median and IQR
            nrends = length(T);
            binsize = 41;
            assert(mod(binsize,2)==1, 'need to be odd?');
            indsoutliers = [];
            upperlims = [];
            lowerlims = [];
            for r =1:nrends
                
                % ====================== 1) CALC LOCAL MEDIAN AND IQR
                indsthis = r-floor(binsize/2):r+floor(binsize/2);
                
                
                % ---- shift if is at edges
                if indsthis(1)<1
                    indsthis = indsthis +(-indsthis(1)+1);
                elseif indsthis(end)>nrends
                    indsthis = indsthis -(indsthis(end)-nrends);
                end
                
                % ----- If N is smaller than binsize...
                if nrends<binsize
                    indsthis = 1:nrends;
                end
                
                ffmed = median(FF(indsthis));
                ffquarts = prctile(FF(indsthis), [25 75]);
                ffiqr = iqr(FF(indsthis));
                
                lim1 = ffquarts(1) - 2*ffiqr;
                lim2 = ffquarts(2) + 2*ffiqr;
                
                % ======================= 2) IS THIS OUTLIER>?
                ffthis = FF(r);
                if ffthis<lim1 | ffthis>lim2
                    % then is outlier
                    indsoutliers = [indsoutliers r];
                end
                
                % =================== COLLECT UPER AND LOWER LIMS FOR PLOT
                upperlims = [upperlims lim2];
                lowerlims = [lowerlims lim1];
                
            end
            
            % ================ PLOT ALL DATA ALONG WITH OUTLIERS
            if plotRawFF==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-' sylname]);
                plot(T, FF, 'ok');
                plot(T, upperlims, '--b');
                plot(T, lowerlims, '--b');
                lt_plot(T(indsoutliers), FF(indsoutliers), {'Color', 'r'});
                %              plot(T(indsoutliers), FF(indsoutliers), 'or');
                axis tight;
            end
            
            % =================== COLLECT ACROSS ALLS YLS
            AllNumRends = [AllNumRends; nrends];
            AllNumRemoved = [AllNumRemoved; length(indsoutliers)];
            
            
            % =================== REMOVE FROM TRIALSTRUCT
            disp(['REMOVED ' num2str(length(indsoutliers)) '/' num2str(nrends) ' outliers for ' birdname '-' exptname '-' sylname]);
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals(indsoutliers) = [];
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals_datenum(indsoutliers) = [];
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals(indsoutliers) = [];
            
            if isfield(TrialStruct.birds(i).exptnum(ii).sylnum(ss), 'isWNhit')
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).isWNhit(indsoutliers) = [];
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).isCatch(indsoutliers) = [];
            end
        end
    end
end

lt_figure; hold on;
xlabel('numrends');
ylabel('num outliers');
title('all syls');
plot(AllNumRends, AllNumRemoved, 'xk');
lt_plot_makesquare_plot45line(gca, 'b', -1);
