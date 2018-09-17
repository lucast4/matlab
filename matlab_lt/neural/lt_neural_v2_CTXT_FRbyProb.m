%%

pretime = 0.15; % in sec, time to keep for plotting.
posttime = 0.1;
minN = 7;
bregions = {'LMAN'};

%%

motif_predur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;


%%

apos =1 ;
% i = 1;
% bb = 1;
% nn = 1;


%%

OUTSTRUCT.All_classnum = [];
OUTSTRUCT.All_birdnum = [];
OUTSTRUCT.All_branchnum = [];
OUTSTRUCT.All_neurnum = [];

OUTSTRUCT.All_N = [];
OUTSTRUCT.All_frmat = {};
OUTSTRUCT.All_frx = {};

% ################################################
numbirds = length(ALLBRANCH.alignpos(apos).bird);
for i=1:numbirds
    numbranches = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for bb=1:numbranches
        
        numneurons = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron);
        for nn=1:numneurons
            DAT = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn);
            
            SummaryStruct.birds(i).neurons(nn)
            if isempty(DAT.FR)
                continue
            end
            
            nclasses = length(DAT.FR.classnum);
            
            for cc =1:nclasses
                
                N = size(DAT.FR.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                frmat = DAT.FR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                t = DAT.FR.classnum(cc).FRsmooth_xbin_CommonTrialDur;
                
                if N<minN
                    continue
                end
                
                % ================ clip time base
                [frmat, t] = lt_neural_QUICK_GetFRmat(frmat, t, motif_predur, pretime, posttime);
                
                if isempty(frmat)
                    continue
                end
                   
                % ========================== output
                OUTSTRUCT.All_classnum = [OUTSTRUCT.All_classnum; cc];
                OUTSTRUCT.All_birdnum = [OUTSTRUCT.All_birdnum; i];
                OUTSTRUCT.All_branchnum = [OUTSTRUCT.All_branchnum; bb];
                OUTSTRUCT.All_neurnum = [OUTSTRUCT.All_neurnum; nn];
                
                OUTSTRUCT.All_N = [OUTSTRUCT.All_N; N];
                OUTSTRUCT.All_frmat = [OUTSTRUCT.All_frmat; frmat];
                OUTSTRUCT.All_frx = [OUTSTRUCT.All_frx; t'];
            end
        end
    end
end


%% ======= plot
figcount=1;
% subplotrows=4;
% subplotcols=6;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ===========================
maxbirds = max(OUTSTRUCT.All_birdnum);
maxbranchnum = max(OUTSTRUCT.All_branchnum);
maxneur = max(OUTSTRUCT.All_neurnum);

for i=1:maxbirds
    for ii = 1:maxbranchnum
        for nn = 1:maxneur
            
            indsthis = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_branchnum==ii ...
                & OUTSTRUCT.All_neurnum==nn;
            
            if ~any(indsthis)
                continue
            end
            
            classlist = OUTSTRUCT.All_classnum(indsthis);
            Nlist = OUTSTRUCT.All_N(indsthis);
            frmatlist = OUTSTRUCT.All_frmat(indsthis);
            frxlist = OUTSTRUCT.All_frx(indsthis);
            
            
            % ========================== plot fr, sorted in color by
            % frequency
            [Nlist, indstmp] = sort(Nlist);
            classlist = classlist(indstmp);
            frmatlist = frmatlist(indstmp);
            frxlist = frxlist(indstmp);
            
            plotcols = lt_make_plot_colors(length(Nlist), 1, [1 0 0]);
            
            % ======= 1) PLOT FIRING RATE
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            for j=1:length(classlist)
                
                x = frxlist{j};
                y = mean(frmatlist{j}, 2);
                ysem = lt_sem(frmatlist{j}');
                
                if length(classlist)>3
                    plot(x, y, 'Color', plotcols{j}, 'LineWidth', 2);
                else
                    shadedErrorBar(x, y, ysem, {'Color', plotcols{j}}, 1);
                end
                
                % --- line
                line([0 0], ylim, 'Color', 'b');
                axis tight
                
                YLIM = ylim;
                ylim([0 YLIM(2)]);
                
            end
            
            % ======== 2) PLOT HISTOGRAM OF FREQUENCIES
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            xlabel('classnum');
            ylabel('freqs');
            
            scatter(1:length(classlist), Nlist, [], cell2mat(plotcols'))
            set(gca, 'XTick', 1:length(classlist), 'XTickLabel', classlist);
            
            YLIM = ylim;
            ylim([0 YLIM(2)]);
            
        end
    end
end


