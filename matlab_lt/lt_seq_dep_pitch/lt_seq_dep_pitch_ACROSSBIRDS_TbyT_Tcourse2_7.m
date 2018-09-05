function [SlopesAll, Learn_And_FFdevAll] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_7(DATTMP, ...
    use_nminus2, plotON, twind, plotWNhits)

% ==
% twind = [0 2]; % seconds, min and max time devs to consider.
% plotON = 1;
locallearnDirs = {'neg', 'pos'};
minrends = 10; % for regression

%% 
if ~exist('twind', 'var')
    twind = [];
end

    
%%

fnames = fieldnames(DATTMP);

figcount=1;
subplotrows=4;
if plotWNhits==1
    subplotcols = 2;
else
subplotcols=length(fnames);
end
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


SlopesAll = nan(length(DATTMP.Train), 6); % syl x [negTr, negCa, negBase, posTr, posCa, posBase]
% MeanFFdevAll= nan(length(DATTMP.Train), 6); % syl x [negTr, negCa, negBase, posTr, posCa, posBase]
Learn_And_FFdevAll = cell(length(DATTMP.Train), 6);

%%
for j=1:length(DATTMP.Train)
    
    colcount = 1;
    for jj=1:length(locallearnDirs)
        
        learndirthis = locallearnDirs{jj};
        
        % ============ what is median threshold of local learning?
        medthresh =[];
        for fname = fnames'
            fname = fname{1};
            
            try
                tmp = DATTMP.(fname)(j).learnlocal;
                medthresh = [medthresh; nanmedian(tmp)];
            catch err
            end
        end
        %         disp(medthresh);
        medthresh = mean(medthresh);
        
        for ff =1:length(fnames)
            fname = fnames{ff};
            
            % ================
            if plotON==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([fname]);
                hsplots = [hsplots hsplot];
                %                 title([bname '-' ename ',' fname]);
                xlabel('local learn (n minus n-1)');
                ylabel('ff dev (n+1 minus n-1)');
            end
            
            if j<=length(DATTMP.(fname))
                x = DATTMP.(fname)(j).learnlocal;
                if use_nminus2==1
                    y = DATTMP.(fname)(j).ffdev_first + DATTMP.(fname)(j).learnlocal;
                else
                    y = DATTMP.(fname)(j).ffdev_first;
                end
                t = DATTMP.(fname)(j).tdev_first;
                % ---- extract info about WN hits
                wnhits = DATTMP.(fname)(j).wnhits;
                wnmiss = DATTMP.(fname)(j).wnmiss;
                
                % =================== COLLECT DESIRED RENDS
                if isempty(twind)
                    if strcmp(learndirthis, 'neg')
                        indtmp = x<medthresh;
                    elseif strcmp(learndirthis, 'pos')
                        indtmp = x>medthresh;
                    end
                else
                    if strcmp(learndirthis, 'neg')
                        indtmp = x<medthresh & t>=twind(1) & t<=twind(2);
                    elseif strcmp(learndirthis, 'pos')
                        indtmp = x>medthresh & t>=twind(1) & t<=twind(2);
                    end
                end
                x = x(indtmp); 
                y = y(indtmp);
                wnhits = wnhits(indtmp);
                wnmiss = wnmiss(indtmp);
                
                
                % ================== PLOT
                if length(x)>=minrends
                    % ----------------- 1) SLOPE
                    if plotON==1
                        plot(x, y, 'x');
                        lt_plot_makesquare_plot45line(gca, 'b');
                    end
                    [b, bint] = lt_regress(y, x, 0, 0, plotON, plotON, 'r', plotON);                    
                else
                    b = nan(2,1);
                    bint = nan(2,2);
                    x = [];
                    y = [];
                end
                
                % ------- PLOT WN HITS
                if plotON==1
                    if plotWNhits==1
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        %                 title([bname '-' ename ',' fname]);
                        ylabel('hits(r), miss(b)');
                        plot(x, wnhits, 'xr');
                        plot(x, wnmiss, 'xb');
                        
                        % -- plot binned
                        [~, indtmp] = sort(x);
                        wnhits_sm = lt_running_stats(wnhits(indtmp), 10);
                        wnmiss_sm = lt_running_stats(wnmiss(indtmp), 10);
                        x_sm = lt_running_stats(x(indtmp), 10);
                        lt_plot(x_sm.Mean(1:2:end), wnhits_sm.Mean(1:2:end), ...
                            {'Errors', wnhits_sm.SEM(1:2:end), 'Color', 'r'});
                        lt_plot(x_sm.Mean(1:2:end), wnmiss_sm.Mean(1:2:end), ...
                            {'Errors', wnmiss_sm.SEM(1:2:end), 'Color', 'b'});
                    end
                end
            else
                b = nan(2,1);
                bint = nan(2,2);
                x = [];
                y = [];
            end
            
            % =================== SAVE OUTPUT
            SlopesAll(j, colcount) = b(2);
            %             MeanFFdevAll(j, colcount) = [];
            Learn_And_FFdevAll{j, colcount} = [x y];
            
            colcount = colcount+1;
            
        end
        
    end
    assert(colcount ==(length(locallearnDirs)*length(fnames))+1, 'did not iterate over all...')
end
if plotON==1
    linkaxes(hsplots, 'xy');
end