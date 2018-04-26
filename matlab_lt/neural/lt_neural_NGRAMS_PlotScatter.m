function lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg)
%% LT 4/24/18 - plots scatter (i.e. mean for pairtype 2 vs. 1)


%% INPUTS

% plottype = 'absfrdiff'; % oneminusrho or absfrdiff
% plotON=0; % only works for absfrdiff
% PairtypesToplot = {...
%     '1  1  1', ... % xaxis
%     '1  0  0'}; % yaxis
%
% % ----------- params for one minus rho, specifically
% dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
% sanitycheckuseneg = 0; % uses negative control data instead of data


%% pull out data
% ========================== FOR COMPATIBILITY WITH OLD CODE, EXTRACT ALL
% FIELDS
fnamesthis = fieldnames(OUTSTRUCT);
for j=1:length(fnamesthis)
    eval([fnamesthis{j} ' = OUTSTRUCT.' fnamesthis{j} ';']);
end

%%
numpairtypes = length(PairtypesToplot);
maxbirds = max(All_birdnum);
maxneur = max(All_neurnum);

pcolors = lt_make_plot_colors(numpairtypes, 0,0);


%%
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% =====================
AllPairs_Means = [];
AllPairs_Bregions = {};
AllPairs_Birdnum = [];
for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
        
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        allmeans = nan(1,numpairtypes);
        
        if strcmp(plottype, 'oneminusrho')
            % ========= plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '[' bregion ']']);
            
            AllRhoNeg = [];
            for k =1:numpairtypes
                
                pairtypethis = PairtypesToplot{k};
                pairtypethis = find(strcmp(PairTypesInOrder, pairtypethis));
                
                inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
                
                % ========= collect data
                pairtype = All_diffsyl_PairType(inds);
                oneminusrho_neg = All_OneMinusRho_NEG(inds);
                oneminusrho = All_OneMinusRho(inds);
                
                if sanitycheckuseneg==1
                    oneminusrho = oneminusrho_neg;
                end
                
                % ----- subtract negative control?
                if dosubtractcontrol==1
                    oneminusrho = oneminusrho - oneminusrho_neg;
                end
                
                % =============== plot histograms
                xcenters = 0.05:0.1:1.95;
                if dosubtractcontrol==1
                    xcenters = -1.15:0.1:1.95;
                end
                lt_plot_histogram(oneminusrho, xcenters, 1, 1, '', 1, pcolors{k});
                
                
                % ================ collect mean to then do scatterplot
                ymean = mean(oneminusrho);
                allmeans(k) = ymean;
                
                % ======== collect negative distribution
                AllRhoNeg = [AllRhoNeg; oneminusrho_neg];
            end
            
            % ======== collect each neuron
            AllPairs_Means = [AllPairs_Means; allmeans];
            AllPairs_Bregions = [AllPairs_Bregions; bregion];
            AllPairs_Birdnum = [AllPairs_Birdnum; i];
            
            % ======== overlay negative distribution
            if dosubtractcontrol==0
                lt_plot_histogram(AllRhoNeg, xcenters, 1, 1, '', 1, [0.7 0.7 0.7]);
            end
            
        elseif strcmp(plottype, 'absfrdiff')
            
            if plotON==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-n' num2str(ii) '[' bregion ']']);
            end
            
            
            for k =1:numpairtypes
                
                pairtypethis = PairtypesToplot{k};
                pairtypethis = find(strcmp(PairTypesInOrder, pairtypethis));
                
                inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
                
                % ========= collect data
                y = All_AbsFRdiff_Zrelshuff(inds);
                
                if plotON==1
                    % =============== plot histograms
                    lt_plot_histogram(y, '', 1, 1, '', 1, pcolors{k});
                end
                % ================ collect mean to then do scatterplot
                ymean = mean(y);
                allmeans(k) = ymean;
            end
            
            % ======== collect each neuron
            AllPairs_Means = [AllPairs_Means; allmeans];
            AllPairs_Bregions = [AllPairs_Bregions; bregion];
            AllPairs_Birdnum = [AllPairs_Birdnum; i];
            
            
        else
            dasfasdfsdf;
        end
        
    end
end

%% =========== PLOT SCATTER COMPARING TWO CLASSES ACROSS ALL NEURONS
if strcmp(plottype, 'oneminusrho')
    lt_figure; hold on;
    % ----- LMAN
    lt_subplot(2,2,1); hold on;
    title('LMAN');
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ok');
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    % ----- RA
    lt_subplot(2,2,2); hold on;
    title('RA');
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ok');
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    % --- COMBINED
    lt_subplot(2,2,3); hold on;
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'og');
    
    % RA
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'or');
    
    % FORMATING
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    
    % ============ SAPRATE BY BIRD
    maxbirds = max(AllPairs_Birdnum);
    for j=1:maxbirds
        lt_figure; hold on;
        
        % LMAN
        indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'og');
        
        % RA
        indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'or');
        
        % FORMATING
        xlabel(PairtypesToplot{1});
        ylabel(PairtypesToplot{2});
        
        xlim([-1 1]);
        ylim([-1 1]);
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        line([-1 1], [-1 1]);
        
        
        
    end
else
    
    maxbirds = max(AllPairs_Birdnum);
    for j=1:maxbirds
        lt_figure; hold on;
        title([SummaryStruct.birds(j).birdname]);
        % LMAN
        indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'og');
        
        % RA
        indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'or');
        
        % FORMATING
        xlabel(PairtypesToplot{1});
        ylabel(PairtypesToplot{2});
        
        %         xlim([-1 1]);
        %         ylim([-1 1]);
        
        lt_plot_makesquare_plot45line(gca, 'k', -2);
    end
end


% =================== COMBINE ALL DATAPOINTS
lt_figure; hold on;
title('all datapoints');
maxbirds = max(AllPairs_Birdnum);
for j=1:maxbirds
    
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'og');
    
    % RA
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'or');
    
    % FORMATING
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    %         xlim([-1 1]);
    %         ylim([-1 1]);
    
    lt_plot_makesquare_plot45line(gca, 'k', -2);
end


%% ========== ancova, 

indstmp = AllPairs_Birdnum==6;
y1 = AllPairs_Means(indstmp, 1);
y2 = AllPairs_Means(indstmp, 2);
group = AllPairs_Bregions(indstmp);

aoctool(y1, y2, group)

%% ====== OVERLAY DISTRIBUTIONS FOR ALL BIRDS
lt_figure; hold on;
title('each bird one fill col (circle(RA), sqiare(L)');
maxbirds = max(AllPairs_Birdnum);
plotcols = lt_make_plot_colors(maxbirds, 0,0);
xlabel(PairtypesToplot{1});
ylabel(PairtypesToplot{2});
for j=1:maxbirds
    
    % ======= fill with color?
    if length(unique(AllPairs_Bregions(AllPairs_Birdnum==j)))>1
        % then color, since multiple brain regions
        pcol = plotcols{j};
    else
        pcol = [1 1 1];
    end
    
    % ============ LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    % get mean and sem in 2 dimensions
    ymean = mean(y,1);
    ysem_x = lt_sem(y(:,1));
    ysem_y = lt_sem(y(:,2));
    
    plot(ymean(1), ymean(2), 's', 'MarkerFaceColor', pcol, ...
        'Color', 'b', 'MarkerSize', 8);
    
    % =============== RA
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    % get mean and sem in 2 dimensions
    ymean = mean(y,1);
    ysem_x = lt_sem(y(:,1));
    ysem_y = lt_sem(y(:,2));
    
    plot(ymean(1), ymean(2), 'o','MarkerFaceColor', pcol, ...
        'Color', 'r', 'MarkerSize', 8);
    
    lt_plot_makesquare_plot45line(gca, 'k', -2);
end


%% ================= FOR EACH NEURON, DIVIDE DATA INTO POS CONTROL

AllPairs_Ratios = AllPairs_Means(:,2)./AllPairs_Means(:,1);

% ====================== PLOT, each bird compare LMAN and RA
maxbirds = max(AllPairs_Birdnum);
for j=1:maxbirds
    lt_figure; hold on;
    title([SummaryStruct.birds(j).birdname]);
    
    Y = {};
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    Y{1} = AllPairs_Ratios(indstmp, :);
    
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    Y{2} = AllPairs_Ratios(indstmp, :);
    
    % === plot
    lt_plot_MultDist(Y, [1 2], 1);
    
    % FORMATING
    ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}])
    xlabel('LMAN -- RA');
    xlim([0 3]);
    
end




%% ====== LINEAR MODEL TESTING WHETHER SLOPE DEPENDS ON BREGION

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1, y2, AllPairs_Birdnum, AllPairs_Bregions, 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
% mdl = 'y1 ~ y2 + bregion:y2 + (-1 + y2|birdnum) + (-1+bregion:y2|birdnum)';  % full model
lme = fitlme(tbl, mdl)


%% =========== SAME, BUT RESTRICTED TO BIRDS WITH BOTH DATA
goodbirds = intersect(AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'LMAN')) ,...
    AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'RA')));
indsgood = ismember(AllPairs_Birdnum, goodbirds');

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
lme = fitlme(tbl, mdl)


%% =========== SAME, BUT RESTRICTED TO BIRDS THAT I RECORDED
% -------- FIND THE GOOD BIRDS
goodbirds = [];
for j=1:length(SummaryStruct.birds)
    if isfield(SummaryStruct.birds(j).neurons(1), 'isRAsobermel')
        continue
    else
        % -- then is my bird
        goodbirds = [goodbirds j];
    end
end
disp({SummaryStruct.birds(goodbirds).birdname});

% ------- INDS FOR THESE BIRDS
indsgood = ismember(AllPairs_Birdnum, goodbirds);

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
lme = fitlme(tbl, mdl)

%% =============== SAME, BUT DO MODEL FOR SINGLE BIRD
birdtodo = 'pu69wh78';
goodbirds = [];
for j=1:length(SummaryStruct.birds)
    if ~strcmp(SummaryStruct.birds(j).birdname, birdtodo)
        continue
    else
        % -- then is my bird
        goodbirds = [goodbirds j];
    end
end
goodbirds = find(strcmp({SummaryStruct.birds.birdname}, birdtodo));

% ------- INDS FOR THESE BIRDS
indsgood = ismember(AllPairs_Birdnum, goodbirds);

% -----
y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1';  % full model
lme = fitlme(tbl, mdl)
