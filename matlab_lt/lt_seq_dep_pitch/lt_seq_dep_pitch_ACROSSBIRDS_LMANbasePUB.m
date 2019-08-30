function lt_seq_dep_pitch_ACROSSBIRDS_LMANbasePUB(DATSTRUCT)
%% lt 5/26/19 - publication figures, effect of musc on baseline pitch

bnames1 = DATSTRUCT.singlesyls.Birdnames(DATSTRUCT.pairedsyls.Pairs_OriginalInds(:,1));
bnames2 = DATSTRUCT.singlesyls.Birdnames(DATSTRUCT.pairedsyls.Pairs_OriginalInds(:,2));
assert(all(strcmp(bnames1, bnames2)));
bnameInd = grp2idx(bnames1);

DATSTRUCT.pairedsyls.bnameInd = bnameInd;
%% ##################### 1) scatter plot - check if slopes are different
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ============== SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME');
xlabel('PBS');
ylabel('MUSC');

inds=DATSTRUCT.pairedsyls.IsSameSyl;
X=abs(DATSTRUCT.pairedsyls.separationPBS_paired(inds));
Y=abs(DATSTRUCT.pairedsyls.separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% ============== DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF');
xlabel('PBS');
ylabel('MUSC');

inds=DATSTRUCT.pairedsyls.IsSameSyl==0;
X=abs(DATSTRUCT.pairedsyls.separationPBS_paired(inds));
Y=abs(DATSTRUCT.pairedsyls.separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'r');


%% ================= LINEAR MODEL, EFFECT OF SYL TYPE ON SLOPE?
matchData = 1; % if 1, then only keeps diff types that are separated less 
% than the max of the same times.
threshold = 375;

Xpbs = abs(DATSTRUCT.pairedsyls.separationPBS_paired);
Ymusc = abs(DATSTRUCT.pairedsyls.separationMUSC_paired);
IsSame = DATSTRUCT.pairedsyls.IsSameSyl;
bnameInd = categorical(DATSTRUCT.pairedsyls.bnameInd);

dat = table(Xpbs', Ymusc', IsSame', bnameInd, ...
    'VariableNames', {'FFpbs', 'FFmusc', 'Issame', 'bnameInd'});

if matchData==1
    if isempty(threshold)
        maxff = max(Xpbs(IsSame==1));
    else
        maxff = threshold;
    end
    indsgood = Xpbs<=maxff;
    dat = dat(indsgood,:);
end


% model = 'FFmusc ~ -1 + FFpbs*Issame + (-1 + FFpbs*Issame|bnameInd)';
model = 'FFmusc ~ -1 + FFpbs*Issame';
lme = fitlme(dat, model);

% ##################### PLOT THINGS

lt_figure; hold on;

% ============= 1) plot model fits
lt_subplot(2,2,1); hold on;
xlabel('PBS');
ylabel('MUSC');

% --- 1) plot raw data
plot(dat.FFpbs(dat.Issame==1), dat.FFmusc(dat.Issame==1), 'ob');
plot(dat.FFpbs(dat.Issame==0), dat.FFmusc(dat.Issame==0), 'or');

% --- 2) overlay fits
F = lme.fitted;

x = dat.FFpbs(dat.Issame==1);
y = F(dat.Issame==1);
[~, ind1] = min(x);
[~, ind2] = max(x);
line([x(ind1) x(ind2)], [y(ind1) y(ind2)], 'Color', 'b');


x = dat.FFpbs(dat.Issame==0);
y = F(dat.Issame==0);
[~, ind1] = min(x);
[~, ind2] = max(x);
line([x(ind1) x(ind2)], [y(ind1) y(ind2)], 'Color', 'r');


% ========== 2)
lt_subplot(2,2,2); hold on;
lt_tools_lme_plotEffects(lme)




%% ================= LINEAR MODEL, EFFECT OF SYL TYPE ON SLOPE?
matchData = 0; % if 1, then only keeps diff types that are separated less 
% than the max of the same times.
threshold = 375;

Xpbs = abs(DATSTRUCT.pairedsyls.separationPBS_paired);
Ymusc = abs(DATSTRUCT.pairedsyls.separationMUSC_paired);
IsSame = DATSTRUCT.pairedsyls.IsSameSyl;
bnameInd = categorical(DATSTRUCT.pairedsyls.bnameInd);

dat = table(Xpbs', Ymusc', IsSame', bnameInd, ...
    'VariableNames', {'FFpbs', 'FFmusc', 'Issame', 'bnameInd'});

if matchData==1
    if isempty(threshold)
        maxff = max(Xpbs(IsSame==1));
    else
        maxff = threshold;
    end
    indsgood = Xpbs<=maxff;
    dat = dat(indsgood,:);
end


% model = 'FFmusc ~ -1 + FFpbs*Issame + (-1 + FFpbs*Issame|bnameInd)';
model = 'FFmusc ~ -1 + FFpbs*Issame';
lme = fitlme(dat, model);

% ##################### PLOT THINGS

lt_figure; hold on;

% ============= 1) plot model fits
lt_subplot(2,2,1); hold on;
xlabel('PBS');
ylabel('MUSC');

% --- 1) plot raw data
plot(dat.FFpbs(dat.Issame==1), dat.FFmusc(dat.Issame==1), 'ob');
plot(dat.FFpbs(dat.Issame==0), dat.FFmusc(dat.Issame==0), 'or');

% --- 2) overlay fits
F = lme.fitted;

x = dat.FFpbs(dat.Issame==1);
y = F(dat.Issame==1);
[~, ind1] = min(x);
[~, ind2] = max(x);
line([x(ind1) x(ind2)], [y(ind1) y(ind2)], 'Color', 'b');


x = dat.FFpbs(dat.Issame==0);
y = F(dat.Issame==0);
[~, ind1] = min(x);
[~, ind2] = max(x);
line([x(ind1) x(ind2)], [y(ind1) y(ind2)], 'Color', 'r');


% ========== 2)
lt_subplot(2,2,2); hold on;
lt_tools_lme_plotEffects(lme)



%% ==== directly compare sametype vs. diff type, [constraining to those that start similar]

Y=abs(DATSTRUCT.pairedsyls.separationMUSC_paired) ...
    - abs(DATSTRUCT.pairedsyls.separationPBS_paired);
IsSame = DATSTRUCT.pairedsyls.IsSameSyl;

% --- get bird index
bnameInd = DATSTRUCT.pairedsyls.bnameInd;

% ----------------- PLOT
lt_figure; hold on;

lt_subplot(2,2,1); hold on;
xlabel('IsSame');
ylabel('FF (musc-pbs)');
plotSpread(Y, 'distributionIdx', IsSame);
% --- plot summary bar
[ymean, ysem] = grpstats(Y, IsSame, {'mean', 'sem'});
x = double(unique(IsSame));
lt_plot_bar(x+1, ymean, {'Errors', ysem});


lt_subplot(2,2,2); hold on;
title('color = bird');
xlabel('IsSame');
ylabel('FF (musc-pbs)');

bnames = unique(bnameInd)';
pcols = lt_make_plot_colors(length(bnames), 0, 0);
for j=bnames
    indthis = bnameInd==j;
    
    plotSpread(Y(indthis), 'distributionIdx', IsSame(indthis), 'distributionColors', pcols{j});
    % --- plot summary bar
    [ymean, ysem] = grpstats(Y, IsSame, {'mean', 'sem'});
    x = double(unique(IsSame));
    lt_plot_bar(x+1, ymean, {'Errors', ysem});
end

lt_subplot(2,2, 3:4); hold on;
title('color = bird');
xlabel('IsSame');
ylabel('FF (musc-pbs)');

bnames = unique(bnameInd)';
pcols = lt_make_plot_colors(length(bnames), 0, 0);
for j=bnames
    
    indthis = bnameInd==j;
    
    % --- plot summary bar
    [ymean, ysem] = grpstats(Y(indthis), IsSame(indthis), {'mean', 'sem'});
    x = 2*j + unique(IsSame(indthis))-1;
    lt_plot_bar(x+1, ymean, {'Errors', ysem, 'Color', 'w'});
    
    % =-- plot datapoints
    x = 2*j+IsSame(indthis);
    plotSpread(Y(indthis), 'distributionIdx', x, 'distributionColors', pcols{j});
end

% ========= LME - control for bird
bnameInd = categorical(bnameInd);
Y = Y';
IsSame = IsSame';
dat = table(Y, IsSame, bnameInd);
model = 'Y ~ IsSame + (-1+ IsSame|bnameInd)';
lme = fitlme(dat, model);

lt_tools_lme_plotEffects(lme);
