function lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub3(DATSTRUCT)
%% lt 7/27/18 - ALL PLOTS LOOKING AT TARGET SYL ONLY
% Relationship between bias and learning
% Learning and CV change

%%


% ====== FIRST, subsample all to only iclude target syls
inds = find(DATSTRUCT.All_Istarg==1);
DATSTRUCT = lt_structure_subsample_all_fields(DATSTRUCT, inds, 1);


% ===== SECOND, extract all fields to qworkspace variables
fnames = fieldnames(DATSTRUCT)';
for fn=fnames
    fn = fn{1};
    eval([fn ' = DATSTRUCT.(fn);'])
end

%% ###########################################################
%% ########################## BIAS AND LEARNING
%% ============= PLOT SUMMARY
NumBirds = max(All_Birdnum);
lt_figure; hold on;

% --------
lt_subplot(3,2,1); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('FF minus base dur learn (dir of learn)');

bias = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
learn_PBS = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
learn_MUSC = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_PBS);
for j=1:length(bias)
    line([bias(j) bias(j)], [learn_PBS(j) learn_MUSC(j)], 'Color', [0.7 0.7 0.7]);
end
plot(bias, learn_PBS, 'ko');
plot(bias, learn_MUSC, 'ro');


% ---------
lt_subplot(3,2,2); hold on;
xlabel('learn (targ dir)');
ylabel('consolidated learning (hz)');

allbias = All_FF_BASE_PBS - All_FF_BASE_MUSC;
x = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
y = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_PBS);

% -- bias in direciton fo learing
indstokeep = sign(allbias) == sign(All_learndir)

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'bo');


% -- bias in direciton opposite to learing
indstokeep = sign(allbias) ~= sign(All_learndir);

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'mo');

% ---
lt_plot_makesquare_plot45line(gca, 'k');

% ---------
lt_subplot(3,2,3); hold on;
ylabel('fraction consolidation');
xlabel('afp bias (in dir of learnig');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = (All_FF_WN_MUSC - All_FF_BASE_PBS)./(All_FF_WN_PBS - All_FF_BASE_PBS);
lt_regress(y, x, 1);
% ---- connect same bird with lines
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
    
    [~, indsort] = sort(xthis);
    xthis = xthis(indsort);
    ythis = ythis(indsort);
    
    for k=1:length(xthis)-1
        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]);
    end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ------
lt_subplot(3,2,4); hold on;
title('no difference in duration from learn day1');
xlabel('afp bias, dir of learn');
ylabel('days (b=1st; r=last)');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
yblue = All_WNdayrange(:,1);
yred = All_WNdayrange(:,2);

% lt_regress(yblue, x, 1, 0, 1, 1, 'b');
% lt_regress(yred, x, 1, 0, 1, 1, 'r');
plot(x, yred+0.1, 'or');
plot(x, yblue-0.1, 'ob');
for i=1:length(x)
    line([x(i) x(i)], [yblue(i)-0.1 yred(i)+0.1], 'Color', [0.6 0.6 0.6]);
end


% ===================== CV CHANGE FROM BASELINE
lt_subplot(3,2,5); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('CV (bu=base; rd=Train)');
title('PBS data');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y_blue = All_CV_BASE_PBS;
y_mag = All_CV_WN_PBS;

plot(x, y_blue, 'ob');
plot(x, y_mag, 'or');
for i=1:length(x)
    if y_blue(i)>y_mag(i)
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'b');
    else
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'r');
    end
end
lt_plot_zeroline;

% -----------------
lt_subplot(3,2,6); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('CV (bu=base; rd=Train)');
title('MUSC data');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y_blue = All_CV_BASE_MUSC;
y_mag = All_CV_WN_MUSC;

plot(x, y_blue, 'ob');
plot(x, y_mag, 'or');
for i=1:length(x)
    if y_blue(i)>y_mag(i)
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'b');
    else
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'r');
    end
end
lt_plot_zeroline;


%% ========== SAME ANALYSES, BUT MUSC AND PBS SUBTRACTING THEIR RESPECTIVE BASELINES
lt_figure; hold on;

% --------
lt_subplot(3,2,1); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('FF minus base dur learn (dir of learn)');
title('MUSC NORM TO MUSC');

bias = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
learn_PBS = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
learn_MUSC = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_MUSC);
for j=1:length(bias)
    line([bias(j) bias(j)], [learn_PBS(j) learn_MUSC(j)], 'Color', [0.7 0.7 0.7]);
end
plot(bias, learn_PBS, 'ko');
plot(bias, learn_MUSC, 'ro');
lt_regress(learn_MUSC, bias, 1);


% ---------
lt_subplot(3,2,2); hold on;
xlabel('against -- towards');
ylabel('FF minus base dur learn (dir of learn)');
title('MUSC NORM TO MUSC');
x = [1 2];
Y = {};
Y{1} = learn_MUSC(bias<0);
Y{2} = learn_MUSC(bias>0);
lt_plot_MultDist(Y, x, 1, 'k');


% ---------
lt_subplot(3,2,3); hold on;
xlabel('learn (targ dir)');
ylabel('consolidated learning (hz)');
title('MUSC NORM TO MUSC');

allbias = All_FF_BASE_PBS - All_FF_BASE_MUSC;
x = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
y = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_MUSC);

% -- bias in direciton fo learing
indstokeep = sign(allbias) == sign(All_learndir);

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'bo');


% -- bias in direciton opposite to learing
indstokeep = sign(allbias) ~= sign(All_learndir);

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'mo');

% ---
lt_plot_makesquare_plot45line(gca, 'k');

% ---------
lt_subplot(3,2,4); hold on;
ylabel('fraction consolidation');
xlabel('afp bias (in dir of learnig');
title('MUSC NORM TO MUSC');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = (All_FF_WN_MUSC - All_FF_BASE_MUSC)./(All_FF_WN_PBS - All_FF_BASE_PBS);
lt_regress(y, x, 1);
% ---- connect same bird with lines
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
    
    [~, indsort] = sort(xthis);
    xthis = xthis(indsort);
    ythis = ythis(indsort);
    
    for k=1:length(xthis)-1
        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]);
    end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ---------
lt_subplot(3,2,5); hold on;
xlabel('against -- towards');
ylabel('fraction consol');
title('MUSC NORM TO MUSC');
Y = {};
Y{1} = y(x<0);
Y{2} = y(x>0);
lt_plot_MultDist(Y, [1 2], 1, 'k');


%% ========== BASELINE BIAS A FUNCTION OF BASELINE PITCH (PBS)

lt_figure; hold on;

xlabel('baseline pitch (PBS)');
ylabel('baseline pitch (MUSC)');

x = All_FF_BASE_PBS;
y = All_FF_BASE_MUSC;

plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'b');



%% ###########################################################
%% ########################## CV AND LEARNING

lt_figure; hold on;

% ============== 1)
lt_subplot(3,2,1); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('change in CV (difference)');
title('PBS');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;

% ==============
lt_subplot(3,2,2); hold on;
xlabel('AGAINST -- TOWARDS');
ylabel('Change in CV');
title('PBS');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
x = x>0;
y = All_CV_WN_PBS - All_CV_BASE_PBS;

% plot(x, y, 'ok');

[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});

lt_plot([0 1]+0.2, ymean, {'Errors', ysem, 'Color', 'r'});

% --- lines between paired experiments
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
    %
    %     [~, indsort] = sort(xthis);
    %     xthis = xthis(indsort);
    %     ythis = ythis(indsort);
    %
    %     for k=1:length(xthis)-1
    %        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]);
    %     end
    % plot with own jitter and color
    xthis = xthis + 0.3*(rand-0.5);
    pcol = 0.9*[rand rand rand];
    lt_plot(xthis, ythis, {'Color', pcol});
end


xlim([-1 2]);
lt_plot_zeroline;
p = ranksum(y(x==0), y(x==1));
lt_plot_pvalue(p, 'ranksum',1);


% ==============
lt_subplot(3,2,3); hold on;
xlabel('AGAINST -- TOWARDS');
ylabel('Change in CV');
title('MUSC');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
x = x>0;
y = All_CV_WN_MUSC - All_CV_BASE_MUSC;

% plot(x, y, 'ok');

[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});

lt_plot([0 1]+0.2, ymean, {'Errors', ysem, 'Color', 'r'});

% --- lines between paired experiments
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
    %
    %     [~, indsort] = sort(xthis);
    %     xthis = xthis(indsort);
    %     ythis = ythis(indsort);
    %
    %     for k=1:length(xthis)-1
    %        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]);
    %     end
    % plot with own jitter and color
    xthis = xthis + 0.3*(rand-0.5);
    pcol = 0.9*[rand rand rand];
    lt_plot(xthis, ythis, {'Color', pcol});
end


xlim([-1 2]);
lt_plot_zeroline;
p = ranksum(y(x==0), y(x==1));
lt_plot_pvalue(p, 'ranksum',1);




% ============== 1)
lt_subplot(3,2,4); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('change in CV (difference)');
title('MUSC');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = All_CV_WN_MUSC - All_CV_BASE_MUSC;
lt_regress(y, x, 1);
lt_plot_zeroline;




% ============== 1)
lt_subplot(3,2,5); hold on;
xlabel('magnitude of AFP bias (during train)');
ylabel('change in CV (train - base)');
title('PBS');
x = All_learndir.*(All_FF_WN_PBS - All_FF_WN_MUSC);
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;



% ===============
lt_subplot(3,2,6); hold on;
xlabel('magnitude of AFP bias, baseline normalized (during train)');
ylabel('change in CV (train - base)');
title('PBS');
x =  All_learndir.* ((All_FF_WN_PBS - All_FF_BASE_PBS) - (All_FF_WN_MUSC - All_FF_BASE_MUSC));
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;

%%  ############# baseline: magntiude of AFP bias correlate with magnitude of CV change?

lt_figure; hold on;

% ----------- SIGNED VALUE
lt_subplot(2,2,1); hold on;
xlabel('base bias');
ylabel('CV (PBS - MUSC)');
bias = All_FF_BASE_PBS - All_FF_BASE_MUSC;
cv = All_CV_BASE_PBS - All_CV_BASE_MUSC;
plot(bias, cv, 'ok')
lt_plot_zeroline;
lt_plot_zeroline_vert;

lt_subplot(2,2,2); hold on;
xlabel('base bias (abs val)');
ylabel('CV (PBS - MUSC)');
bias = abs(All_FF_BASE_PBS - All_FF_BASE_MUSC);
cv = All_CV_BASE_PBS - All_CV_BASE_MUSC;
plot(bias, cv, 'ok');
lt_regress(cv, bias, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% ############## CHANGE IN PITCH CV, USING RE-WINDOWED PC 
% rewindowing was done for pitch contour analysis. I reanalyze here
% becaosem I am more confident in the re-windowed pitch

nsyls = length(DATSTRUCT.All_PitchCont_BASE_PBS);

Yall = nan(nsyls, 2); % PBS CV [baseline, training]; mean of day CVs
for j=1:nsyls
    
    % --------------- BASELINE
    if (0)
        func_cv = @(x)(std(x)/mean(x));
        ccthis = mean(cellfun(func_cv, DATSTRUCT.All_PitchCont_BASE_PBS(j).All_ffvals)); % mean of day CV
        Yall(j,1) = ccthis;
    else
        twind = DATSTRUCT.All_PitchCont_BASE_PBS(j).All_twind(1,:);
        nday = length(DATSTRUCT.All_PitchCont_BASE_PBS(j).All_PCmat);
        cvall = [];
        for dd=1:nday
            pc_all = DATSTRUCT.All_PitchCont_BASE_PBS(j).All_PCmat{dd};
            ff = mean(pc_all(:, twind(1):twind(2)), 2);
            cvthis = std(ff)/mean(ff);
            cvall = [cvall cvthis];
        end
        cvmean = mean(cvall);
        Yall(j,1) = cvmean;
    end
    % -------------- WN
    func_cv = @(x)(std(x)/mean(x));
    ccthis = mean(cellfun(func_cv, DATSTRUCT.All_PitchCont_WN_PBS(j).All_ffvals)); % mean of day CV
    Yall(j,2) = ccthis;
end

% ================
lt_figure; hold on;

% -------------- 1)
lt_subplot(3,2,1); hold on
title('PBS, change in CV (WN minus base)');
xlabel('old CV (using original extracted pitch');
ylabel('new CV (using new time window for PC');

x = All_CV_WN_PBS - All_CV_BASE_PBS;
y = Yall(:,2) - Yall(:,1);

plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');

%% ############# CHANGE IN WIGGLE AS FUNCTION OF DIRECTION OF LEARNING
% ============================= EXTRACT SUMMARY OF WIGGLES
if isfield(DATSTRUCT, 'Wiggle_WN')
    lt_figure; hold on;
    
    % ============== 1)
    lt_subplot(3,2,1); hold on;
    xlabel('baseline AFP bias (pos = dir of learning)');
    ylabel('change in wiggle (difference)');
    title('PBS');
    x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
    y = Wiggle_WN - Wiggle;
   
    lt_regress(y, x, 1);
    lt_plot_zeroline;
    
    % ==============
    lt_subplot(3,2,2); hold on;
    xlabel('AGAINST -- TOWARDS');
    ylabel('Change in wiggle');
    title('PBS');
    x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
    x = x>0;
    y = Wiggle_WN - Wiggle;
     
    [ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
    
    lt_plot([0 1]+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    
    % --- lines between paired experiments
    for j=1:NumBirds
        indsthis = find(All_Birdnum==j);
    
        xthis = x(indsthis);
        ythis = y(indsthis);
        %
        %     [~, indsort] = sort(xthis);
        %     xthis = xthis(indsort);
        %     ythis = ythis(indsort);
        %
        %     for k=1:length(xthis)-1
        %        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]);
        %     end
        % plot with own jitter and color
        xthis = xthis + 0.3*(rand-0.5);
        pcol = 0.9*[rand rand rand];
        lt_plot(xthis, ythis, {'Color', pcol});
    end
    
    
    xlim([-1 2]);
    lt_plot_zeroline;
    p = ranksum(y(x==0), y(x==1));
    lt_plot_pvalue(p, 'ranksum',1);
    
    
%     % =============
%     lt_subplot(3,2,3); hold on;
%     ylabel('CV_change/Wiggle_change');
%     xlabel('baseline AFP bias (pos = dir of learning)');
%     title('PBS');
%     x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
% 
%     y_wiggle = Wiggle_WN - Wiggle;
%     y_cv = All_CV_WN_PBS - All_CV_BASE_PBS;
%     y = y
%     
%     
%     lt_regress(y, x, 1);
%     lt_plot_zeroline;
    
    
    
else
    lt_figure; hold on;
    lt_plot_text(0, 0.5, 'have not extracted Wiggle_WN, so cant look at wiggle change');
end