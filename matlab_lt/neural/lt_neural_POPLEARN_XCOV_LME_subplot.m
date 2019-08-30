%% ======= plot
lt_figure; hold on;
lt_subplot(2,2,1); hold on;

% === actual data
% title('actual data');
% plot(dat.exptnum(dat.istarg==1), dat.Yresponse(dat.istarg==1), 'or');
% plot(dat.exptnum(dat.istarg==0 & dat.issame==1), dat.Yresponse(dat.istarg==0 & dat.issame==1), 'ob');
% plot(dat.exptnum(dat.istarg==0 & dat.issame==0), dat.Yresponse(dat.istarg==0 & dat.issame==0), 'ok');

% === actual data
elist = double(dat.exptnum);
clist = double(dat.chanpairID);

title('actual data');
for kk=1:max(elist)
    cc=0;
    for kkk=1:max(clist)

        indsthis = elist==kk & clist==kkk;
        if ~any(indsthis)
            continue
        end
        
        x=double(kk)+cc;
        cc = cc+0.1;
                   
        plot(x, dat.Yresponse(dat.istarg==1 & indsthis), 'or');
        try
        plot(x, dat.Yresponse(dat.istarg==0 & dat.issame==1 & indsthis), 'ob');
        catch err
        end
        try
        plot(x, dat.Yresponse(dat.istarg==0 & dat.issame==0 & indsthis), 'ok');
        catch err
        end
    end
end


% === predicted data
lt_subplot(2,2,2); hold on;
title(['fitted, ' formula]);
F = lme.fitted;
plot(dat.exptnum(dat.istarg==1), F(dat.istarg==1), 'or');
plot(dat.exptnum(dat.istarg==0 & dat.issame==1), F(dat.istarg==0 & dat.issame==1), 'ob');
plot(dat.exptnum(dat.istarg==0 & dat.issame==0), F(dat.istarg==0 & dat.issame==0), 'ok');

% === residuals
lt_subplot(2,2,3); hold on;
title('residuals');
F = lme.residuals;
plot(dat.exptnum(dat.istarg==1), F(dat.istarg==1), 'or');
plot(dat.exptnum(dat.istarg==0 & dat.issame==1), F(dat.istarg==0 & dat.issame==1), 'ob');
plot(dat.exptnum(dat.istarg==0 & dat.issame==0), F(dat.istarg==0 & dat.issame==0), 'ok');

% === residuals
lt_subplot(2,2,4); hold on;
title('fit response, with random noise');
F = lme.random;
plot(dat.exptnum(dat.istarg==1), F(dat.istarg==1), 'or');
plot(dat.exptnum(dat.istarg==0 & dat.issame==1), F(dat.istarg==0 & dat.issame==1), 'ob');
plot(dat.exptnum(dat.istarg==0 & dat.issame==0), F(dat.istarg==0 & dat.issame==0), 'ok');

%% === effect

lt_tools_lme_plotsummary(lme);
if onlyIfSameType==1
    lt_subplot(3,3,9); hold on;
    lt_plot_text(0, 0.5, '1=diff, 2=targ, 3=same');
end