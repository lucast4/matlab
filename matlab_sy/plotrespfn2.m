function plotrespfn2(rsptimes,meanvals,colors,rawsong,numdiv,stimleng,randvecfin,song_fs,stimnames,xlist,linestyles,time_bnds);
%rsptimes points with a resp to be plotted
%meanvals list of meanvals at those resp, size is RxC, where
%R is the number of Responses, and C are stimuli to be plotted, 1st column
%is assumed to be BOS, and values are normalized to that level.

%overlap of .5s
overlap=.5
plots=[1 3]
%calculate song_bnds to plot
totleng=length(rawsong)/song_fs-2
%assume numdiv is 2, anything else is too complicated for me right now.
if(numdiv>1)
    time_bnds=[2 (totleng+overlap)/numdiv+2 stimleng];
    

end
%quick and dirty technique, plot it twice and then cut it off.
%first plot
clear rat_resp
figure
for plti=1:numdiv
    clear h
    pltcount=plots(plti);
    
    %these are the song plots
    subplot(3*numdiv,1,pltcount+2)
    [sm,sp,t,f]=evsmooth(rawsong,song_fs,0.01);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    hold on;
    ind=find(randvecfin==1);
    xvals=xlist(ind);
    xvals2=xlist(ind+1);
    yarray=2500;
    drawxrasters(xvals',xvals2',yarray,'r')
    axis([time_bnds(1+plti-1) time_bnds(2+plti-1) 2500 8000])
    box off;
    set(gca,'ycolor','w');
    
        Xlabel('Time (S)', 'Fontsize', 14);
    
    %colormap('gray');
    subplot(3*numdiv,1,pltcount:pltcount+1)
    bosresp=meanvals(:,1);
    
    if(plti==1)
        for i=2:length(meanvals(1,:))
            rat_resp(:,i)=log2(meanvals(:,i)./bosresp)
        end
    end
     for i=2:length(meanvals(1,:))
            plot(mean(rsptimes,2), rat_resp(:,i),'O','Color',colors(i-1));
            hold on;
            h(i-1)=plot(mean(rsptimes,2), rat_resp(:,i), 'Linewidth',2,'Color',colors(i-1),'Linestyle',linestyles{i-1})
     end
    
    plot(mean(rsptimes,2),zeros(length(rsptimes(:,1))),'r:','Linewidth',2)
     
     
     if(plti==numdiv)
        legend(h,stimnames);
        legend boxoff;
     end
        maxvl=max(max(rat_resp));
       minvl=min(min(rat_resp));
        set(gca,'xcolor','w');
        axis([time_bnds(1+plti-1)/.005-overlap time_bnds(2+plti-1)/.005 minvl maxvl]);
        box off;
        Ylabel('Log2(stim/bos)','Fontsize',14)
end
end


