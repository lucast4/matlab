figure        
clear ind;
clear conscresp
clear epochsum;
%for g100o55
plotnames={'rev' 'p10' 'm10'  'bos' 'jig' 'jigc' }
%base='site6sum'
stimorder=[1 2 3 4 5 6 ]
%for g100o55e2
%plotnames={'rev' 'p10' 'm10' 'bos' 'm2st' 'jig' 'jigc' 'jig2' 'jig2c'}
%stimorder=[2 5 3 1 4 6 7 8 9]
%find areas where there is a response%
clusts=[1 2 3];
meanspikerat2=[];
numclust=length(clusts);

for cc=1:length(clusts)
    cvl=clusts(cc);
    
        %calculate the baseline firing rate
        % subplot(rast_dens+3,1,plotcount+2)
        for ii=1:length(stimorder)
            edges=[binsize 2*(1/binsize)];
            numtrials=length(stimf(cvl,1).histdist(:,1))
            histout=sum(stimf(cvl,ii).histdist(:,edges(1):edges(2)),2);
            histout=(histout/length(edges(1):edges(2)));
            meanspikerat(ii)=mean(histout)
            errspikerat(ii)=std(histout)
            %errorbar(meanspikerat ,errspikerat,'+','Linewidth',3)
            %set(gca,'xticklabel',structnames','Fontsize',14);
            %set(get(gca,'Ylabel'),'String','Song Spike
            %Rate','Fontsize',14); 
       end
       
        thresh(cvl)=mean(meanspikerat)+2.5*mean(errspikerat)
        %[usr_in]=input('Threshold okay?')
        %if(isempty(ans))
        %else
        %[x]=input();
        %thresh(cvl)=x;
        %end       
end  %now search all the stimuli to find bins during the song where spike rate is 
  %above this threshold.
  
  for cc=1:numclust
      concatarr=[];
      clst=clusts(cc);
    for ii=1:length(stimorder)
      ind{ii}=find(stimf(clst,ii).meanhist>thresh(clst))
      concatarr=[concatarr ind{ii}];  
    end
  
        concatarr=unique(concatarr);
        respind=find(concatarr>2*(1/binsize));
        n=concatarr(respind);
  %further restrict n so that there must be two consecutive bins.
  %Write a function to find consecutive strings of integers
        out_test=getconsec(n);
        diff=out_test(:,2)-out_test(:,1);
        ind2=find(diff>0);
        conscresp{clst}=out_test(ind2,:);
  
  
  end
  %now in those bins, calculate the sum and standard error
  %each row of meanspikerat, is a stimulus, each column is particular
  %stimulus segment
  for cc=1:numclust
      clst=clusts(cc);  
      for ii=1:length(stimorder)
        stimf(clst,ii).meanspikerat=[];
        stimf(clst,ii).histout=[];
        stimf(clst,ii).errspikerat=[];
     
          %write a loop because anything else is too confusing right now
      for jj=1:length(conscresp{clst}(:,1))
          histout=sum(stimf(clst,ii).histdist(:,conscresp{clst}(jj,1):conscresp{clst}(jj,2)),2);
          numtrials=length(stimf(clst,ii).histdist(:,1));
          stimf(clst,ii).histout{jj}=histout/(conscresp{clst}(jj,2)-conscresp{clst}(jj,1))/binsize;
          stimf(clst,ii).meanspikerat(jj)=mean(stimf(clst,ii).histout{jj});
          stimf(clst,ii).errspikerat(jj)=std(stimf(clst,ii).histout{jj})/sqrt(numtrials);
      end
  end

  end  
for cc=1:numclust  
    clst=clusts(cc);
    cmb_dt=[];
    for jj=1:length(conscresp{clst}(:,1))
        meanspikerat2=[];
        errspikerat2=[];
      
      for ii=1:length(stimorder)
        meanspikerat2=[meanspikerat2 stimf(clst,stimorder(ii)).meanspikerat(jj)];
        
      end
      cmb_dt=[cmb_dt;meanspikerat2];
  end
  epochsum(clst).means=cmb_dt;
  epochsum(clst).rsplst=conscresp{clst}
  epochsum(clst).names=plotnames;
end
%save  
 % eval(['save ',base,'.mat epochsum']);
 save sum.mat epochsum