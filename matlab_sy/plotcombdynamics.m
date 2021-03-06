%rewritten extensively 4.16.09,
%to incorporate new figure requirements.

%takes sumdyn,
%and a ps (plotstruct).
%ps.minx - minimum x value to include
%ps.maxx - maximum x value to include.
%ps.col,  col={ 'k' 'r','c'}
%ps.colvec
%ps.type, (pct or off or dis)
%ps.addx
%ps.excludebs
%ps.plotavg=1
%ps.comb
%ps.normeff=0
%ps.plot_type, 1 is for shift, 2 is for asymp, 3 is for rev.


function [outvlaczcomb,outvlmuzcomb]=plotcombdynamics(sumdynin,sumbs,ps)
tm_matrix=-4:.2:ps.maxx
axes(ps.ax);
axis(ps.plotbnds);
% This is to separate each individual bsvl, i.e. plotseparate...
if(ps.plotsep)
    bsnumlist=[]
    for ii=1:length(sumdynin)
        bsnumlist=[bsnumlist sumdynin(ii).bsnum]
        
    end
    bsvls=unique(bsnumlist);
  
    for ii=1:length(bsvls)
      plotlist{ii}=find(bsnumlist==bsvls(ii));
  end
else
    plotlist{1}=1:length(sumdynin);
end

 
    
for pnum=1:length(plotlist)
    sumdyn=sumdynin(plotlist{pnum});
%     ax(pnum)=subplot(1, length(plotlist), pnum);
    sdynln=length(sumdyn)
    bsindcomb=[];
    sumdyncomb=[];
    outstruct=[];   
    ln=length(tm_matrix);
    outvlaczcomb=zeros(sdynln,ln);
    outvlmuzcomb=zeros(sdynln,ln);
    ct=1;
    for ii=1:length(sumdyn) 
        smcr=sumdyn(ii);
        crsumbs=sumbs(smcr.bsnum);
        
        if(~isempty(smcr.acz))
                if(ps.usepct)
                    aczvl=smcr.ac_pct;
                    muzvl=smcr.mu_pct;
                else
                    aczvl=smcr.acz
                    muzvl=smcr.muz
                end
        
%    
%             if(smcr.acz(1)>0)
%                 smcr.drxn=1;
%             else
%                 smcr.drxn=0;
            if(ps.flip)
                if(smcr.drxn=='do'&~ps.usepct)
                    aczvl=-aczvl
                    muzvl=-muzvl
                end
            end
        end
if(~isempty(smcr.exadjtms))
 if((max(smcr.exadjtms)>ps.minx)&&(min(smcr.exadjtms)<ps.minx))
    [acvls]=interp1(smcr.exadjtms,aczvl,tm_matrix);
    [muvls]=interp1(smcr.exadjtms,muzvl,tm_matrix);
    
%                 slope=calcslope(smcr.tms(indx),yvl(indx));
%                 equalflag=1;
  
%     if(length(acvls)>ps.maxx)
%         ln=ps.maxx+1;
%     else
%         ln=length(acvls);
%     end
    
        outvlaczcomb(ct,1:ln)=acvls(1:ln)
        outvlmuzcomb(ct,1:ln)=muvls(1:ln)
        ct=ct+1;
   
    

            if(ps.plotraw)
              
%             plot([0:ln-1],outvlaczcomb(ct-1,:),'Color',ps.ac_col,'Linewidth',2)
                tmsind=find(smcr.exadjtms<=ps.maxx);
                [out,sortind]=sort(smcr.exadjtms(tmsind));
                tmsind=tmsind(sortind);
                plot(tm_matrix,acvls,'Color',ps.exacfillcol,'Linewidth',2)
                hold on;
         
%                 plot(tm_matrix,acvls,'o','MarkerSize',5,'MarkerFaceColor',ps.ac_col)
%                 plot(tm_matrix,muvls,'o','MarkerSize',5,'MarkerFaceColor',ps.mu_col)
%             hold on;
                plot(tm_matrix,muvls,'Color',ps.exmufillcol,'Linewidth',2)
%             plot([0:ln-1],outvlmuzcomb(ct-1,:),'Color',ps.mu_col,'Linewidth',2)
             end
%     end
        
      end
end
    end
end
        

 inds=[1:length(outvlaczcomb(:,1))]
% figure
 [outacmn,outacstd]=calcmeanstder2(outvlaczcomb(inds,:));
[outmumn,outmustd]=calcmeanstder2(outvlmuzcomb(inds,:));

% figure
rvacmn=[outacmn]
rvmumn=[outmumn]
rvacer=[outacstd]
rvmuer=[outmustd]

if(isfield(ps,'plotbnds'))
    inds=find(tm_matrix>=ps.analbnds(1)&tm_matrix<=ps.analbnds(2));
else
    inds=1:length(tm_matrix);
end    
    xvec=[tm_matrix(inds);]
    fillx=[xvec xvec(end:-1:1)]

    muvec=[rvmumn(inds)+rvmuer(inds)]
    muvec2=[rvmumn(inds)-rvmuer(inds)]
    filly=[muvec muvec2(end:-1:1)]

    acvec=[rvacmn(inds)+rvacer(inds)]
acvec2=[rvacmn(inds)-rvacer(inds)]
filly2=[acvec acvec2(end:-1:1)]

if(~isfield(ps,'mufillcol'))
    mufillcol=[0.82 0.82 0.82]
    acfillcol=[0.92 0.96 0.98]
else
    mufillcol=ps.mufillcol
    acfillcol=ps.acfillcol
end

% %     mupts=mumean(1:length(indvl))
%     acpts=acmean(1:length(indvl))
%     yvec=[avls.initmean{notevl}*ones(length(mupts),1);mupts(end:-1:1)']
 
% figure
% figure
    hold on;
%     fill(xvec,yvec,acfillcol);
filly([1 end])=0
filly2([1 end])=0
indy=find(~isnan(filly)); 
indy2=find(~isnan(filly2)); 
if(ps.ploter)
    fill(fillx(indy),filly(indy),mufillcol,'edgecolor','w');
    fill(fillx(indy2),filly2(indy2),acfillcol,'edgecolor','w');
end
    plot([xvec], rvacmn(inds),'Color',ps.ac_col,'Linewidth',3)
hold on;
plot([xvec],rvmumn(inds),'Color',ps.mu_col,'Linewidth',3)


function [normout]=normvls(vls,norm)
    if(norm)
        normout=vls./vls(1);
    else
        normout=vls;
    end

function [yvl,yvl2,equalflag]=geteffvl(smcr,indx,ps);
     if(ps.comb&ps.normeff)
          %first normalize each
          normtargeff=norm_onevls(smcr.targeff(indx));
          normctrleff=norm_onevls(smcr.contreff(indx));
          yvl=mean([normtargeff;normctrleff]);
          yvl2=yvl;
          equalflag=1;
     elseif(~ps.comb&~ps.normeff)
            yvl=smcr.targeff(indx);
            yvl2=smcr.contreff(indx);
            equalflag=0;
     elseif(ps.comb&~ps.normeff)
            yvl=mean([smcr.targeff(indx); smcr.contreff(indx)]);
            yvl2=yvl;
            equalflag=1;
     elseif(~ps.comb&ps.normeff)
             yvl=norm_onevls(smcr.targeff(indx));
             yvl2=norm_onevls(smcr.contreff(indx));
             equalflag=0;
     end
                
function [yout]=norm_onevls(yin);
    yout=yin/mean(yin);
    
    
    function [slope]= calcslope(xin,yin);
        
        s=polyfit(xin,yin',1);
        slope=s(1);