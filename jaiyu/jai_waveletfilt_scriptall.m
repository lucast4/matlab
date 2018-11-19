

%% do filter

% define frequency
filter_frequencies=2:2:250;
% load filters
filter_params=filter_load(filter_frequencies);



[filtered_sig,~,~]=signal_filter_looper(curreeg_raw.trace,filter_params,timefilter);







function [out] = filter_load(filter_f)

out=cellfun(@(x)load(sprintf('/home/jai/Src/Filter/cfcampfilt%s.mat',num2str(x))),...
    num2cell(filter_f,1),'uniformoutput',0);
out=cellfun(@(x) x.cfcampfilt,out,'uniformoutput',0);

end


function [outz,out,out_phase,out_filt]=signal_filter_looper(signal,filters,timefilter)
if nargin<3
    timefilter=ones(numel(signal),1)==1;
end
% signal - vector of values
% signal_t - timestamps of signal
% filters - n x 1 cell each with filter function
% triggers - vector of times to align signal
% win - 1 x 2 vector seconds (+ve values) relative to
% triggers to align
out=cell(numel(filters),1);
outz=cell(numel(filters),1);
out_phase=cell(numel(filters),1);
out_filt=cell(numel(filters),1);
parfor gg = 1:numel(filters)
    %Compute amplitude time series
    tmpfiltf = filtfilt(filters{gg},1,signal);
    tmpfiltfh=hilbert(tmpfiltf);
    tmpfilth = abs(tmpfiltfh);
    tmpfilthphi=angle(tmpfiltfh);
    
    
    tmpfilt_z=(tmpfilth-mean(tmpfilth(timefilter)))/std(tmpfilth(timefilter));
    outz{gg}=tmpfilt_z;
    out{gg}=tmpfilth;
    out_phase{gg}=tmpfilthphi;
    out_filt{gg}=tmpfiltf;
    
end
end