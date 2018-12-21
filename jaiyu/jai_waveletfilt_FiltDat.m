function [tmpfiltfh] = jai_waveletfilt_FiltDat(lfp, filter_frequencies)
%% lt 12/13/18 - takes in lfp (one trial), and outputs spectrogram
% filters data using presaved filter kernels


% filter_frequencies=20:2:250;

% tmpfiltfh, cell each corresponding toa  freq, containing hilvert
% transform.
%% load filters
filter_params=filter_load(filter_frequencies);

%% filter data
% save only the hilvert transforms.
[~, ~, ~, ~, tmpfiltfh]=signal_filter_looper(double(lfp),filter_params);


end

function [out] = filter_load(filter_f)

out=cellfun(@(x)load(sprintf('/bluejay5/lucas/analyses/neural/FILTER/JaiYu/cfcampfilt%s.mat',num2str(x))),...
    num2cell(filter_f,1),'uniformoutput',0);
out=cellfun(@(x) x.cfcampfilt,out,'uniformoutput',0);

end


function [outz,out,out_phase,out_filt, out_hilbert]=signal_filter_looper(signal,filters,timefilter)
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
out_hilbert =cell(numel(filters),1);
% parfor gg = 1:numel(filters)
for gg = 1:numel(filters)
    %Compute amplitude time series
    tmpfiltf = filtfilt(filters{gg},1,signal);
    
    % --- lt sanity checke
    if (0)
        figure; plot(signal)
        hold on; plot(filtfilt(filters{12}, 1, signal))
        
        tmp = 20:2:250;
        tmp(12)
        
        datfilt = lt_neural_filter(signal, 1500, 0, 40, 44);
        plot(datfilt, 'm')
    end
    
    tmpfiltfh=hilbert(tmpfiltf);
    tmpfilth = abs(tmpfiltfh);
    tmpfilthphi=angle(tmpfiltfh);
    
    if(0)
        figure; hold on;
        subplot(2,2,1); hold on;
        plot(signal, 'k');
        plot(tmpfiltf, 'b');
        plot(tmpfilth, 'r');
        
        subplot(2,2,2); hold on;
        plot(tmpfilthphi, 'k');
        
        subplot(2,2,3); hold on
        title('orig + Re(hilbert) + imag(hilbert)');
        plot(tmpfiltf, '--k');
        plot(real(tmpfiltfh), ':r');
        plot(imag(tmpfiltfh), '-b');
        plot(abs(tmpfiltfh), '-m');
    end
    
    tmpfilt_z=(tmpfilth-mean(tmpfilth(timefilter)))/std(tmpfilth(timefilter));
    outz{gg}=tmpfilt_z;
    out{gg}=tmpfilth;
    out_phase{gg}=tmpfilthphi;
    out_filt{gg}=tmpfiltf;
    out_hilbert{gg} = tmpfiltfh;
end
end