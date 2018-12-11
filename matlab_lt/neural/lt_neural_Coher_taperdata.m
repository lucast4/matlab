function [data1OUT, data2OUT, data_proj1, data_proj2, tapersOUT, J12, f] = ...
    lt_neural_Coher_taperdata(data1,data2,params)
%% lt modified 12/10/18 to extract tapered data that would be used for coherence.

%%
% Multi-taper coherency,cross-spectrum and individual spectra - continuous process
%
% Usage:
% [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
% Input:
% Note units have to be consistent. See chronux.m for more information.
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms:
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1).
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...).
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional.
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       C (magnitude of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi (phase of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S12 (cross spectrum -  frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S1 (spectrum 1 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S2 (spectrum 2 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - theoretical/jackknife (depending on err(1)=1/err(1)=2) standard deviation for phi.
%                Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
%                bands for phi - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 2; error('Need data1 and data2'); end;
data1=change_row_to_column(data1);
data2=change_row_to_column(data2);
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave]=getparams(params);
if nargout > 8 && err(1)~=2;
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
% if nargout > 6 && err(1)==0;
%     %   Errors computed only if err(1) is nonzero. Need to change params and run again.
%     error('When errors are desired, err(1) has to be non-zero.');
% end;
N=check_consistency(data1,data2);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass);
tapers=dpsschk(tapers,N,Fs); % check tapers

[data1OUT, tapersOUT, data_proj1] =lt_neural_Coher_taperdata2(data1,tapers,nfft,Fs);
[data2OUT, tapersOUT, data_proj2] =lt_neural_Coher_taperdata2(data2,tapers,nfft,Fs);

    J1=mtfftc(data1,tapers,nfft,Fs);
    J2=mtfftc(data2,tapers,nfft,Fs);
    J1=J1(findx,:,:); J2=J2(findx,:,:);
    J12 = conj(J1).*J2;
if (0)
    S12=squeeze(mean(J12, 2));
    S1=squeeze(mean(conj(J1).*J1,2));
    S2=squeeze(mean(conj(J2).*J2,2));
    if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
    C12=S12./sqrt(S1.*S2);
    C=abs(C12);
    phi=angle(C12);
    if nargout>=9;
        [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave);
    elseif nargout==8;
        [confC,phistd]=coherr(C,J1,J2,err,trialave);
        % else
        %     confC = [];
        %     phistd = [];
        %     Cerr = [];
    end
end

% ====== plot amplitude and phase of all fft output for each taper
if (0) % LT, plots fft (amplitude) and pjhaes
    figure;
    hsplots = [];
    
    hsplot = subplot(4,2,1); hold on;
    hsplots = [hsplots hsplot];
    title('fft ampl (1)');
    plot(f, abs(J1), '-');
    
    hsplot = subplot(4,2,2); hold on;
    hsplots = [hsplots hsplot];
    title('fft ampl (2)');
    plot(f, abs(J2), '-');
    
    hsplot = subplot(4,2,3); hold on;
    hsplots = [hsplots hsplot];
    title('angle of fft (1)');
    angleAll = [];
    for j=1:size(J1,2)
        jthis = J1(:,j);
        
        theta = atan2(imag(jthis), real(jthis));
        angleAll = [angleAll theta];
    end
    plot(f, angleAll, '-');
    
    hsplot = subplot(4,2,4); hold on;
    hsplots = [hsplots hsplot];
    title('angle of fft (2)');
    angleAll = [];
    for j=1:size(J2,2)
        jthis = J2(:,j);
        
        theta = atan2(imag(jthis), real(jthis));
        angleAll = [angleAll theta];
    end
    plot(f, angleAll, '-');
    
    hsplot = subplot(4,2,6); hold on;
    hsplots = [hsplots hsplot];
    title('S12(k), S1(b), S2(r)');
    plot(f, abs(S12), 'k');
    plot(f, abs(S1), 'b');
    plot(f, abs(S2), 'r')'
    linkaxes(hsplots, 'x');
    
    
    hsplot = subplot(4,2,5);
    hsplots = [hsplots hsplot];
    title('coherence');
    plot(f, C);
    
    linkaxes(hsplots, 'x');
    
    
    % ==== plot rose plot for a given time bin across tapers
    figure; hold on;
    indf = 4;
    j1this = J1(indf,:);
    j2this = J2(indf,:);
    
    for k=1:length(j1this)
        subplot(4,2,k); hold on;
        title('one time slice, J at multiple tapers');
        plot(real(j1this(k)), imag(j1this(k)), 'o');
        plot(real(j2this(k)), imag(j2this(k)), 'o');
        
        xlim([-10 10]);
        ylim([-10 10]);
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
    end
    
    % ==== plot magnitude and angle correlations across tapers
    figure;
    flist = [1 4 5 6 7];
    for kkk=1:length(flist)
        fff = flist(kkk);
        j1this = J1(fff,:);
        j2this = J2(fff,:);
        
        % -- magnitude
        subplot(5,2,kkk*2-1); hold on;
        title(['ampl [f=' num2str(f(fff))]);
        plot(abs(j1this), 'ok');
        plot(abs(j2this), 'ok');
        lt_plot_zeroline;
        ylim([0 10])
        
        % -- angle
        subplot(5,2,kkk*2); hold on;
        title('angle');
        plot(angle(j1this), 'ok');
        plot(angle(j2this), 'ok');
        line(xlim, [0 0]);
    end
    
    
end

