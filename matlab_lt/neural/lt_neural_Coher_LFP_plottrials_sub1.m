% --- 1) LFP (overlay the channels)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['[' plotname ']' motiftoplot], 'Color', titlecol);
ylabel(['LFP, b-r' bregionpair]);

plot(t_lfp, lfpmat1(indtrial,:), '-b');
plot(t_lfp, lfpmat2(indtrial,:), '-r');

axis tight;
line([0 0], ylim);

% --- 2) LFP (filtered)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['LFP(' num2str(lfpband_lo) '-' num2str(lfpband_hi) '), b-r' bregionpair]);

plot(t_lfp, lt_neural_filter(double(lfpmat1(indtrial,:)), fs, 0, lfpband_lo, lfpband_hi), '-b');
plot(t_lfp, lt_neural_filter(double(lfpmat2(indtrial,:)), fs, 0, lfpband_lo, lfpband_hi), '-r');
axis tight;
ylim([-100 100]);
line([0 0], ylim);

if plotonlyLFP==0
    % --- 2) COHERENCE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('coherence');
    %                 clim = [-0.5 0.5];
    %                 ylabel(['LFP(10-40hz), b-r' bregionpair]);
    imagesc(t_coh, ff_coh, cohmat(:,:, indtrial));
    colorbar
    axis tight;
    ylim([0 ffhilim]);
    
    % --- 3/4) SPECTROGRAM (2 CHANS)
    SpectrumAll = cell(1,2);
    % --- chan1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('chan1')
    [S,t,f] = mtspecgramc(lfpmat1(indtrial,:), movingwin, params);
    t = t-predur_lfp;
    imagesc(t, f, 10*log10(S)');
    colorbar
    axis tight;
    ylim([0 ffhilim]);
    % GET SPECTRUM
    spectrum = mean(S(t>twind_spectrum(1) & t<twind_spectrum(2), :),1);
    SpectrumAll{1} = spectrum;
    
    % --- chan2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('chan2');
    [S,t,f] = mtspecgramc(lfpmat2(indtrial,:), movingwin, params);
    t = t-predur_lfp;
    imagesc(t, f, 10*log10(S)');
    colorbar
    axis tight;
    ylim([0 ffhilim]);
    % GET SPECTRUM
    spectrum = mean(S(t>twind_spectrum(1) & t<twind_spectrum(2), :),1);
    SpectrumAll{2} = spectrum;
    
    
    % ----- get power spectrums at specific time window
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('spectra');
    plot(f, 10*log10(SpectrumAll{1}), '-b');
    plot(f, 10*log10(SpectrumAll{2}), '-r');
    axis tight;
end