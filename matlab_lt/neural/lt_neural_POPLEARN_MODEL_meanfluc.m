%% lt 2/14/19 - for spiking xcov, should mean flucturation from trial to trial
% lead to localized peaks (in lag space)?

ntrials = 1000;
t = 1:100;
lag = 2;

% ---- add certain types of noise.
% ==== 1) trial by trial mean fluctuation.
noise_tbyt = 0; % std of standard normal.

% ==== 2) shared jitter (i.e. simulates syl segementation errors);
noise_jitter = 0;

% === decoup[le x and y
nocorr = 0;

% ==== make y flat? (instead of sinusoidal)
makeYflat = 0;

% ---- shared noise, 2 versions
sharednoisever =1; % tyhen jitters t that goies into sin... [DEFAULT]
sharednoisever =2; % tyhen adds scalar to random positions, then lags it to creat the other signal.
% also will make y just flat except bumps where x happens to get 1s.

% NOTE:
% TO test whether t by t fluictioantion of mean fralone can lead to peaks, use:
% noise_tbyt = 5; % std of standard normal.
% noise_jitter = 0;
% nocorr = 1;

% TO TEST whether t by t affects structure of peak (assuming one existsa nd
% is real), use:
% noise_tbyt = 1; % std of standard normal.
% noise_jitter = 0;
% nocorr = 0;

% TO TEST whether shared jitter (i.e. segmentation issues) can affect structure of peak (assuming one existsa nd
% is real), use:
% noise_tbyt = 0; % std of standard normal.
% noise_jitter = 1;
% nocorr = 0;

% TO TEST whether jitter alone (assuming no corr) can lead to peak, use
% follwionog:
% noise_tbyt = 0; % std of standard normal.
% noise_jitter = 1;
% nocorr = 1;


% =
% noise_tbyt = 0; % std of standard normal.
% noise_jitter = 0;
% nocorr = 0;
% makeYflat = 1;

%% --- generate signals across trials
t_input = repmat(t, ntrials, 1);

% --- add some noise to t (i.e. shared noise)
if sharednoisever==1
t_input = t_input + 0.5*randn(size(t_input));
end

% ============= ADD ACROSS TRIAL JITTER (I.E. SEGMENTATION NOISE)
ep = noise_jitter*randn(size(t_input,1),1);
ep = repmat(ep, 1, size(t_input,2));
t_input = t_input+ep;

x = sin(t_input);

% ========== add discrete shared noise
if sharednoisever==2
    tmp = 5*median(abs(x(:)))*(rand(size(x))>0.5);
    x = x+tmp;
end

% --- generate y as lagged version of x.
if nocorr==1
    t_input2 = repmat(t, ntrials, 1);
    t_input2 = t_input2 + ep;
    y = sin(t_input2);
else
    if sharednoisever==2
        y = ones(size(x));
        y = y+tmp;
        y = y(:, lag:end);
    elseif sharednoisever==1
        y = x(:, lag:end);
    end
        x = x(:,1:end-lag+1);
    t = t(lag:end);
end
% y = sin(t_input+5);


if makeYflat==1
   y = ones(size(x));
end

% ========== ADD NOISE
% ------- 1) trial by trial mean fluctuation
ep = noise_tbyt*randn(size(y,1),1);
ep = repmat(ep, 1, size(y,2));
y = y+ep;
x = x+ep;
% y = y.*ep;
% x = x.*ep;


%% --- plot individual trial.
if (0)
    lt_figure; hold on;
    plot(t(1,:), x(1,:), 'k');
    plot(t(1,:), y(1,:), 'r');
    
    [cc, lags] = xcorr(x(1,:), y(1,:), 20, 'unbiased');
    lt_figure; hold on;
    plot(lags, cc, '-k');
end

%% ================ get xcorr

% === go thru all trials
ntrials = size(x,1);

CCall = [];
CCall_shift = [];
for i=1:ntrials
    
    % --data
    xthis = x(i,:);
    ythis = y(i,:);
    
    [cc, lags] = xcorr(xthis, ythis, 20, 'unbiased');
    CCall = [CCall; cc];
    
    % =================== AUTOCOVARIANCE
    
    
    % --- shift predictor.
    if i<ntrials
        ythis = y(i+1,:);
    else
        ythis = y(i-1,:);
    end
    [cc, lags] = xcorr(xthis, ythis, 20, 'unbiased');
    CCall_shift = [CCall_shift; cc];
end


% plot
lt_figure; hold on;

% ---
subplot(3,2,1); hold on;
title('shuff(k), dat(r)');
plot(lags, mean(CCall,1), 'r');
plot(lags, mean(CCall_shift,1), 'k');

subplot(3,2,2); hold on;
title('dat minus shuff');
plot(lags, mean(CCall,1)-mean(CCall_shift,1), 'b');


% --- convert to zcorte
tmpmean = mean(CCall_shift,1);
tmpstd =  std(CCall_shift,1);
CCall_z = (CCall - tmpmean)./tmpstd;
subplot(3,2,3); hold on;
plot(lags, mean(CCall_z), 'b');
title('z vs shuff');

%% ===== frequency domain
lt_switch_chronux(1);

% == take fourier transform of xcov function
ccdat = mean(CCall,1);

tmp = fft(ccdat);
tmp2 = tmp./(tmp.*tmp);

ccdat_post = ifft(tmp2);
figure; hold on; plot(lags, ccdat_post)