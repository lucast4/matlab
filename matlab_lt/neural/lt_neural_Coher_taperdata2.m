function [data, tapers, data_proj] =lt_neural_Coher_taperdata2(data,tapers,nfft,Fs)
%% lt 12/10/18 - used to extract tapered data for troubleshooting.

%%
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if nargin < 4; error('Need all input arguments'); end;
data=change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
if (0) % LT, plots each dat pre and post taper
    figure; for j=1:5, subplot(3,2,j); hold on; plot(data(:,j), '-k'); plot(data_proj(:,j), 'b');end
end
J=fft(data_proj,nfft)/Fs;   % fft of projected data
if (0) % LT, plots fft (amplitude) and pjhaes
    ninds = size(J,1)/2;
    figure; 
    plot(abs(J(2:ninds, :)), '-');
    
    figure; hold on;
    for j=1:size(J,2)
        
        jthis = J(2:ninds,j);
        
        theta = atan2(imag(jthis), real(jthis));
        plot(theta, '-');
    end
    
end
