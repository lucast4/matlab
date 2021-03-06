function SegExtract =lt_neural_QUICK_GetRawNeural(SegExtract, SummaryStruct, birdnum, neurnum, ...
    extrapad, motifpredur, flipdat)
%% extrapad adds time before onset and offset (equal)
if ~exist('extrapad', 'var')
    extrapad = 0;
end
if ~exist('flipdat', 'var')
    flipdat=0;
end
%% lt 1/24/18 - extracts raw neural (unfiltered) aligned to data for each trial in segextract

%%



%% go and load neural dat for this seg

[SongDat, ~, ~] = lt_neural_ExtractDat2(SummaryStruct, birdnum, neurnum);
% cd(SummaryStruct.birds(birdnum).neurons(neurnum).dirname);
% tmp = load('data.mat');

Datthis = load([SummaryStruct.birds(birdnum).neurons(neurnum).dirname '/data.mat']);

%% extract neural dat for period in all trials of SegExtract
fs = SegExtract(1).fs;
% motifpredur = SegExtract(1).

for i=1:length(SegExtract)
    disp(['trial ' num2str(i)]);
    ons = SegExtract(i).global_ontime_motifInclFlank;
    off = SegExtract(i).global_offtime_motifInclFlank;
    %    fs = SegExtract(i).fs;
    
    % ---- ADD EXTRPAD
    ons = ons-extrapad;
    off = off+extrapad;
    
    % --- expected duration [based on first trial]
    if i==1
        nsampmax = ceil((off-ons)./(1/SegExtract(i).fs));
    end
    
    % -- sanity check
    indtmp = SegExtract(i).global_tokenind_DatAlignedToOnsetOfThis;
    assert(abs(ons - (SongDat.AllOnsets(indtmp)-motifpredur-extrapad))<0.01, 'probably problem...');
    
    % -- find corresponding samples
    ons_samp = floor(fs*ons);
    off_samp = floor(fs*off);
    
    % ---- make same duerations
    nsamp = off_samp-ons_samp+1;
    if nsamp > nsampmax
        assert(nsamp-nsampmax<10);
        off_samp = off_samp-(nsamp - nsampmax);
    elseif nsamp <nsampmax
        assert(nsampmax-nsamp<10)
        off_samp = off_samp + nsampmax-nsamp;
    end
        
    datseg = single(Datthis.data(ons_samp:off_samp));
    %    datseg = single(datseg);
    
    % get timing onset and offset
    tOnOff_reltoken = [-(motifpredur+extrapad) (off-ons)-(motifpredur+extrapad)];    
    
    if flipdat==1
        datseg = datseg';
        tOnOff_reltoken = tOnOff_reltoken';
    end
    
    % ========== output
    SegExtract(i).neural_rawdat = datseg;
    SegExtract(i).neural_rawdat_tOnOff_reltok = tOnOff_reltoken;
end

