clear all; close all;

%% PARAMS [MODIFIABLE]
binomialp = 0.09; % for simulating origianl dataset (101010111, etc).
contiglength = 12; % cycles in each contig.
ncycle = contiglength*1000; % totla number of cycles is multipe of contig period length.

minregular = 4; % number of regular cycling periods (minimum to call a contig period "reguylar")
maxregular = contiglength;

Nshuffs = 1000; % for getting p-val;/

keepwithinbounds=1; % entire cycling period must be within bounds to be considered regular.

%% PARAMS [autmatic]

regstr = ['1{' num2str(minregular) ',' num2str(maxregular) '}'];

% ==== contig boungs
contigons = 1:contiglength:ncycle;
contigoffs = contiglength:contiglength:ncycle;

%% ==== craete a fake datset and find all regular cycling period

dat = rand([1, ncycle]) > binomialp;

%% === mark each contig as 1 (has reg cycl) or 0.

tmp = [0 abs(diff(dat))];
tmp = num2str(tmp);
tmp(tmp==' ') = [];

[regonsets, regoffsets] = regexp(tmp, regstr, 'start', 'end');
regonsets = regonsets-1;
regoffsets = regoffsets-1;

% === find the contig periods that contain regular cycles
% goodcont = [];
% for i=regonsets
%     goodcont = [goodcont; find(contigons<=i & contigoffs>i)];
% end

if keepwithinbounds ==1 % keeps entire cycling period within bounds.
    contigs_good = []; % 1, 0, for each contig, indicating if found reg cycling in it...
    for i=1:length(contigons)
        
        isgood = any(find(contigons(i)<=regonsets & contigoffs(i)>regonsets ...
            & contigons(i)<=regoffsets & contigoffs(i)>regoffsets)); % make sure both onset and offset is within bounds of contig period.
        
        contigs_good = [contigs_good isgood];
    end
else
    contigs_good = []; % 1, 0, for each contig, indicating if found reg cycling in it...
    for i=1:length(contigons)
        
        
        
        contigs_good = [contigs_good any(find(contigons(i)<=regonsets & contigoffs(i)>regonsets))];
    end
end

%% ===== shuffle many times and collect - each time dtermine the "good" contig periods.

contigs_good_AllShuff = []; % shuff x contigs

for n=1:Nshuffs
    
    datshuff = dat(randperm(length(dat)));
    
    %% ==== get good conrtigs for this shuffled dat.
    tmp = [0 abs(diff(datshuff))];
    tmp = num2str(tmp);
    tmp(tmp==' ') = [];
    
    [regonsets, regoffsets] = regexp(tmp, regstr, 'start', 'end');
    regonsets = regonsets-1;
    regoffsets = regoffsets-1;
    
    
    % === find the contig periods that contain regular cycles
    % goodcont = [];
    % for i=regonsets
    %     goodcont = [goodcont; find(contigons<=i & contigoffs>i)];
    % end
    
    if keepwithinbounds ==1 % keeps entire cycling period within bounds.
        contigs_good_shuff = []; % 1, 0, for each contig, indicating if found reg cycling in it...
        for i=1:length(contigons)
            
            isgood = any(find(contigons(i)<=regonsets & contigoffs(i)>regonsets ...
                & contigons(i)<=regoffsets & contigoffs(i)>regoffsets)); % make sure both onset and offset is within bounds of contig period.
            
            contigs_good_shuff = [contigs_good_shuff isgood];
        end
    else
        contigs_good_shuff = []; % 1, 0, for each contig, indicating if found reg cycling in it...
        for i=1:length(contigons)
            
            
            
            contigs_good_shuff = [contigs_good_shuff any(find(contigons(i)<=regonsets & contigoffs(i)>regonsets))];
        end
    end
    
    
    %     contigs_good_shuff = []; % 1, 0, for each contig, indicating if found reg cycling in it...
    %     for i=1:length(contigons)
    %
    %         contigs_good_shuff = [contigs_good_shuff any(find(contigons(i)<=regonsets & contigoffs(i)>regonsets))];
    %     end
    
    contigs_good_AllShuff = [contigs_good_AllShuff; contigs_good_shuff];
end

%% ##################################### 2 WAYS TO CALCULATE P-VALUE

%% =========== 1 == [KAY VERSION] for each origianl good contig, ask how manuy shuffle cases are also good.
% one p-val for each good contig.

tmp = contigs_good_AllShuff(:, logical(contigs_good));
pvals = mean(tmp,1);
figure; hold on;
xlabel('pval (one for each regular cycling conting period');
title('dsitribution of cycle period"s pvals');

lt_plot_histogram(log10(pvals));
line([log10(0.05) log10(0.05)], ylim);

nsig = sum(pvals<0.05)/length(pvals);
lt_plot_annotation(1, ['frac sig: ' num2str(nsig)], 'r');

%% =========== 2 == [altenrative method] for each shuffle, count numebr of good contigs
% gets one p-val for each epoch.

ngood = sum(contigs_good);
ngood_shuff = sum(contigs_good_AllShuff,2);

figure; hold on;
xlabel('number regular cycles detected');
title('one pval per epoch, based on number of regular cycles detected');
lt_plot_histogram(ngood_shuff);
line([ngood ngood], ylim);

p = (sum(ngood_shuff>ngood)+1)/(length(ngood_shuff)+1);
lt_plot_annotation(1, ['pval = ' num2str(p)]);