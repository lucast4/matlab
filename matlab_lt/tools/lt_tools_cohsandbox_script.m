
%% ====== show an example
phithis = 0.7; 
nthis = 100
        coh = lt_tools_cohsandbox(phithis, nthis, 1);

%%  how much does coherence vary with changing sample size and phi consistency?

phifraclist = 0:0.1:1;
nlist = [5 10 20 40 75 100 500 1000];

Datout = nan(length(phifraclist), length(nlist));

for j=1:length(phifraclist)
    phithis = phifraclist(j);
    for jj=1:length(nlist)
        nthis = nlist(jj);
        
        coh = lt_tools_cohsandbox(phithis, nthis, 0);
        Datout(j, jj) = coh;
    end
end

figure; hold on;
xlabel('phi consitentcy (1=unform, 0 = one val)');
ylabel('coh');
title('each line diff N');
for j=1:size(Datout,2)
    pcol = [rand rand rand];
    datthis = Datout(:,j);
    x = phifraclist;
    plot(x, datthis, '-o', 'Color', pcol);
    lt_plot_text(x(end), datthis(end), ['n=' num2str(nlist(j))], pcol);
end

%% how much does power correlation influence coherence? (all else fixed)
% THIS INDICATES TAHT IS IT POWER CORRELATION THAT INFLUENCE COH (JUST
% POWER ALONE DOES NOT).
n = 200;
phifrac = 0.25;
mixcoefflist = 0:0.025:1;

cohall = [];
powercorrall =[];
power_a_all = [];
power_b_all = [];

for i=1:length(mixcoefflist)
   mixthis = mixcoefflist(i);
   
   [coh, powercorr, power_a, power_b] = lt_tools_cohsandbox(phifrac, n, 0, mixthis);
   
   cohall = [cohall; coh];
   powercorrall = [powercorrall; powercorr];
   power_a_all = [power_a_all; power_a];
    power_b_all = [power_b_all; power_b];

end
figure; hold on;
subplot(2,2,1); hold on;
xlabel('power corr');
ylabel('coh');
title('power corr gotten by changing mixture of power timecourses');
plot(powercorrall, cohall, '-ok');

subplot(2,2,2); hold on;
xlabel('mix coeff');
ylabel('coh');
plot(mixcoefflist, cohall, 'ok');

subplot(2,2,3); hold on
xlabel('power(a or b)');
ylabel('coh');
plot(power_a_all, cohall, '-ob');
plot(power_b_all, cohall, '-or');

subplot(2,2,4); hold on;
xlabel('mix coeff');
ylabel('power (a or b)');
plot(mixcoefflist, power_a_all, 'ob');
plot(mixcoefflist, power_b_all, 'or');
