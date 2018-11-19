%% play around with phaseoff to see how lag affects phi.

phaseoff = 15;

%%
t = 1:100;
% x = sin(t)+rand(size(t));
% y = sin(t)+rand(size(t));
x = sin(t-phaseoff);
y = sin(t);

[C,phi, ~, ~, ~, f] = coherencyc(x, y);

figure;
subplot(3,2,1); hold on;
title('blue(1), red(2)');
plot(t, x, '-b');
plot(t, y, '-r');

subplot(3,2,2); hold on;
ylabel('coh');
xlabel('f');
plot(f, C);

subplot(3,2,4); hold on;
ylabel('phi');
xlabel('f');
plot(f, phi);