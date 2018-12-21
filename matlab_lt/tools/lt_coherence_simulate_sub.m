%% generate a summed signal given some components
function y = lt_coherence_simulate_sub(freqlist, amplist, indthis, t)

yall = [];
for j=1:length(freqlist{indthis})
   
    f = freqlist{indthis}(j);
    a = amplist{indthis}(j);
    
    y = a*sin(2*pi*f*t);
    yall = [yall; y];
end
y = sum(yall,1);
