
function AEE_N = MovingAverage(BPM_Ori,BPM0,window)
BPM0 = BPM0';
N = length(BPM_Ori);
AEE = mean(abs(BPM_Ori - BPM0));
for i = 1:N - window + 1
    BPM_filt(i) = mean(BPM_Ori(i:i+window-1));
end
for j = N - window + 2:N
    BPM_filt(j) = BPM_Ori(j);
end
AEE_N = mean(abs(BPM_filt - BPM0));