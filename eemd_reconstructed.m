
function bpm_track = eemd_reconstructed(x1,bpm_prev)
fs = 125;
allmode = eemd(x1,0.05,5); % 0.05,5;0,1
allmode = allmode(:,2:9);
allmode = allmode'; 

y1=ssa_application(allmode(1,:)',150);
y2=ssa_application(allmode(2,:)',150);
y3=ssa_application(allmode(3,:)',150);
y4=ssa_application(allmode(4,:)',150);
y5=ssa_application(allmode(5,:)',150);
y6=ssa_application(allmode(6,:)',150);
y7=ssa_application(allmode(7,:)',150);
y8=ssa_application(allmode(8,:)',150);
allmode_re = [y1';y2';y3';y4';y5';y6';y7';y8']; 
mean_bpm = [];
for m = 1:8
    number = 0;    
    data_len = length(allmode_re(m,:));
    for n = 1:data_len-1
        if (allmode_re(m,n) <= 0 && allmode_re(m,n+1) >= 0)
            number = number + 1;
        end
        if (allmode_re(m,n) >= 0 && allmode_re(m,n+1) <= 0)
            number = number + 1;
        end
    end
    mean_bpm = [mean_bpm number/2/8*60];
end

diff_mean_bpmprev = abs(mean_bpm - bpm_prev); 
mode_num = find(diff_mean_bpmprev == min(diff_mean_bpmprev));
if (length(mode_num) > 1)
    mode_num = mode_num(find(mode_num>1));
    mode_num = max(mode_num);
end
onemode_sig = allmode_re(mode_num,:);

window = boxcar(1000);
nfft = 4096;
[Pxx_reconstructed f_reconstructed] = periodogram(onemode_sig,window,nfft,fs);
Pxx_reconstructed = Pxx_reconstructed(21:110); 
bpm_reconstructed = f_reconstructed(21:110)*60;

bpm_track = bpm_reconstructed(find(Pxx_reconstructed == max(Pxx_reconstructed)));