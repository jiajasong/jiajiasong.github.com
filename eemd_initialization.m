function bpm_track = eemd_initialization(x1)
fs = 125;
allmode = eemd(x1,0.05,5); % ,0,1
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


window = boxcar(1000); 
nfft = 4096;
[Pxx_reconstructed f_reconstructed] = periodogram(y4,window,nfft,fs);
Pxx_reconstructed = Pxx_reconstructed(21:110); 
bpm_reconstructed = f_reconstructed(21:110)*60;
bpm_track = bpm_reconstructed(find(Pxx_reconstructed == max(Pxx_reconstructed)));