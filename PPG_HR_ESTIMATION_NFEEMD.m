%notes: the results will be different due to added random nosie and we
%choose the best results.
clc;
clear all;
load('DATA_08_TYPE02'); % For the calculation of the first 12 datasets
load('DATA_08_TYPE02_BPMtrace');
tic
fs = 125;
MADET = [];                              
BPM_Ori = [];                            
prev_count_1 = 0;                       
i_set_1 = [];                            
prev_count_2 = 0;                        
i_set_2 = [];                            
HR_Cal_Start = 0;
test = 1;
[row,column] = size(sig);
num_window = (column - 1000)/250 + 1;
for i = 1:num_window; 
    Pxx1_peak = [];
    Pxx_peak = [];
    bpm1_peak = [];
    bpm_peak = [];
    Pxx_peak_acc = [];
    Pxx_peak_acc_2 = [];
    Pxx_peak_acc_3 = [];
    Pxx_peak_acc_4 = [];
    bpm_peak_acc = [];
    bpm_peak_acc_2 = [];
    bpm_peak_acc_3 = [];
    bpm_peak_acc_4 = [];
    Pxx_peak_tf2 = [];
    bpm_peak_tf2 = [];
    p_loc_remove = []; 
    start = (i*2-1-1)*125+1; % start sampling point of signal
    stop = start+999; % stop sampling point of signal
    % notes: single-channel PPG and tri-axis acceleration data
    x1 = sig(2,start:stop); % PPG signal
    x2 = sig(4,start:stop); % ACCX signal
    x3 = sig(5,start:stop); % ACCY signal
    x4 = sig(6,start:stop); % ACCZ signal
    window = boxcar(1000);
    nfft = 4096; 
    
    lf = 21;
    hf = 110; % 0.6Hz-3.5Hz
    
    % PPG periodogram
    [Pxx1 f1] = periodogram(x1,window,nfft,fs);
    Pxx1 = Pxx1(lf:hf); 
    bpm1 = f1(lf:hf)*60;
    
    % spectrum peaks
    dPxx1 = diff(Pxx1);
    for r = 1:length(dPxx1)-1
        if (dPxx1(r)>0 && dPxx1(r+1)<0 && Pxx1(r+1)>max(Pxx1)*0.3)
            Pxx1_peak = [Pxx1_peak Pxx1(r+1)];
            bpm1_peak = [bpm1_peak bpm1(r+1)];
        end
    end
    
    % initialization
    if (test == 1)
        if (length(bpm1_peak) == 1 && bpm1_peak > 60)
            HR_Cal_Start = 1;
            Index_Start = i;
            test = 0;
        else
            BPM_Ori = [BPM_Ori 0];
        end
    else
        test = 0;
    end
    
    if (HR_Cal_Start == 1)
       %% acc periodogram
        [Pxx2 f2] = periodogram(x2,window,nfft,fs);
        Pxx2 = Pxx2(lf:hf); 
        bpm2 = f2(lf:hf)*60;
        [Pxx3 f3] = periodogram(x3,window,nfft,fs);
        Pxx3 = Pxx3(lf:hf); 
        bpm3 = f3(lf:hf)*60;
        [Pxx4 f4] = periodogram(x4,window,nfft,fs);
        Pxx4 = Pxx4(lf:hf); 
        bpm4 = f4(lf:hf)*60;

       %% pre-processing
        x1_rbw = x1 - Lowpassfilter(x1,fs*2); % removing base wandering
        x1_rbw_rhf = Lowpassfilter(x1_rbw,fs/3);   % removing high-frequency coponents
        [Pxx f] = periodogram(x1_rbw_rhf,window,nfft,fs);
        Pxx = Pxx(lf:hf); 
        bpm = f(lf:hf)*60;

       %% judge the presence of MAs
        amp(1) = max(Pxx2);
        amp(2) = max(Pxx3);
        amp(3) = max(Pxx4);
        if (amp(1)<0.1 && amp(2)<0.1 && amp(3)<0.1)
            ma_detect = 0;
            MADET = [MADET ma_detect];
        else
            ma_detect = 1;
            MADET = [MADET ma_detect];
        end

       %% the previous HR value
        if (i == Index_Start)
            bpm_prev = bpm1(find(Pxx1 == max(Pxx1))); 
            if (bpm_prev < 60)
                bpm_prev = [];
                bpm_prev = eemd_initialization(x1);
            end
        else
            bpm_prev = BPM_Ori(i-1);
        end

       %% thresholds
        th_acc_peaks = 5;
        th_NF_1 = 15;
        th_EEMD_1 = 15;
        th_NF_2 = 10;
        th_EEMD_2 = 10;
        th_off = 30;
        th_pass = 8;
        th_HRcalibration2 = 40;
        
       %% the case of no MAs
        if (ma_detect == 0)
            if (i == 1)
                bpm_track = bpm_prev;
                BPM_Ori = [BPM_Ori bpm_track];
            else
                bpm_track = bpm(find(Pxx == max(Pxx)));
                if (bpm_track < 60 || bpm_track > 190)
                    bpm_track = bpm1(find(Pxx1 == max(Pxx1)));
                    if (bpm_track < 60 || bpm_track > 190)
                        bpm_track = BPM_Ori(i-1);
                        BPM_Ori = [BPM_Ori bpm_track];
                        prev_count_1 = prev_count_1 + 1;
                        i_set_1 = [i_set_1 i];
                    else
                        BPM_Ori = [BPM_Ori bpm_track]; 
                    end
                else
                    if (abs(bpm_track - bpm_prev) > th_pass) 
                        bpm_track = [];
                        bpm_track = BPM_Ori(i-1); 
                        BPM_Ori = [BPM_Ori bpm_track];
                        prev_count_1 = prev_count_1 + 1;
                        i_set_1 = [i_set_1 i];
                    else
                        BPM_Ori = [BPM_Ori bpm_track];
                    end
                end
            end
            % HR Calibration 1
            if (prev_count_1 == 2)
                if (max(i_set_1) - min(i_set_1) == 1)
                    bpm_track = bpm1(find(Pxx1 == max(Pxx1))); 
                    if (bpm_track > 60 && bpm_track < 190 && abs(bpm_track-bpm_prev) < th_pass) 
                        BPM_Ori(i) = bpm_track;
                    else
                        BPM_Ori(i) = BPM_Ori(i-1);
                    end
                    prev_count_1 = 0;
                    i_set_1 = [];
                end
            end
            if (prev_count_1 >= 2)
                if (mean(diff(i_set_1)) ~= 1)
                    prev_count_1 = 0;
                    i_set_1 = [];
                end
            end
        end
 

       %% the case of containong MAs 
        if (ma_detect == 1)
            % accx
            dPxx2 = diff(Pxx2);
            for r = 1:length(dPxx2)-1
                if (dPxx2(r)>0 && dPxx2(r+1)<0 && Pxx2(r+1)>max(Pxx2)*0.3)
                    Pxx_peak_acc_2 = [Pxx_peak_acc_2 Pxx2(r+1)];
                    bpm_peak_acc_2 = [bpm_peak_acc_2 bpm(r+1)]; 
                end
            end
            bpm_peak_acc_2_max = bpm_peak_acc_2(find(Pxx_peak_acc_2 == max(Pxx_peak_acc_2)));
            % accy
            dPxx3 = diff(Pxx3);
            for r = 1:length(dPxx3)-1
                if (dPxx3(r)>0 && dPxx3(r+1)<0 && Pxx3(r+1)>max(Pxx3)*0.3)
                    Pxx_peak_acc_3 = [Pxx_peak_acc_3 Pxx3(r+1)];
                    bpm_peak_acc_3 = [bpm_peak_acc_3 bpm(r+1)]; 
                end
            end
            bpm_peak_acc_3_max = bpm_peak_acc_3(find(Pxx_peak_acc_3 == max(Pxx_peak_acc_3)));
            % accz
            dPxx4 = diff(Pxx4);
            for r = 1:length(dPxx4)-1
                if (dPxx4(r)>0 && dPxx4(r+1)<0 && Pxx4(r+1)>max(Pxx4)*0.3)
                    Pxx_peak_acc_4 = [Pxx_peak_acc_4 Pxx4(r+1)];
                    bpm_peak_acc_4 = [bpm_peak_acc_4 bpm(r+1)]; 
                end
            end
            bpm_peak_acc_4_max = bpm_peak_acc_4(find(Pxx_peak_acc_4 == max(Pxx_peak_acc_4)));
            % th_acc_peaks is tunable.
            if (length(bpm_peak_acc_2) > th_acc_peaks || length(bpm_peak_acc_3) > th_acc_peaks || length(bpm_peak_acc_4) > th_acc_peaks) 
                bpm_track = bpm_prev;
                BPM_Ori = [BPM_Ori bpm_track];
                prev_count_2 = prev_count_2 + 1;
                i_set_2 = [i_set_2 i];
            else
                bpm_peak_acc_max = [bpm_peak_acc_2_max bpm_peak_acc_3_max bpm_peak_acc_4_max];
                bpm_peak_acc_max = unique(bpm_peak_acc_max);

                Pxx_peak_acc = [Pxx_peak_acc_2 Pxx_peak_acc_3 Pxx_peak_acc_4];
                bpm_peak_acc = [bpm_peak_acc_2 bpm_peak_acc_3 bpm_peak_acc_4];
                bpm_peak_acc = unique(bpm_peak_acc);

                dPxx = diff(Pxx);
                for r = 1:length(dPxx)-1
                    if (dPxx(r)>0 && dPxx(r+1)<0 && Pxx(r+1)>max(Pxx)*0.3)
                        Pxx_peak = [Pxx_peak Pxx(r+1)];
                        bpm_peak = [bpm_peak bpm(r+1)]; 
                    end
                end
                num_peak_ori = length(bpm_peak);
                for p = 1:num_peak_ori
                    for q = 1:length(bpm_peak_acc)
                        if (abs(bpm_peak(p)-bpm_peak_acc(q))<=8)
                            p_loc_remove = [p_loc_remove p];
                        end
                    end
                end
                p_loc_remove = unique(p_loc_remove);
                bpm_peak(p_loc_remove) = [];
                bpm_peak = unique(bpm_peak);
                num_peak_rev = length(bpm_peak); 

                m = 21;
                n = 110;

                if (num_peak_rev == 0)
                    if (length(bpm_peak_acc) <= 3)
                        tf2 = x1; 
                        for j = 1:length(bpm_peak_acc)
                            tf2 = trapfilter(tf2,bpm_peak_acc(j));
                        end
                    else 
                        tf2 = x1;
                        for j = 1:length(bpm_peak_acc_max)
                            tf2 = trapfilter(tf2,bpm_peak_acc_max(j));
                        end
                    end
                    [Pxx_tf2 f_tf2] = periodogram(tf2,window,nfft,fs); 
                    Pxx_tf2 = Pxx_tf2(m:n); 
                    bpm_tf2 = f_tf2(m:n)*60;
                    dPxx_tf2 = diff(Pxx_tf2);
                    for k = 1:length(dPxx_tf2)-1
                        if (dPxx_tf2(k)>0 && dPxx_tf2(k+1)<0)
                            Pxx_peak_tf2 = [Pxx_peak_tf2 Pxx_tf2(k+1)];
                            bpm_peak_tf2 = [bpm_peak_tf2 bpm_tf2(k+1)]; 
                        end
                    end
                    bpm_track = bpm_peak_tf2(find(Pxx_peak_tf2 == max(Pxx_peak_tf2)));
                    bpm_track = unique(bpm_track);
                    trap_value = bpm_track; 
                    if (abs(bpm_track - bpm_prev) > th_NF_1 || bpm_track > 190) 
                        bpm_track = [];
                        bpm_track = eemd_reconstructed(x1,bpm_prev);
                        eemd_value = bpm_track; 
                        if (abs(bpm_track - bpm_prev) > th_EEMD_1) 
                            if (abs(trap_value-eemd_value) < 5 && trap_value < 190 && eemd_value < 190)
                                bpm_track = (trap_value + eemd_value)/2;
                                if (abs(bpm_track-bpm_prev) < th_off)
                                    BPM_Ori = [BPM_Ori bpm_track];
                                else
                                    bpm_track = BPM_Ori(i-1);
                                    BPM_Ori = [BPM_Ori bpm_track];
                                    prev_count_2 = prev_count_2 + 1;
                                    i_set_2 = [i_set_2 i];
                                end
                            else
                                bpm_track = [];
                                if (i ~= 1)
                                    bpm_track = BPM_Ori(i-1);
                                    BPM_Ori = [BPM_Ori bpm_track];
                                    prev_count_2 = prev_count_2 + 1;
                                    i_set_2 = [i_set_2 i];
                                else
                                    bpm_track = bpm_prev;
                                    BPM_Ori = [BPM_Ori bpm_track];
                                    prev_count_2 = prev_count_2 + 1;
                                    i_set_2 = [i_set_2 i];
                                end
                            end
                        else
                            BPM_Ori = [BPM_Ori bpm_track];
                        end
                    else
                        BPM_Ori = [BPM_Ori bpm_track];
                    end
                end

                if (num_peak_rev == 1)
                    if (abs(bpm_peak(1) - bpm_prev) < th_pass)
                        bpm_track = bpm_peak(1);
                        BPM_Ori = [BPM_Ori bpm_track];
                    else
                        bpm_track = eemd_reconstructed(x1,bpm_prev);
                        if (abs(bpm_track - bpm_prev) > th_pass) 
                            bpm_track = [];
                            if (i ~= 1)
                                bpm_track = BPM_Ori(i-1); 
                            else
                                bpm_track = bpm_prev;
                            end
                            BPM_Ori = [BPM_Ori bpm_track];
                            if (i ~= 1)
                                prev_count_2 = prev_count_2 + 1;
                                i_set_2 = [i_set_2 i];
                            end
                        else
                            BPM_Ori = [BPM_Ori bpm_track];
                        end
                    end
                end

                if (num_peak_rev > 1)
                    bpm_peak_choose = bpm_peak(find(abs(bpm_peak - bpm_prev) < th_pass)); 
                    if (length(bpm_peak_choose) == 0)
                        if (length(bpm_peak_acc) <= 3)
                            tf2 = x1;
                            for j = 1:length(bpm_peak_acc)
                                tf2 = trapfilter(tf2,bpm_peak_acc(j));
                            end
                        else 
                            tf2 = x1;
                            for j = 1:length(bpm_peak_acc_max)
                                tf2 = trapfilter(tf2,bpm_peak_acc_max(j));
                            end
                        end
                        [Pxx_tf2 f_tf2] = periodogram(tf2,window,nfft,fs); 
                        Pxx_tf2 = Pxx_tf2(m:n); 
                        bpm_tf2 = f_tf2(m:n)*60;
                        dPxx_tf2 = diff(Pxx_tf2);
                        for k = 1:length(dPxx_tf2)-1
                            if (dPxx_tf2(k)>0 && dPxx_tf2(k+1)<0)
                                Pxx_peak_tf2 = [Pxx_peak_tf2 Pxx_tf2(k+1)];
                                bpm_peak_tf2 = [bpm_peak_tf2 bpm_tf2(k+1)]; 
                            end
                        end
                        bpm_track = bpm_peak_tf2(find(Pxx_peak_tf2 == max(Pxx_peak_tf2)));
                        bpm_track = unique(bpm_track);
                        trap_value = bpm_track; 
                        if (abs(bpm_track - bpm_prev) > th_NF_2 || bpm_track > 190)
                            bpm_track = [];
                            bpm_track = eemd_reconstructed(x1,bpm_prev);
                            eemd_value = bpm_track; 
                            if (abs(bpm_track - bpm_prev) > th_EEMD_2) 
                                if (abs(trap_value-eemd_value) < 5 && trap_value < 190 && eemd_value < 190)
                                    bpm_track = (trap_value + eemd_value)/2;
                                    if (abs(bpm_track-bpm_prev) < th_off)
                                        BPM_Ori = [BPM_Ori bpm_track];
                                    else
                                        bpm_track = BPM_Ori(i-1);
                                        BPM_Ori = [BPM_Ori bpm_track];
                                        prev_count_2 = prev_count_2 + 1;
                                        i_set_2 = [i_set_2 i];
                                    end
                                else
                                    bpm_track = [];
                                    if (i ~= 1)
                                        bpm_track = BPM_Ori(i-1); 
                                        BPM_Ori = [BPM_Ori bpm_track];
                                        prev_count_2 = prev_count_2 + 1;
                                        i_set_2 = [i_set_2 i];
                                    else
                                        bpm_track = bpm_prev;
                                        BPM_Ori = [BPM_Ori bpm_track];
                                        prev_count_2 = prev_count_2 + 1;
                                        i_set_2 = [i_set_2 i];
                                    end
                                end
                            else
                                BPM_Ori = [BPM_Ori bpm_track];
                            end
                        else
                            BPM_Ori = [BPM_Ori bpm_track];
                        end
                    elseif (length(bpm_peak_choose) == 1)
                        bpm_track = bpm_peak_choose(1);
                        if (abs(bpm_track - bpm_prev) > th_pass) 
                            bpm_track = [];
                            if (i ~= 1)
                                bpm_track = BPM_Ori(i-1); 
                                BPM_Ori = [BPM_Ori bpm_track];
                            else
                                bpm_track = bpm_prev;
                                BPM_Ori = [BPM_Ori bpm_track];
                            end
                        else
                            BPM_Ori = [BPM_Ori bpm_track];
                        end
                        if (i ~= 1)
                            if (BPM_Ori(i) == BPM_Ori(i-1))
                                prev_count_2 = prev_count_2 + 1;
                                i_set_2 = [i_set_2 i];
                            end
                        end
                    else
                        diff_value = abs(bpm_peak_choose - bpm_prev); 
                        bpm_track = bpm_peak_choose(find(diff_value == min(diff_value))); 
                        if (length(bpm_track) > 1)
                            if (BPM_Ori(i-1)-BPM_Ori(i-2) >= 0)
                                bpm_track = bpm_track(find(bpm_track == max(bpm_track)));
                            else
                                bpm_track = bpm_track(find(bpm_track == min(bpm_track)));
                            end
                        end
                        if (abs(bpm_track - bpm_prev) > th_pass) 
                            bpm_track = [];
                            if (i ~= 1)
                                bpm_track = BPM_Ori(i-1); 
                                BPM_Ori = [BPM_Ori bpm_track];
                                prev_count_2 = prev_count_2 + 1;
                                i_set_2 = [i_set_2 i];
                            else
                                bpm_track = bpm_prev;
                                BPM_Ori = [BPM_Ori bpm_track];
                                prev_count_2 = prev_count_2 + 1;
                                i_set_2 = [i_set_2 i];
                            end
                        else
                            BPM_Ori = [BPM_Ori bpm_track];
                        end
                    end
                end
            end
            m = 21;
            n = 110;
            % HR Calibration 2
            if (i >= 31)
                bpm_prediction = gm_1_1_using(BPM_Ori(i-30:i-1),1);
            else
                bpm_prediction = bpm_prev+30;
            end
            if (prev_count_2 >= 3 && abs(bpm_prediction-bpm_prev) > th_pass) 
                if (max(i_set_2) - min(i_set_2) == 2)
                    if (length(bpm_peak_acc) <= 3)
                        tf2 = x1;
                        for j = 1:length(bpm_peak_acc)
                            tf2 = trapfilter(tf2,bpm_peak_acc(j));
                        end
                    else 
                        tf2 = x1;
                        for j = 1:length(bpm_peak_acc_max)
                            tf2 = trapfilter(tf2,bpm_peak_acc_max(j));
                        end
                    end
                    [Pxx_tf2 f_tf2] = periodogram(tf2,window,nfft,fs); 
                    Pxx_tf2 = Pxx_tf2(m:n); 
                    bpm_tf2 = f_tf2(m:n)*60;
                    dPxx_tf2 = diff(Pxx_tf2);
                    for k = 1:length(dPxx_tf2)-1
                        if (dPxx_tf2(k)>0 && dPxx_tf2(k+1)<0)
                            Pxx_peak_tf2 = [Pxx_peak_tf2 Pxx_tf2(k+1)];
                            bpm_peak_tf2 = [bpm_peak_tf2 bpm_tf2(k+1)]; 
                        end
                    end
                    bpm_track = bpm_peak_tf2(find(Pxx_peak_tf2 == max(Pxx_peak_tf2)));
                    bpm_track = unique(bpm_track);
                end
                
                if (bpm_track < BPM_Ori(i-1) || bpm_track > BPM_Ori(i-1) + th_HRcalibration2)
                    bpm_track = eemd_reconstructed(x1,bpm_prev);
                    if (bpm_track < BPM_Ori(i-1) || bpm_track > BPM_Ori(i-1) + th_HRcalibration2)
                        bpm_track = BPM_Ori(i-1);
                        BPM_Ori(i) = bpm_track;
                        prev_count_2 = 0;
                        i_set_2 = [];
                    else
                        BPM_Ori(i) = bpm_track;
                        prev_count_2 = 0;
                        i_set_2 = [];
                    end
                else
                    BPM_Ori(i) = bpm_track;
                    prev_count_2 = 0;
                    i_set_2 = [];
                end
            elseif (prev_count_2 > 3)
                prev_count_2 = 0;
                i_set_2 = [];
            end
            if (prev_count_2 >= 2)
                if (mean(diff(i_set_2)) ~= 1)
                    prev_count_2 = 0;
                    i_set_2 = [];
                end
            end
        end
       %%
    end
    BPM_Ori
end
BPM_Ori
toc
index = find(BPM_Ori ~= 0);
for k = 1:length(index)
    BPM_Ori_Cal(k) = BPM_Ori(index(k));
    BPM0_Cal(k) = BPM0(index(k));
end
BPM0_Cal = BPM0_Cal';
for i = 1:5
    AEE(i) = MovingAverage(BPM_Ori_Cal,BPM0_Cal,i);
end
AEE_MIN = min(AEE);
AEE_MIN