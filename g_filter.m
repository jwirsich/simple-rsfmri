% filter fMRI timecourse
% conn toolbox filter used (https://www.nitrc.org/projects/conn)
% 17-02-2017 Jonathan Wirsich / Connectlab
function xfilt = g_filter(ts, TR)
%     filter data
%     1-pole low-pass filter:
%     where a = T/τ, T = the time between samples, and τ (tau) is the filter time constant.

    %Frequencies proposed by [Power, 2014, NeuroImage]
    lowfreq = 0.009;
    highfreq = 0.08;
    
%     T = TR;
%     
%     tau_low = 1/(2*pi*lowfreq); %0.01Hz -> 1/(2*pi*frequency) ~ 15.915494
%     tau_high = 1/(2*pi*highfreq); %0.01Hz -> 1/(2*pi*frequency) ~ 15.915494
%     a = T/tau_low;
%     xfilt = filter(a, [1 a-1], ts);
            
%     xfilt = transpose(reproj);
%     a2= T/tau_high;
%     xfilt = filter([1-a2 a-1], a, ts);

%     order    = 2;
%     fcutlow  = 0.009;
%     fcuthigh = 0.08;
%     %sampling rate
%     fs = 1/TR;
%     [b,a]    = butter(order,[fcutlow,fcuthigh]/(fs/2), 'bandpass');
%     xfilt        = filter(b,a,ts);

    dim = size(ts);
    xfilt = zeros(dim);
    %TODO can probablu be simplified
    for i = 1:dim(1)
        filt = [lowfreq,highfreq];
        [y,fy]=conn_filter(TR,filt,ts(i,:)');
        xfilt(i,:) = y';
    end
    
end