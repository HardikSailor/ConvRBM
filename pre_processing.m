function [final] = pre_processing(s,fs)
% function for pre-processing for speech signals
     preemph = [1 -0.97];
      final = filter(preemph,1,s);
% xs = rmsilence(s,fs); % silence removal
%    s = activlev(s,fs,'n'); % level-normalization
% Signal normalization
%    ms = mean(s1); ss = std(s1); 
%    sn = (s1-ms)./ss;
%    [segments, fs] = detectVoiced1(sn,fs);
%    xs = segments{1,1};
%            xs = sn;
%          xs = rmsilence(sn,fs); % silence removal
%    fNq = fs/2;
%    Wn = [100/fNq 7500/fNq];
%    [b,a] = butter(10,Wn);
% % bandpass filtering (300 Hz - 3200Hz)
% %    [b1,a1]=potsband(fs); 
%       final=filter(b,a,s);
%            final = filter(preemph,1,final);

%    final = formant_enhance(final,fs);
   %%% preemphesis filter
%     final = filter([1 -0.97],1,xf);
% [z,p,k] = butter(10,800/fs1,'low');
% [sos,g] = zp2sos(z,p,k);	     % Convert to SOS form
% % Hd = dfilt.df2tsos(sos,g);
% xf = filtfilt(sos,g,xs);
%  
% [b,a]=butter(30,800/(fs/2),'low');
% % 
%   xf = filter(b,a,xff);
 
end
