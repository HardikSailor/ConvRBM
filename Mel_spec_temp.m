function [feats] = Mel_spec_temp(sig,fs,winlen,winhop,nfilter)
x_frames = enframe(sig,hamming(winlen),winhop);
STFT = abs(fft(x_frames',512));
STFT = STFT(1:257,:);
% STFT = (abs(spectrogram(sig,winlen,winhop,512,fs)));
% STFT = (STFT);
 melbank = filtbankm(nfilter,512,fs,0,fs/2,'m');
 Melspec = STFT'*melbank';
 Melspec = log(Melspec);
 feats = (Melspec');
% feats = STFT;
end