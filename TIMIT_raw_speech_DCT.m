%%% TIMIT raw speech log and dct features %%%%
clear all; close all; clc;

filelist = dir('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTest_40128fbankTEOAM_htk/*.rawf');

 mkdir('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTest_40128fbankTEOAMCC_htk');

for index=1:length(filelist)
   fprintf('Processing %s\n',filelist(index).name);
   fname = strcat('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTest_40128fbankTEOAM_htk/',filelist(index).name);
    htkname = strcat('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTest_40128fbankTEOAMCC_htk/',filelist(index).name);
   fbank = readhtk(fname);
   fbank = fbank';
   CC = cep_coeff(fbank,40,13,22);
%    a1 = bsxfun(@minus,CC,mean(CC));
    a1 = CC';
% a2 = bsxfun(@rdivide,a1,std(dct_data));
  % allfeats{index}=data;
%     data = data';
%     dct_data = dct(data);
%     dct_data = dct_data(1:13,:);
%     dct_data = dct_data';
   writehtk(htkname,a1,0.01,9);
    
end

clear all; close all; clc;

filelist = dir('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTrain_40128fbankTEOAM_htk/*.rawf');

 mkdir('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTrain_40128fbankTEOAMCC_htk');

for index=1:length(filelist)
   fprintf('Processing %s\n',filelist(index).name);
   fname = strcat('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTrain_40128fbankTEOAM_htk/',filelist(index).name);
    htkname = strcat('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiTrain_40128fbankTEOAMCC_htk/',filelist(index).name);
   fbank = readhtk(fname);
   fbank = fbank';
   CC = cep_coeff(fbank,40,13,22);
%    a1 = bsxfun(@minus,CC,mean(CC));
   a1 = CC';
% a2 = bsxfun(@rdivide,a1,std(dct_data));
  % allfeats{index}=data;
%     data = data';
%     dct_data = dct(data);
%     dct_data = dct_data(1:13,:);
%     dct_data = dct_data';
   writehtk(htkname,a1,0.01,9);
    
end