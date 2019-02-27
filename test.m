clear all
clc
%cd E:\SA_TIMIT\MFCC_0DA_BN\
cd E:\TIMIT_feats_2015\MelFB_13dct_htk\
x1=load('dr1_faks0_sa1_te.mfs.mat');
x2=load('dr1_fcjf0_sa1_tr.mfs.mat');
x3=load('dr1_fdac1_sa1_te.mfs.mat');
x4=load('dr1_fdaw0_sa1_tr.mfs.mat');
x5=load('dr1_fdml0_sa1_tr.mfs.mat');
x6=load('dr1_fecd0_sa1_tr.mfs.mat');
 x1 = x1.Xrec';
 x2 = x2.Xrec';
 x3 = x3.Xrec';
 x4 = x4.Xrec';
 x5 = x5.Xrec';
 x6 = x6.Xrec';

[p1,~]=size(x1);
[p2,~]=size(x2);
[p3,~]=size(x3);
[p4,~]=size(x4);
[p5,~]=size(x5);
[p6,~]=size(x6);
% x1 = dct(x1');
% x1 = x1';
% x2 = dct(x2');
% x2 = x2';
% x3 = dct(x3');
% x3 = x3';
% x4 = dct(x4');
% x4 = x4';
% x5 = dct(x5');
% x5 = x5';
% x6 = dct(x6');
% x6 = x6';

X=([x1;x2;x3;x4;x5;x6]);


Z1=[ones(p1,1);2*ones(p2,1);3*ones(p3,1);4*ones(p4,1);5*ones(p5,1);6*ones(p6,1)];

ydata = tsne(X,[]);


scatter(ydata(:,1), ydata(:,2), 15, Z1,'filled');