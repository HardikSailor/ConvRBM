clear all
clc;
close all;

x1 = load('F:\Deep Learning\testing_TSNE\SA\dr1\dr1_mwbt0_sa1_te.mfc');

x2 = load('F:\Deep Learning\testing_TSNE\SA\dr1\dr1_mwbt0_sa2_te.mfc');

x1=log10(x1');
x2 = log10(x2');
x3 = [x1 x2];
ydata1 = tsne(x3','o');
%  ydata2 = tsne(x2,[],2,39);
 %%
 
   scatter(ydata1(:,1), ydata1(:,2), 9, 'b','o');
%  hold on
%  scatter(ydata2(:,1), ydata2(:,2), 9, 'b','x');
%  hold off
 
