%%% generating revcor filters with specific frequencies
clear all; clc; close all;
load W3_WSJ1_60_filters;
W = W3;
% W = reshape(W(end:-1:1,:),[128,60]);
% multiple_subplots(10,6,W,60);
% figure;
%%% 86065-19 for 1406 Hz
load('C:\Users\hardik\Google Drive\database\86065u19-1-80');
x1 = Y_thief(2,:);
% plot(x1); axis tight;
%   set(gca,'Visible','off');
%    set(gca,'XTickLabel','');set(gca,'YTickLabel','');
%    tightfig;
% figure;
 plot(W(:,25));axis tight;
  set(gca,'Visible','off');
   set(gca,'XTickLabel','');set(gca,'YTickLabel','');
   tightfig;