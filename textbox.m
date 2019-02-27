clc;clear all;close all;
[x Fs]=wavread('dr1_faks0_sa1_te');

% y=[
%     0 9640 
% 9640 11240 
% 11240 12783 
% 12783 14078 
% 14078 16157 
% 16157 16880 
% 16880 17103
% 17103 17587
% 17587 18760 
% 18760 19720 
% 19720 19962
% 19962 21514 
% 21514 22680
% 22680 23800 
% 23800 24104
% 24104 26280 
% 26280 28591 
% 28591 29179 
% 29179 30337
% 30337 31880 
% 31880 32500 
% 32500 33170
% 33170 33829 
% 33829 35150 
% 35150 37370
% 37370 38568 
% 38568 40546
% 40546 42357 
% 42357 45119 
% 45119 45624 
% 45624 46855 
% 46855 48680 
% 48680 49240 
% 49240 51033 
% 51033 52378 
% 52378 54500 
% 54500 55461 
% 55461 57395 
% 57395 59179 
% 59179 60600 
% 60600 63440 
% ];

[a b c] = textread('dr1_faks0_sa1_te.phn', ...
    '%d %d %s')

figure=figure(1)
plot(((1:length(x))/length(x)),x);
set(gca, 'Position', get(gca, 'OuterPosition') - ...
get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
for i=1:length(a)
% axis tight
% annotation(figure,'textbox',...
%     [((a(i)+b(i))/2)/length(x) 0.771825396825397 0.02 0.04],...
%     'String',{c{i,1}},...
%     'FitBoxToText','off');
% %     'LineStyle','none');

annotation(figure,'textbox',...
    [((a(i)+b(i))/2)/length(x) 0.92 0.02 0.04],...
    'String',{c{i,1}},...
    'FitBoxToText','off',...
    'LineStyle','none');


end

% t = annotation('textbox');
% s = t.FontSize;
% t.FontSize = 8;
% t.String='a';
% t.LineStyle='none';
% t.Position=[0.9 0.3 0.1 0.1];

%     [0.54109324009324 0.668085106382979 0.0351305361305361 0.0617021276595746],...