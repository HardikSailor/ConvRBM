function multiple_subplots(rows,columns,input,bases)
ha = tight_subplot(rows,columns,[0.01,0.01],[0.01,0.01],[0.01,0.01]);
if bases~=size(input,2)
    input = input';
end
for i=1:bases; axes(ha(i));plot(input(:,i),'k'); axis tight; end
set(ha,'XTickLabel','');set(ha,'YTickLabel','');
%     set(ha,'LooseInset',get(ha,'TightInset'));
set(ha,'Visible','off');
% set(gca,'XTickLabel','');set(gca,'YTickLabel','');
end