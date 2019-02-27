function Y = windowing_feats(XF,N,SR,TWIN,THOP)
  nwin = round(TWIN*SR);
% Always use rectangular window for now
%  if USEHANN == 1
     window = hann(nwin)';
%  else
%     window = ones(1,nwin);
%  end
   window = window/sum(window);
   XE = [zeros(N,round(nwin/2)),XF,zeros(N,round(nwin/2))];
%    XE = [XF.^2];

  hopsamps = round(THOP*SR);

  ncols = 1 + floor((size(XE,2)-nwin)/hopsamps);

  Y = zeros(N,ncols);

  winmx = repmat(window,N,1);

  for i = 1:ncols
    Y(:,i) = sqrt(mean(winmx.*XE(:,(i-1)*hopsamps + [1:nwin]),2));
%     Y(:,i) = sqrt(mean(XE(:,(i-1)*hopsamps + [1:nwin]),2));
  end
end
