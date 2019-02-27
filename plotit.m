function plotit(values)
hist(values(:));
title(sprintf('mm = %g',...
    mean(mean(abs(values)))));
end