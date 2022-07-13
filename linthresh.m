function y = linthresh(x, x0)
  y = max(x-x0, 0)./(1-x0);
end

% x0 = 0.2;
% x = 0:0.01:1.0;
% y = linthresh(x, x0);
% plot(x, y);
