function y = softplus(x, x0, beta)
  y = 1.0/((1-x0)*beta) .* log(1 + exp(beta.*(x-x0)));
end

% x0 = 0.2;
% beta = 100.0;
% x = 0:0.01:1.0;
% y = softplus(x, x0, beta);
% plot(x, y);
