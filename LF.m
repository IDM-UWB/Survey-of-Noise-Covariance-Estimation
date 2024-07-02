function [innov]  =  LF(sys,xp,K,z)
% LF(sys,xp,K,z) linear filter with constant gain K
  [nz,N] = size(z);
  innov = zeros(nz,N);
  for i = 1:N
    innov(:,i) = z(:,i)-sys.H*xp;
    xf = xp+K*innov(:,i);           % measurement update / filtering
    xp = sys.F*xf;                 % time update / prediction
  end
end

