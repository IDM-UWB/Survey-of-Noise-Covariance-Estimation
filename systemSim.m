function [z]  =  systemSim(sys,N)
  % systemSim(sys,N) simulates the system from initial state [0;...;0] for N time instants
  nx = size(sys.F,1);
  nz = size(sys.H,1);
  v = chol(sys.R)'*randn(nz,N);
  w = chol(sys.Q)'*randn(nx,N);
  x = zeros(nx,N);
  z = zeros(nz,N);
  for i = 1:N
    z(:,i) = sys.H*x(:,i)+v(:,i);
    if i ~= N,x(:,i+1) = sys.F*x(:,i)+w(:,i);end
  end
end

