function [est] = TOMBM(sys,z,param)
  %TOMBM (sys,z,param) Third Order Moment Bayesian Method
  %
  % 3OMBM - Section 4.3
  %
  % based on:
  % D. M. Wiberg, T. D. Powel, D. Ljungquist, "An on-line parameter
  % estimator for quick convergence and time-varying linear systems",
  % IEEE Transactions on Automatic Control, vol. 45, no. 10, pp. 1854-1863,
  % 2000.
  %
  % estimates Q
  % SYS.F, SYS.H nad SYS.R are system matrices
  % Z is nz/N matrix of measurements from N time instants
  % PARAM.XP, PARAM.PP describes initial estimate of the state and its variance
  % PARAM.RHO, PARAM.PIF initial weights for basis matrices and their variance
  % PARAM.QA basis matrices for Q
  
  N = size(z,2); % obtain number of measurements
  nx = size(sys.F,2); % obtain state dimension
  
  aN = length(param.Qa);
  
  %Initial parameters
  xp = param.xp;                %x0
  Pp = param.Pp;                %Px0
  hp = zeros(nx,aN);
  rho = param.rho; % parameter values
  pif = param.pif; % parameter values variance
  Wp = cell(1,aN);
  Wf = cell(1,aN);
  for j = 1:aN
    Wp{j} = eye(nx);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP I - COMPUTE MEANS, CMs and 3rd ORDER MOMENTS
  for i = 1:N
    e = z(:,i) - sys.H*xp; %innovation
    V = sys.R+sys.H*Pp*sys.H'; %innovation covariance matrix
    K = Pp*sys.H'/V; %filter gain
    xf = xp + K*e; %state estimate measurement update
    Pf = (eye(nx)-K*sys.H)*Pp; % state CM measurement update
    for j = 1:aN
      % filtering
      rho(j) = rho(j)+hp(:,j)'*sys.H'/V*e; % parameter values
      pif(j) = pif(j)-hp(:,j)'*sys.H'/V*sys.H*hp(:,j);%parameter values variance
      hp(:,j) = (eye(nx)-K*sys.H)*(hp(:,j)+Wp{j}*sys.H'/V*e); % par. values gain
      Wf{j} = (eye(nx)-K*sys.H)*Wp{j}*(eye(nx)-K*sys.H)'; % 3rd order moment
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP II - ESTIMATE Q
    Qh = zeros(nx);
    for j = 1:aN
      Qh = Qh+param.Qa{j}*rho(j); %merge estimate Q
    end
    % prediction
    xp = sys.F*xf; %state estimate predictive update
    Pp = sys.F*Pf*sys.F' + Qh; % predictive covariance update
    % parameters time update
    for j = 1:aN
      hp(:,j) = sys.F*hp(:,j);
      Wp{j} = sys.F*Wf{j}*sys.F'+param.Qa{j}*pif(j);
    end
  end;
  est.Q = Qh;

  
  
