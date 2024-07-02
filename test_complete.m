% examplary test file for initialization and execution of noise CM estimation
% algorithms introduced in the considered paper

N = 1e3+1;  % number of simulation steps
sys.R = [3 0;0 2];
sys.Q = [2 -0.5;-0.5 1];
sys.F = [0.9 0;-0.3 0.8];
sys.H = [1 0;0 1];

% initial estimate
param.xp = [0;0]; % initial state estimate
param.lags = 1; % # of time lags
param.K = [0.8 0;0 0.8]; % stable linear filter gain


% Section 3.1 - ICM PARAMETERS 
%--------------------------------------------------
ICMpar.xp = param.xp; % initial state estimate
ICMpar.lags = param.lags; % time lag of autocovariance function
ICMpar.K = param.K; % stable linear filter gain

% Section 3.2 - IOCM PARAMETERS
%--------------------------------------------------
IOCMpar.B0  =  zeros(1,4); % initial condition for estimating B

% Section 3.3 - WCM PARAMETERS
%--------------------------------------------------
WCMpar.xp = param.xp; % initial state estimate
WCMpar.Neq = param.lags+1; % # of time lags (1), i.e. # of equations (2)
WCMpar.K = param.K; % stable linear filter gain

[nz,nx] = size(sys.H); 

% basis matrices for Q
Qapr = zeros(2,2,8);
Qapr(:,:,1) = [1 0;0 0];
Qapr(:,:,2) = [0 1;0 0];
Qapr(:,:,3) = [0 0;1 0];
Qapr(:,:,4) = [0 0;0 1];
Qapr(:,:,5) = zeros(nx);
Qapr(:,:,6) = zeros(nx);
Qapr(:,:,7) = zeros(nx);
Qapr(:,:,8) = zeros(nx);

% basis matrices for R
Rapr = zeros(2,2,8);
Rapr(:,:,1) = zeros(nz);
Rapr(:,:,2) = zeros(nz);
Rapr(:,:,3) = zeros(nz);
Rapr(:,:,4) = zeros(nz);
Rapr(:,:,5) = [1 0;0 0];
Rapr(:,:,6) = [0 1;0 0];
Rapr(:,:,7) = [0 0;1 0];
Rapr(:,:,8) = [0 0;0 1];

WCMpar.Qapr = Qapr;
WCMpar.Rapr = Rapr;

% Section 3.4 - MACM PARAMETERS
%--------------------------------------------------
% implementation suited only for 2 dimensional models as in the paper

MACMpar.W1 = eye(2*nx*1); MACMpar.W2 = eye(3*nx*1); % weight matrices

% Section 3.5 - MLM PARAMETERS
%--------------------------------------------------
MLMpar.M = 50; % maximum number of EM steps
MLMpar.Qa = [3.43  0.75; 0.75 2.73]; % initial value of Q
MLMpar.Ra = [2.47 -0.58;-0.58 3.37]; % initial value of R
MLMpar.xp = param.xp; % initial state estimate
MLMpar.Pp = zeros(2); % initial state CM estimate

% Section 3.6 - DCM PARAMETERS
%--------------------------------------------------
DCMpar.xp = param.xp; % initial state estimate
DCMpar.Neq = param.lags+1; % # of time lags (1), i.e. # of equations (2)
DCMpar.K = param.K; % stable linear filter gain

% Section 3.7 - MDCM PARAMETERS
%--------------------------------------------------
MDCMpar.lags = param.lags; % time lag of autocovariance function

% Section 4.1 - CMM PARAMETERS
%--------------------------------------------------
% initial state estimate
CMMpar.xp = param.xp;
CMMpar.Pp = eye(2);

% initial estimates of Q and R
CMMpar.Q = eye(2);
CMMpar.R = eye(2);

% initial time instant for matrices estimation
CMMpar.erq = floor(N/2);

% Section 4.2 - GMBM PARAMETERS
%--------------------------------------------------
% initial state estimates
GMBMpar.xp = param.xp;                                                                                                                                                        
GMBMpar.Pp = eye(2);

% qunatised matrices for Q and R
Qquant = zeros(2,2,4);
Qquant(:,:,1) = sys.Q;                                                                                                                                                       
Qquant(:,:,2) = [2 0;0 5];                                                                                                                                                       
Qquant(:,:,3) = [2 -1;-1 2];                                                                                                                                                     
Qquant(:,:,4) = [2 1;1 3];                                                                                                                                                                                                                                                                                                                                
GMBMpar.Qquant = Qquant;

Rquant = zeros(2,2,2);
Rquant(:,:,1) = [3 0;0 1.5];                                                                                                                                                       
Rquant(:,:,2) = sys.R;                                                                                                                                                        
GMBMpar.Rquant = Rquant; 

% Section 4.3 - 3OMBM PARAMETERS
%--------------------------------------------------
% basis matrices for Q
Q1 = [1 0;0 0];
Q2 = [0 0;0 1];
Q3 = [0 -0.5;-0.5 0];
TOMBMpar.Qa = {Q1,Q2,Q3};

% initial state estimate
TOMBMpar.xp = param.xp;
TOMBMpar.Pp = eye(2);
TOMBMpar.rho = [2;1;1]; % initial weights for basis matrices
TOMBMpar.pif = [1;1;1]; % weights variance

% Section 4.4 - VBM PARAMETERS
%--------------------------------------------------
VBMpar.rho = @(N)1-min(50./(N),1e-1); % behaviour of precision factor

% parameters of inverse Gamma distribution
VBMpar.alfa = [0;0]; 
VBMpar.beta = [1;1];

VBMpar.itCnt = 2; % # of filtering iterations

% initial state estimate
VBMpar.xp = param.xp;
VBMpar.Pp = zeros(2);
param.rho = @(N)1-min(50./(N),1e-1);

% single run test

% #1  3.1 - ICM    (Q,R)
% #2  3.2 - IOCM   (Q,R)
% #3  3.3 - WCM    (Q,R)
% #4  3.4 - MACM   (Q,R)
% #5  3.5 - MLM    (Q,R)
% #6  3.6 - DCM    (Q,R)
% #7  3.7 - MDCM   (Q,R)
% #8  4.1 - CMM    (Q,R)
% #9  4.2 - GMBM   (Q,R)
% #10 4.3 - 3OMBM  (Q, )
% #11 4.4 - VBM    ( ,R)

MC = 1e0;
cntMethods = 11;
est = cell(cntMethods,MC);
for imc = 1:MC
  z = systemSim(sys,N); % simulate system for N time instants
  for i = 1:cntMethods
    switch i
      case 1
        est{i,imc} = ICM(sys,z,ICMpar);
      case 2
        est{i,imc} = IOCM(sys,z,IOCMpar);
      case 3
        est{i,imc} = WCM(sys,z,WCMpar);
      case 4
        est{i,imc} = MACM(sys,z,MACMpar);
      case 5
        est{i,imc} = MLM(sys,z,MLMpar);
      case 6
        est{i,imc} = DCM(sys,z,DCMpar);
      case 7
        est{i,imc} = MDCM(sys,z,MDCMpar);
      case 8
        est{i,imc} = CMM(sys,z,CMMpar); 
      case 9
        est{i,imc} = GMBM(sys,z,GMBMpar); 
      case 10
        est{i,imc} = TOMBM(sys,z,TOMBMpar);
        est{i,imc}.R = sys.R*NaN;
      case 11
        est{i,imc} = VBM(sys,z,VBMpar);
        est{i,imc}.Q = sys.Q*NaN;
        est{i,imc}.R = diag(est{i}.R(:,end));
    end
  end
end

methods = {'ICM','IOCM','WCM','MACM','MLM','DCM','MDCM','CMM','GMBM','3OMBM','VBM'};

for i=1:cntMethods
  fprintf('Method %s\n',methods{i});
  Q = est{i,1}.Q
  R = est{i,1}.R
end
