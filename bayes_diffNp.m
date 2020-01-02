
filepath = "E:\ketizu\code\matlab\flim\TCSPC-image-simulation-master\data\已确定的数据\论文data\数据\";
filename = "不同光子数1000组\sim_DiffNp25-1000_tau2_1000组20191220.mat";
DiffNpData = load(filepath + filename);

Nbins = DiffNpData.Nbins;
IRFmu = DiffNpData.IRFmu;
IRFsigma = DiffNpData.IRFsigma;
Nphotons = DiffNpData.Nphotons;
Ap = DiffNpData.Nphotons;
photonCountDecay = DiffNpData.photonCountDecay;
T = DiffNpData.T;
time = linspace(0,T,Nbins)';

[~,Photon_number_IRF] = MySimTCSPCimage(0, 0, IRFmu, IRFsigma, T, Nbins,0);

%% setup parameter
t_start = 0.1;  
t_end = T;            % t_start and t_end define the region of FLIM curve you wish to analyze
nexpo = 1;            % 1 : mono-expon; 2 : bi-expon


% Set up the parameter grid
if nexpo == 1
    p_min = [1,  0.8, 0.1]';   % irf shift, photon porition, tau
    p_max = [1,  0.99,  5]';      % parameter minimum, maximum value (row vectors)
    dp = [1,0.01,0.01]';      % step:主要是求寿命
elseif nexpo == 2
    p_min = [1,0.2,tau1assym,0.01,Eassym]';  % irf_shift,a1_porition,tau1,a2,tau2_porition
    p_max = [1,1,tau1assym,1,Eassym]';       % irf_shift,a1_porition,tau1,a2,tau2_porition
    dp = [1,0.0025,0.01,0.01,0.01]';         % 
end

dt = time(2) - time(1);             % size of time bin
fit_start = round(t_start/dt);      % define roi in terms of time bin #
fit_end = round(t_end/dt);

prior = 2;       % 1 : constant prior; 2 : hyperparameter
fitmethod = 3;   % 1:LS; 2:bayes(Gibbs Sampler); 3:Bayes fit (Grid)
fitvar = 2;      % if fitvar == 1, fit E (tau2/tau1) while if fitvar == 2, fit tau2

noisefrom = 1;   noiseto = 1;   %not important for bayes 
Nsample = 1;     % 拟合个数 %don't need for Grid Sampling

% NumNp = length(Nphotons);
NumNp = 1;
Ndecay = 1000;    
Taulist = zeros(NumNp, Ndecay);
bayes_meanTau = zeros(NumNp, 1);
bayes_stdTau = bayes_meanTau;

tic 
for i = 1:1
    for j = 1:1000
%     [photonCountDecay,~] = MySimTCSPCimage(tau, Nphotons, IRFmu, IRFsigma, T, Nbins,10);
%     photonCountDecay = poissrnd(double(photonCountDecay));
    
        decay = photonCountDecay(2, :, j)';
        alpha0 = 1/(sum(time.*decay)/sum(decay));                   % 平均光子数到达时间  
        alpha = 1/alpha0;                                           % 指数先验

        [pavg,sigp,pvec,post,margpost,mle] = bayes_fit(time,decay,dp,p_min,p_max,...
                                           nexpo,prior,fit_start,fit_end,alpha);
        est_tau = pavg(3);
        Taulist(i,j) = est_tau;
        bayes_meanTau(i) = mean(Taulist(i,:));
        bayes_stdTau(i) = std(Taulist(i,:));    
    end           
end
toc
