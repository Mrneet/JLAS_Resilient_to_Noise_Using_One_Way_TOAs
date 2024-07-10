clc;
clear; 
%close all; clc;

rng('default');
warning off

% ----- simulation setting -----
noiseOrRange = 'nse';

an=[0,0;
    0,800;
    500,800;
    700,600;
    900,400;
    700,200;
    500,0;
    0,400;
    250,800;
    250,0]';

un=[400,400]';

switch noiseOrRange
    case 'nse'
        % ******* vs. noise power config *******
        nsePwrDB = 0:5:35;
        sigma = 0.5;
        errLvlDB = 0.5;
        std_deviation_range = 10.^(nsePwrDB/10);
        NumEnsembles = 1000;
    case 'err'
        % ******* vs. sensor error config *******
        nsePwrDB = 0; 
        errLvlDB = -10:5:30;
        sigma = sqrt(10.^(errLvlDB/10)); 
        std_deviation_range = 10.^(nsePwrDB/10);
        NumEnsembles = 1000;
end

[N,M] = size(an);           % N is the dimension
K = length(nsePwrDB);       % number of noise levels
R = length(sigma);          % number of ranges

%standardized data
nse = randn(M,NumEnsembles);
err = randn(N,M,NumEnsembles);
nse = nse - mean(nse,2) *ones(1,size(nse,2));
err = err - repmat(mean(err,3),[1,1,size(err,3)]);

%calculate the movement speed v of the UN
theta_r=2*pi*randn();
vd=50*rand;
v=[vd*cos(theta_r),vd*sin(theta_r)]';

% ti
ti = 0:0.005:0.005*(M-1);

iternum=100000;

time_CFPS = zeros(R, K,NumEnsembles);
time_RIS = zeros(R, K,NumEnsembles);
pos1 = zeros(N,NumEnsembles);
pos2 = zeros(N,NumEnsembles);
v1 = zeros(N,NumEnsembles);
v2 = zeros(N,NumEnsembles);

% Calculate CRLB
crlb_un=zeros(K,1);
crlb_v=zeros(K,1);
for is = 1:R   % loop through ranges  
    disp(['Sensor position error: ',num2str(errLvlDB(is)),', ',num2str(is),'/',num2str(R),' ...']);
    
    for in = 1:K    % loop through noise powers 
        disp(['Noise power (10log(\sigma^2): ',num2str(nsePwrDB(in)),', ',num2str(in),'/',num2str(K),' ...']);
        
        Qs=sigma^2*eye(M);

        crlb_matrix=CRLB(an, un, v, std_deviation_range(in),M,sigma);
        crlb_un(in)=trace(crlb_matrix(1:2,1:2));
        crlb_v(in)=trace(crlb_matrix(3:4,3:4));

        for i = 1:NumEnsembles
            an_m = an + err(:,:,i)*sqrtm(Qs);
            beta_i= (2 * 10^-5) * rand(M,1) - 10^-5;

            tao_i = zeros(M, 1);

            beta = (2 * 10^-5) * rand() - 10^-5;

            omega = (40 * rand() - 20) * 10^-6;
           
            noise_matrix = zeros(M, length(std_deviation_range));

            for ie = 1:length(std_deviation_range)
                noise_matrix(:, ie) = sqrt(std_deviation_range(ie)) * randn(M, 1);
            end

            for ii=1:M
                tao_i(ii)=norm(un+v*ti(ii)-an(:,ii))+beta+omega*ti(ii)-beta_i(ii)+noise_matrix(ii,in);
            end

            alpha = tao_i + beta_i;
         
            tic;
            phi1_CFPS = CFPS(an_m, alpha,std_deviation_range(in), M,sigma);
            time_CFPS(is, in,i) = toc;
            pos1(:,i) = phi1_CFPS(1:2);
            v1(:,i)=phi1_CFPS(3:4);
            tic;
            [theta_RIS,it_RIS] = RIS(an_m, tao_i, alpha, std_deviation_range(in), M,iternum,sigma,beta,omega,beta_i);
            time_RIS(is, in,i) = toc;
            pos2(:,i) = theta_RIS(1:2);
            v2(:,i) = theta_RIS(3:4);
        end

        rmse_1(is,in) = sqrt(mean(sum((pos1 - repmat(un,1,NumEnsembles)).^2,1)));
        rmse_2(is,in) = sqrt(mean(sum((pos2 - repmat(un,1,NumEnsembles)).^2,1)));

        vrmse_1(is,in) = sqrt(mean(sum((v1 - repmat(v,1,NumEnsembles)).^2,1)));
        vrmse_2(is,in) = sqrt(mean(sum((v2 - repmat(v,1,NumEnsembles)).^2,1)));

        Bia_1(is,in) = mean(sqrt(sum((pos1 - repmat(un,1,NumEnsembles)).^2,1)));
        Bia_2(is,in) = mean(sqrt(sum((pos2 - repmat(un,1,NumEnsembles)).^2,1)));

        vBia_1(is,in) = mean(sqrt(sum((v1 - repmat(v,1,NumEnsembles)).^2,1)));
        vBia_2(is,in) = mean(sqrt(sum((v2 - repmat(v,1,NumEnsembles)).^2,1)));
    end
    
end

total_time_CFPS = sum(time_CFPS(:)); 
total_time_RIS = sum(time_RIS(:)); 

disp(['CFPS time：', num2str(total_time_CFPS)]);
disp(['RIS time：', num2str(total_time_RIS)]);

figure;

%CFPS
semilogy(nsePwrDB, rmse_1, 's-', 'LineWidth', 1.5);
grid on;hold on;
%RIS
semilogy(nsePwrDB, rmse_2, '^-', 'LineWidth', 1.5); 
grid on;hold on;
%crlb
semilogy(nsePwrDB, sqrt(crlb_un), '--', 'LineWidth', 1.5);
legend('CFPS', 'RIS', 'CRLB');
xlabel('10log(\sigma^2(m^2))','fontsize',12); ylabel('RMSE(p) (m)','fontsize',12);

figure;
semilogy(nsePwrDB, Bia_1, 's-', 'LineWidth', 1.5);
grid on;hold on;
semilogy(nsePwrDB, Bia_2, '^-', 'LineWidth', 1.5); 
grid on;hold on;
xlabel('10log(\sigma^2(m^2))','fontsize',12); ylabel('bias(p) (m)','fontsize',12);
legend('CFPS','RIS');

figure;

%CFPS
semilogy(nsePwrDB, vrmse_1, 's-', 'LineWidth', 1.5);
grid on;hold on;
%RIS
semilogy(nsePwrDB, vrmse_2, '^-', 'LineWidth', 1.5); 
grid on;hold on;
%crlb
semilogy(nsePwrDB, sqrt(crlb_v), '--', 'LineWidth', 1.5);
legend('CFPS', 'RIS', 'CRLB');
xlabel('10log(\sigma^2(m^2))','fontsize',12); ylabel('RMSE(v) (m/s)','fontsize',12);

figure;
semilogy(nsePwrDB, vBia_1, 's-', 'LineWidth', 1.5);
grid on;hold on;
semilogy(nsePwrDB, vBia_2, '^-', 'LineWidth', 1.5); 
grid on;hold on;
xlabel('10log(\sigma^2(m^2))','fontsize',12); ylabel('bias(v) (m/s)','fontsize',12);
legend('CFPS','RIS');








