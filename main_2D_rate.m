%% This is for the 2D outage/rate code
clear
%figure
%close all
f = 3.5e9; % 28 GH 
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
wavelength = c/f;
d = 3; %waveguide height
sigmanoise = 10^(-11); %noise power
ct = 500000;

M=4; % the number of users
Rtr = M*3; %target data rate 
D_w = 10;
D_all = 40; %overall area D_all X D_all  

snrdbm = [0: 10 : 40];

parfor isnr = 1:  length(snrdbm)
    P = 10^((snrdbm(isnr)-30)/10);
    eps = (2^(Rtr/2)-1)/(P/sigmanoise)/eta;
    sum_test=0;
    rate_ca_ict = zeros(ct,1);
    rate_pa_ict = zeros(ct,1);
    sum1 = 0; sum2=0; sum3=0; sum4=0;
    out_app_temp = zeros(ct,1); out_conv_app_temp = zeros(ct,1); out_ana_temp = zeros(ct,1);
     out_app_xx = zeros(ct,1); tempxx = zeros(ct,1); rate_pa_rand_ict = zeros(ct,1);
     rate_ca_rand_ict = zeros(ct,1);
    %z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
    for ict = 1 : ct 
        z_U = [rand(M, 1)*D_all-D_all/2 rand(M, 1)*D_w-D_w/2]; %locations of M users

        %pinching antennas
        PA = [z_U(:,1) zeros(M,1)];  %antenna position

        D_pa = sqrt(d^2 + (z_U(:,1)-PA(:,1).').^2 + (z_U(:,2)-PA(:,2).').^2);%pairwise distance

        H_fad = complex(1/sqrt(2)*randn(M,M),1/sqrt(2)*randn(M,M)); %complex gaussian fading
        H_pa = H_fad./D_pa.^2;

        %[upa,epa] = eig(H_pa'*H_pa);
        %lambda_pa = diag(epa); 
        %lambda_pa = sort(lambda_pa,'descend'); 
        H_pa_inv = inv(H_pa);
        beta = sqrt(1/(trace(H_pa_inv*H_pa_inv')));
        %betaxxx = sqrt(prod(lambda_pa)/sum(lambda_pa));        
        rate_pa_ict(ict) = log2(1+eta*P*beta^2/sigmanoise); 

        %random beamforming or without beamforming
        rate_pa_rand_ict(ict) = log2(1+eta*P/M* abs(H_pa(1,1))^2 ...
            /(sigmanoise+ sum(eta*P/M* abs(H_pa(1,2:end)).^2 ))); 
                
        %conventional antennas
        CAtemp = ((1:M) - (M+1)/2) * wavelength;% locations of conventional antennas
        CA = [CAtemp' zeros(M,1)];  

        D_ca = sqrt(d^2 + (z_U(:,1)-CA(:,1).').^2 + (z_U(:,2)-CA(:,2).').^2);%pairwise distance
        H_fad_ca = complex(1/sqrt(2)*randn(M,M),1/sqrt(2)*randn(M,M)); %complex gaussian fading
        H_ca = H_fad_ca./D_ca.^2;

        %[uca,eca] = eig(H_ca'*H_ca);
        %lambda_ca = diag(eca);  
        %lambda_ca = sort(lambda_ca,'descend');
        H_ca_inv = inv(H_ca);
        beta_ca = sqrt(1/(trace(H_ca_inv*H_ca_inv')));
        %beta_ca = sqrt(prod(lambda_ca)/sum(lambda_ca));
        rate_ca_ict(ict) = log2(1+eta*P*beta_ca^2/sigmanoise); 
        rate_ca_rand_ict(ict) = log2(1+eta*P/M* abs(H_ca(1,1))^2 ...
            /(sigmanoise+ sum(eta*P/M* abs(H_ca(1,2:end)).^2 )));

        %outage
        if M*rate_pa_ict(ict)<Rtr
            sum1 = sum1 + 1;
        end         
        if M*rate_ca_ict(ict)<Rtr
            sum2 = sum2 + 1;
        end    
 
    
    end 
    %mean(tempxx)
    rate_pa(isnr) = mean(rate_pa_ict);
    rate_conv(isnr) = mean(rate_ca_ict);  
    rate_ca_rand(isnr) = mean(rate_ca_rand_ict);
    rate_pa_rand(isnr) = mean(rate_pa_rand_ict);

     out_pa(isnr) = sum1/ct;
     out_conv(isnr) = sum2/ct;  
end

%plot( snrdbm,M*rate_conv, snrdbm,M*rate_pa )

plot( snrdbm,M*rate_conv, snrdbm,M*rate_pa, snrdbm,M*rate_ca_rand, snrdbm,M*rate_pa_rand )
   
%figure
%semilogy( snrdbm,out_conv,  snrdbm,min(1,out_pa)    )
 