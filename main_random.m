%% 
clear
%figure
close all
f = 28e9; % 28 GH 
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
wavelength = c/f;
d = 3; %waveguide height
eps = 1; %factional power parameter
%P = 0.1; % transmit power
sigmanoise = 10^(-12); %noise power
ct = 50000;
Rtr = 2*2; %target data rate 
 
D_all = 40; %overall area D_all X D_all  

snrdbm = [0: 10 : 40];

 
for isnr = 1: length(snrdbm)
    P = 10^((snrdbm(isnr)-30)/10);
    eps = (2^(Rtr/2)-1)/(P/sigmanoise)/eta;
    sum_test=0;
    rate_ca_ict = zeros(ct,1);
    rate_pa_ict = zeros(ct,1);
    sum1 = 0; sum2=0; sum3=0; sum4=0;
    out_app_temp = zeros(ct,1); out_conv_app_temp = zeros(ct,1); out_ana_temp = zeros(ct,1);
     out_app_xx = zeros(ct,1); rate_pa_ran_ict = zeros(ct,1);
    %z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
    for ict = 1 : ct 
        z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
        
        %pinching antennas
        pa1 = z_U(1); pa2 = z_U(2); %antenna position
        d11 = sqrt(d^2+(pa1-z_U(1))^2); d12 = sqrt(d^2+(pa1-z_U(2))^2);
        d21 = sqrt(d^2+(pa2-z_U(1))^2); d22 = sqrt(d^2+(pa2-z_U(2))^2); 

        D_pa = [d11 d21;d12 d22];
        H_fad = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %complex gaussian fading
        H_pa = H_fad./D_pa.^2;

        temp_mat = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %random gaussian fading
        [utt,ett] = eig(temp_mat'*temp_mat);
        rate_pa_ran_ict(ict) = 0; minpa=[];
        for n = 1 : 2
            temp1 = eye(2);
            beam = temp1(:,n);%utt(:,n);
            other = mod(n,2)+1;
            channelgain =  abs(H_pa(n,:)*beam)^2;
            interfegain = abs(H_pa(other,:)*beam)^2;       
            ratet = log2(1+eta*P*channelgain/(eta*P*interfegain+sigmanoise));
            rate_pa_ran_ict(ict) = rate_pa_ran_ict(ict) + ratet; 
            minpa = [minpa ratet];
        end
        rate_pa_min_ict(ict) = min(minpa);
        

        [upa,epa] = eig(H_pa'*H_pa);
        lambda_pa = diag(epa);
        beta = sqrt(prod(lambda_pa)/sum(lambda_pa));
        rate_pa_ict(ict) = log2(1+eta*P*beta^2/sigmanoise); 
                
        %conventional antennas
        ca1 = -wavelength/4; ca2 = wavelength/4; %antenna position
        d11 = sqrt(d^2+(ca1-z_U(1))^2); d12 = sqrt(d^2+(ca1-z_U(2))^2);
        d21 = sqrt(d^2+(ca2-z_U(1))^2); d22 = sqrt(d^2+(ca2-z_U(2))^2); 

        D_ca = [d11 d21;d12 d22];
        H_fad_ca = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %complex gaussian fading
        H_ca = H_fad_ca./D_ca.^2;

        [uca,eca] = eig(H_ca'*H_ca);
        lambda_ca = diag(eca);
        beta_ca = sqrt(prod(lambda_ca)/sum(lambda_ca));
        rate_ca_ict(ict) = log2(1+eta*P*beta_ca^2/sigmanoise); 
        
        %random/no beamforming
        rate_ca_ran_ict(ict) = 0;minca=[];
        for n = 1 : 2
            temp1 = eye(2);
            beam = temp1(:,n);%utt(:,n);
            other = mod(2,n)+1;
            channelgain =  abs(H_ca(n,:)*beam)^2;
            interfegain = abs(H_ca(other,:)*beam)^2;
            ratet = log2(1+eta*P*channelgain/(eta*P*interfegain+sigmanoise));
            rate_ca_ran_ict(ict) = rate_ca_ran_ict(ict) + ratet; 
            minca = [minca ratet];
        end
        rate_ca_min_ict(ict) = min(minca);

        %outage
        if 2*rate_pa_ict(ict)<Rtr
            sum1 = sum1 + 1;
        end         
        if 2*rate_ca_ict(ict)<Rtr
            sum2 = sum2 + 1;
        end    

        %approximation 
        d1 = D_pa(1,1); d2 = D_pa(1,2);
        out_app_temp(ict) = 2*d1^4*d2^4*eps/(d1^4+d2^4)*...
            (1+8*d1^4*d2^4/(d2^8-d1^8)*log(d2/d1));
        d3 = D_ca(1,1); d4 = D_ca(2,2);
        out_conv_app_temp(ict) = (d3^4+d4^4)*eps;
 
        %test
        out_app_xx(ict) = 2*d1^4*d2^4*eps/(d1^4+d2^4);
    
    end 
    
    rate_pa(isnr) = mean(rate_pa_ict);  
    rate_pa_min(isnr) = mean(rate_pa_min_ict);
    rate_ca_min(isnr) = mean(rate_ca_min_ict);
    rate_pa_ran(isnr) = mean(rate_pa_ran_ict); 
    rate_conv(isnr) = mean(rate_ca_ict); 
    rate_ca_ran(isnr) = mean(rate_ca_ran_ict); 
    out_app(isnr) = mean(out_app_temp);
    out_conv_app(isnr) = mean(out_conv_app_temp);
    out_app_x(isnr) = mean(out_app_xx);

     out_pa(isnr) = sum1/ct;
     out_conv(isnr) = sum2/ct;  
end

plot( snrdbm,2*rate_conv, snrdbm, rate_ca_ran,...
    snrdbm,2*rate_pa,snrdbm,rate_pa_ran)
 
figure

plot( snrdbm,rate_conv,   snrdbm,rate_ca_min,...
    snrdbm,rate_pa, snrdbm,rate_pa_min )
%figure
%semilogy( snrdbm,out_conv, snrdbm,min(1,out_conv_app), snrdbm,out_pa   ,snrdbm, out_app,snrdbm,out_app_x )
 