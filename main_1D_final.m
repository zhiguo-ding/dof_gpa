%% This is for the 1D outage/rate code
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
Rtr = 2*4; %target data rate 
 
D_all = 40; %overall area D_all X D_all  

snrdbm = [0: 10 : 40];

 
parfor isnr = 1:  length(snrdbm)
    P = 2*10^((snrdbm(isnr)-30)/10);
    eps = (2^(Rtr/2)-1)/(P/sigmanoise)/eta;
    sum_test=0;
    rate_ca_ict = zeros(ct,1);
    rate_pa_ict = zeros(ct,1);
    sum1 = 0; sum2=0; sum3=0; sum4=0;
    out_app_temp = zeros(ct,1); out_conv_app_temp = zeros(ct,1); out_ana_temp = zeros(ct,1);
     out_app_xx = zeros(ct,1); tempxx = zeros(ct,1);
    %z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
    for ict = 1 : ct 
        z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
        
        %need to order the users
        if abs(z_U(1))>abs(z_U(2))
            z_U = [z_U(2) z_U(1)];
        end 

        %pinching antennas
        pa1 = z_U(1); pa2 = z_U(2); %antenna position
        d11 = sqrt(d^2+(pa1-z_U(1))^2); d12 = sqrt(d^2+(pa1-z_U(2))^2);
        d21 = sqrt(d^2+(pa2-z_U(1))^2); d22 = sqrt(d^2+(pa2-z_U(2))^2); 

        D_pa = [d11 d21;d12 d22];
        H_fad = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %complex gaussian fading
        H_pa = H_fad./D_pa.^2;

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

        %outage
        if 2*rate_pa_ict(ict)<Rtr
            sum1 = sum1 + 1;
        end         
        if 2*rate_ca_ict(ict)<Rtr
            sum2 = sum2 + 1;
        end    

        %approximation 
        d1 = D_pa(1,1); d2 = D_pa(1,2);
        if d1>d2
            dfdf=0;
        end
        out_app_temp(ict) = min(1,2*d1^4*d2^4*eps/(d1^4+d2^4)*...
            (1+8*d1^4*d2^4/(d2^8-d1^8)*log(d2/d1)));
        d3 = D_ca(1,1); d4 = D_ca(2,2);
        out_conv_app_temp(ict) = min(1,(d3^4+d4^4)*eps);
 
        %test
        out_app_xx(ict) = min(1,2*d1^4*d2^4*eps/(d1^4+d2^4));
        %tempxx(ict) = (1+8*d1^4*d2^4/(d2^8-d1^8)*log(d2/d1));
    
    end 
    %mean(tempxx)
    rate_pa(isnr) = mean(rate_pa_ict);
    rate_conv(isnr) = mean(rate_ca_ict); 
    out_app(isnr) = mean(out_app_temp);
    out_conv_app(isnr) = mean(out_conv_app_temp);
    out_app_x(isnr) = mean(out_app_xx);

     out_pa(isnr) = sum1/ct;
     out_conv(isnr) = sum2/ct;  
end

plot( snrdbm,2*rate_conv, snrdbm,2*rate_pa )
   
%figure
%semilogy( snrdbm,out_conv, snrdbm,min(1,out_conv_app), snrdbm,min(1,out_pa)   ,snrdbm, min(1,out_app),snrdbm,min(1,out_app_x) )
 