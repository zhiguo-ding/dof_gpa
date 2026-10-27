%% This is for 1D DoF
clear
%figure
%close all
f = 3.5e9; % 28 GH 
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
wavelength = c/f;
d = 3; %waveguide height
ct = 500000;
D_all = 40; %overall area D_all X D_all  

deltavec = [0.01 : 0.01 : 0.05];

parfor isnr = 1: length(deltavec)
    delta = deltavec(isnr); 
    sum1 = 0; sum2 = 0; prank_ca_anat = zeros(ct,1);
    prank_pa_anat = zeros(ct,1);prank_pa_appt = zeros(ct,1);
    for ict = 1 : ct 
        z_U = rand(2, 1) * D_all - D_all/2; %locations of two users
        
        %z_U = [0 D_all/2];
        %pinching antennas
        pa1 = z_U(1); pa2 = z_U(2); %antenna position
        d11 = sqrt(d^2+(pa1-z_U(1))^2); d12 = sqrt(d^2+(pa1-z_U(2))^2);
        d21 = sqrt(d^2+(pa2-z_U(1))^2); d22 = sqrt(d^2+(pa2-z_U(2))^2); 

        D_pa = [d11 d21;d12 d22];
        H_fad = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %complex gaussian fading
        H_pa = H_fad./D_pa.^2;

        [upa,epa] = eig(H_pa'*H_pa);
        lambda_pa = diag(epa); 
        lambda_pa = sort(lambda_pa,'descend');
        if lambda_pa(2)/lambda_pa(1) < delta
            sum1 = sum1 + 1;
        end        
                
        %conventional antennas
        ca1 = -wavelength/4; ca2 = wavelength/4; %antenna position
        d11 = sqrt(d^2+(ca1-z_U(1))^2); d12 = sqrt(d^2+(ca1-z_U(2))^2);
        d21 = sqrt(d^2+(ca2-z_U(1))^2); d22 = sqrt(d^2+(ca2-z_U(2))^2); 

        D_ca = [d11 d21;d12 d22];
        H_fad_ca = complex(1/sqrt(2)*randn(2,2),1/sqrt(2)*randn(2,2)); %complex gaussian fading
        H_ca = H_fad_ca./D_ca.^2;

        [uca,eca] = eig(H_ca'*H_ca);
        lambda_ca = diag(eca);  
        lambda_ca = sort(lambda_ca,'descend');    
        if lambda_ca(2)/lambda_ca(1) < delta
            sum2 = sum2 + 1;
        end

        %analytical results
        d1 = D_pa(1,1); d2 = D_pa(1,2);
        prank_pa_anat(ict) = delta*(2*(d1^8+6*d1^4*d2^4+d2^8)/(d1^4+d2^4)^2 ...
            +16*d1^8*d2^8/(d2^4-d1^4)/(d1^4+d2^4)^3*log(d2^4/d1^4));
        prank_pa_appt(ict) = 2*delta*(1+4*d1^4/d2^4);
    
        temp = [D_ca(1,1) D_ca(2,1)]; %for conventional case, need to sort
        d3 = min(temp); d4 = max(temp);
        t2 = d3/d4;
        prank_ca_anat(ict) = 1 - t2^8*(1-delta)^3*(1+delta)...
            /(t2^4+delta)^2/(1+t2^4*delta)^2;
    end 
    prank_ca(isnr) = sum2/ct;
    prank_pa(isnr) = sum1/ct;

    prank_ca_ana(isnr) = mean(prank_ca_anat);
    prank_pa_ana(isnr) = mean(prank_pa_anat);
    prank_pa_app(isnr) = mean(prank_pa_appt);
    

     
end

plot( deltavec,prank_ca, deltavec,prank_ca_ana...
    ,deltavec,prank_pa, deltavec,prank_pa_ana, deltavec,prank_pa_app)
 
  