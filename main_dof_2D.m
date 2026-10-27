%% 
clear
%figure
%close all
f = 3.5e9; % 28 GH 
c = 3e8; % speed of light
eta = (c/4/pi/f)^2;
wavelength = c/f;
d = 3; %waveguide height
eps = 1; %factional power parameter
%P = 0.1; % transmit power
sigmanoise = 10^(-12); %noise power
ct = 500000;
Rtr = 2*2; %target data rate 
M = 4;% number of users
D_all = 40; %overall area D_all X D_all  
D_w = 10;

snrdbm = [0: 10 : 40];

deltavec = [0.01 : 0.01 : 0.05];

 
parfor isnr = 1: length(deltavec)
    delta = deltavec(isnr); 
    sum1 = 0; sum2 = 0;
    for ict = 1 : ct 
        z_U = [rand(M, 1)*D_all-D_all/2 rand(M, 1)*D_w-D_w/2]; %locations of M users

        %pinching antennas
        PA = [z_U(:,1) zeros(M,1)];  %antenna position

        D_pa = sqrt(d^2 + (z_U(:,1)-PA(:,1).').^2 + (z_U(:,2)-PA(:,2).').^2);%pairwise distance

        H_fad = complex(1/sqrt(2)*randn(M,M),1/sqrt(2)*randn(M,M)); %complex gaussian fading
        H_pa = H_fad./D_pa.^2;

        [upa,epa] = eig(H_pa'*H_pa);
        lambda_pa = diag(epa); 
        lambda_pa = sort(lambda_pa,'descend');
        if lambda_pa(M)/lambda_pa(1) < delta
            sum1 = sum1 + 1;
        end        
                
        %conventional antennas
        CAtemp = ((1:M) - (M+1)/2) * wavelength;% locations of conventional antennas
        CA = [CAtemp' zeros(M,1)];  

        D_ca = sqrt(d^2 + (z_U(:,1)-CA(:,1).').^2 + (z_U(:,2)-CA(:,2).').^2);%pairwise distance
        H_fad_ca = complex(1/sqrt(2)*randn(M,M),1/sqrt(2)*randn(M,M)); %complex gaussian fading
        H_ca = H_fad_ca./D_ca.^2;

        [uca,eca] = eig(H_ca'*H_ca);
        lambda_ca = diag(eca);  
        lambda_ca = sort(lambda_ca,'descend');    
        if lambda_ca(M)/lambda_ca(1) < delta
            sum2 = sum2 + 1;
        end
 
    end 
    prank_ca(isnr) = sum2/ct;
    prank_pa(isnr) = sum1/ct; 
    

     
end

plot( deltavec,prank_ca, deltavec,prank_pa )
 
  