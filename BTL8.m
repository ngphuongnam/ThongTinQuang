clear all
P = -10:2:20;
chaoslen = [10 100 500];

error = 0;
sumbit = 0;
BER = zeros(length(chaoslen),length(P));
loop = 100;
for k = 1:length(P)
P_LED = 10^(P(k)/10);  %transmitted optical power by individual LED
    for m = 1:length(chaoslen)
    len = chaoslen(m); %length of chaos sequence
        for i=1:loop
            for i=1:1
                x1(i) = randi(4)-1;
                x2(i) = randi(4)-1;
                x3(i) = randi(4)-1;
                x4(i) = randi(4)-1;
                x5(i) = randi(4)-1;
                x6(i) = randi(4)-1;
                x7(i) = randi(4)-1;
                x8(i) = randi(4)-1;
            end
            %*************QPSK MODULATION ******************
            B1 = qammod(double(x1),4,-pi/2);
            B2 = qammod(double(x2),4,-pi/2);
            B3 = qammod(double(x3),4,-pi/2);
            B4 = qammod(double(x4),4,-pi/2);
            B5 = qammod(double(x5),4,-pi/2);
            B6 = qammod(double(x6),4,-pi/2);
            B7 = qammod(double(x7),4,-pi/2);
            B8 = qammod(double(x8),4,-pi/2);

            %*************CHAOSTIC GEN ******************
            chaosticSeq1 = zeros(1,len);
            chaosticSeq2 = zeros(1,len);
            chaosticSeq3 = zeros(1,len);
            chaosticSeq4 = zeros(1,len);
            chaosticSeq5 = zeros(1,len);
            chaosticSeq6 = zeros(1,len);
            chaosticSeq7 = zeros(1,len);
            chaosticSeq8 = zeros(1,len);
            
            chaosticSeq1(1) = 0.1;
            chaosticSeq2(1) = 0.2;
            chaosticSeq3(1) = 0.3;
            chaosticSeq4(1) = 0.4;
            chaosticSeq5(1) = 0.6;
            chaosticSeq6(1) = 0.7;
            chaosticSeq7(1) = 0.8;
            chaosticSeq8(1) = 0.9;
            r = 4; % r parameter for chaotic regime

            for j=1:len-1
               chaosticSeq1(j+1) = r*chaosticSeq1(j)*(1-chaosticSeq1(j));
               chaosticSeq2(j+1) = r*chaosticSeq2(j)*(1-chaosticSeq2(j));
               chaosticSeq3(j+1) = r*chaosticSeq3(j)*(1-chaosticSeq3(j));
               chaosticSeq4(j+1) = r*chaosticSeq4(j)*(1-chaosticSeq4(j));
               chaosticSeq5(j+1) = r*chaosticSeq5(j)*(1-chaosticSeq5(j));
               chaosticSeq6(j+1) = r*chaosticSeq6(j)*(1-chaosticSeq6(j));
               chaosticSeq7(j+1) = r*chaosticSeq7(j)*(1-chaosticSeq7(j));
               chaosticSeq8(j+1) = r*chaosticSeq8(j)*(1-chaosticSeq8(j));
            end

            for j=1:length(B1)
                C1 = B1(j)*chaosticSeq1;
                C2 = B2(j)*chaosticSeq2;
                C3 = B3(j)*chaosticSeq3;
                C4 = B4(j)*chaosticSeq4;
                C5 = B5(j)*chaosticSeq5;
                C6 = B6(j)*chaosticSeq6;
                C7 = B7(j)*chaosticSeq7;
                C8 = B8(j)*chaosticSeq8;
            end
            %*************SUM ******************
            D = C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8;

            %*************INTERLEAVING ******************
            randintrlv(D,2048);

            %********** LED **************
            %Calc Gain of room
            theta = 70;% semi-angle at half power
            ml=-log10(2)/log10(cosd(theta));%Lambertian order of emission
            Adet=1e-4; %detector physical area of a PD
            lx=5; ly=5; lz=3;% room dimension in meter
            h=2.15;%the distance between source and receiver plane
            [XT,YT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4]);% position of LED; it is assumed all LEDs are located at same point for
            Nx=lx*5; Ny=ly*5;% number of grid in the receiver plane
            x=linspace(-lx/2,lx/2,Nx);
            y=linspace(-ly/2,ly/2,Ny);
            [XR,YR]=meshgrid(x,y);
            D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);% distance vector from source 
            cosphi_A1=h./D1;% angle vector
            %Calc transmit power
            nLED=60;% number of LED array nLED*nLED
            P_total=nLED*nLED*P_LED;%Total transmitted power

            %***********VLC CHANNEL ***************
            H_A1=(ml+1)*Adet.*cosphi_A1.^(ml+1)./(2*pi.*D1.^2); % FINAL GAIN OF CHANNEL
            %***********VLC CHANNEL ***************


            %***********PD***************
            %Calc reiceive power
            Ts=1;%gain of an optical filter; ignore if no filter is used
            index=1.5;%refractive index of a lens at a PD; ignore if no lens is used
            FOV=70;%FOV of a receiver
            G_Con=(index^2)/(sind(FOV).^2);%gain of an optical concentrator; ignore if no lens is used
            P_rec=P_total.*H_A1.*Ts.*G_Con;

            n = awgn(D,P_rec(9,9)/10);
            N = sqrt(P_rec(9,9)/10)*(randn(1,len) +1j*randn(1,len));
            D = D.*P_rec(9,9) + N;

            %*************DEINTERLEAVING ******************
            randdeintrlv(D,2048);

            %*************CHAOSTIC GEN ******************
            E1 = D.*chaosticSeq1;
            E2 = D.*chaosticSeq2;
            E3 = D.*chaosticSeq3;
            E4 = D.*chaosticSeq4;
            E5 = D.*chaosticSeq5;
            E6 = D.*chaosticSeq6;
            E7 = D.*chaosticSeq7;
            E8 = D.*chaosticSeq8;
            
            E1 = sum(E1)/len;
            E2 = sum(E2)/len;
            E3 = sum(E3)/len;
            E4 = sum(E4)/len;
            E5 = sum(E5)/len;
            E6 = sum(E6)/len;
            E7 = sum(E7)/len;
            E8 = sum(E8)/len;
            %*************QPSK DEMODULATION ******************
            S1 = qamdemod(double(E1),4,-pi/2);
            S2 = qamdemod(double(E2),4,-pi/2);
            S3 = qamdemod(double(E3),4,-pi/2);
            S4 = qamdemod(double(E4),4,-pi/2);
            S5 = qamdemod(double(E5),4,-pi/2);
            S6 = qamdemod(double(E6),4,-pi/2);
            S7 = qamdemod(double(E7),4,-pi/2);
            S8 = qamdemod(double(E8),4,-pi/2);

            %*************CALCULATE BER ******************
            sumbit = sumbit + 8;
            trans = [de2bi(x1,2) de2bi(x2,2) de2bi(x3,2) de2bi(x4,2) de2bi(x5,2) de2bi(x6,2) de2bi(x7,2) de2bi(x8,2)];
            recei = [de2bi(S1,2) de2bi(S2,2) de2bi(S3,2) de2bi(S4,2) de2bi(S5,2) de2bi(S6,2) de2bi(S7,2) de2bi(S8,2)];
            diff = numel(find(recei~=trans));
            error = error + diff;
        end
    BER(m,k) = error/sumbit;
    end
end
hold on;
plot(P,BER(1,:),'b--o',P,BER(2,:),'g',P,BER(3,:),'c*');
legend('ChaosSeq=10','ChaosSeq=100','ChaosSeq=500');
hold off;