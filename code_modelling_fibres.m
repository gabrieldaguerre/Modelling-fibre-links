
subplot(3,1,1)
t = 0e-9 : 0.01e-9 : 5e-9;
t01 = 2.5e-9;
tau0=0.05e-9; %width of the pulses (0.5 for b=1 ; 0.05 for b=10)
b= 10e9 ; %bit rate in bit/s
tsep=1/b; % time between each pulse in s


R=0.7; % detector responsitivity mA/mW
NEP = 2e-9; % NEP detector mW/sqrt(Hz)

e = 1.6e-19; % electronic charge
D = -20e-12; %Dispersion
Delta_Lambda = 0.2; %source bandwidth nm
l1 = 10; %km
l2 = 50; %km
l3 = 100; %km
P0 = 1; %mW
Attenuation = 0.3;%dB/km

Pulse_In = P0 * (exp(-(t-t01).^2 / (2*tau0^2)));%mW


Delta_Tau1 = D * Delta_Lambda * l1; %dispersion L1=10km
Delta_Tau2 = D * Delta_Lambda * l2; %dispersion L2=50km
Delta_Tau3 = D * Delta_Lambda * l3; % dispersion L3=100km

Tau1 = sqrt( (tau0).^2 +  (Delta_Tau1).^2 );
Tau2 = sqrt( (tau0).^2 +  (Delta_Tau2).^2 );
Tau3 = sqrt( (tau0).^2 +  (Delta_Tau3).^2 );

A1 = tau0 / Tau1;
A2 = tau0 / Tau2;
A3 = tau0 / Tau3;


%10 km 3bit
%compute each possibility
%(0,0,0)
Pulse_Out0_10km = 0* P0 * A1 * (exp(-(t-t01).^2 / (2*Tau1^2))); %(0,0,0)
Pulse_Out0_10km_dBm = 10*log10(Pulse_Out0_10km);%dBm
Pulse_Out0_dBm_10km_afterAtt = Pulse_Out0_10km_dBm - Attenuation * l1;
Pulse_Out0_10km_mW_final = 10.^(Pulse_Out0_dBm_10km_afterAtt / 10);

Pulse_Out0_10km_mA_final = Pulse_Out0_10km_mW_final * R;

Thermal_Noise_Current_0_10km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out0_10km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out0_10km = sqrt(e*b*Pulse_Out0_10km_mA_final).*randn(size(Pulse_Out0_10km_mA_final));  % Shot noise current mA

Total_Noise_Out0_10km_mA = Thermal_Noise_Current_0_10km + Shot_Noise_Current_Out0_10km; 
Total_Noise_Out0_10km_mW = Total_Noise_Out0_10km_mA/R;

Pulse_Out0_10km_mW_final_noise = (Pulse_Out0_10km_mA_final + Total_Noise_Out0_10km_mA)/R;

SNR_Out0_10km = Pulse_Out0_10km_mW_final_noise ./ Total_Noise_Out0_10km_mW ;

% pulse (0,0,1)
Pulse_Out3_10km = P0 * A1 * (exp(-(t-t01-tsep).^2 / (2*Tau1^2))); %(0,0,1)
Pulse_Out3_10km_dBm = 10*log10(Pulse_Out3_10km);%dBm
Pulse_Out3_dBm_10km_afterAtt = Pulse_Out3_10km_dBm - Attenuation * l1;
Pulse_Out3_10km_mW_final = 10.^(Pulse_Out3_dBm_10km_afterAtt / 10);

Pulse_Out3_10km_mA_final = Pulse_Out3_10km_mW_final * R;

Thermal_Noise_Current_3_10km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out3_10km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out3_10km = sqrt(e*b*Pulse_Out3_10km_mA_final).*randn(size(Pulse_Out3_10km_mA_final));  % Shot noise current mA

Total_Noise_Out3_10km_mA = Thermal_Noise_Current_3_10km + Shot_Noise_Current_Out3_10km; 
Total_Noise_Out3_10km_mW = Total_Noise_Out3_10km_mA/R;

Pulse_Out3_10km_mW_final_noise = (Pulse_Out3_10km_mA_final + Total_Noise_Out3_10km_mA)/R;

SNR_Out3_10km = Pulse_Out3_10km_mW_final_noise ./ Total_Noise_Out3_10km_mW;

%pulse (0,1,0) central pulse
Pulse_Out1_10km = P0 * A1 * (exp(-(t-t01).^2 / (2*Tau1^2))); %(0,1,0)
Pulse_Out1_10km_dBm = 10*log10(Pulse_Out1_10km);%dBm
Pulse_Out1_dBm_10km_afterAtt = Pulse_Out1_10km_dBm - Attenuation * l1;
Pulse_Out1_10km_mW_final = 10.^(Pulse_Out1_dBm_10km_afterAtt / 10);

Pulse_Out1_10km_mA_final = Pulse_Out1_10km_mW_final * R;

Thermal_Noise_Current_1_10km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out1_10km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out1_10km = sqrt(e*b*Pulse_Out1_10km_mA_final).*randn(size(Pulse_Out1_10km_mA_final));  % Shot noise current mA

Total_Noise_Out1_10km_mA = Thermal_Noise_Current_1_10km + Shot_Noise_Current_Out1_10km; 
Total_Noise_Out1_10km_mW = Total_Noise_Out1_10km_mA / R;

Pulse_Out1_10km_mW_final_noise = (Pulse_Out1_10km_mA_final + Total_Noise_Out1_10km_mA)/R ;

SNR_Out1_10km = Pulse_Out1_10km_mW_final_noise ./ Total_Noise_Out1_10km_mW;

%pulse (1,0,0)
Pulse_Out2_10km = P0 * A1 * (exp(-(t-t01+tsep).^2 / (2*Tau1^2))); %(1,0,0)
Pulse_Out2_10km_dBm = 10*log10(Pulse_Out2_10km);%dBm
Pulse_Out2_dBm_10km_afterAtt = Pulse_Out2_10km_dBm - Attenuation * l1;
Pulse_Out2_10km_mW_final = 10.^(Pulse_Out2_dBm_10km_afterAtt / 10);

Pulse_Out2_10km_mA_final = Pulse_Out2_10km_mW_final * R;
Thermal_Noise_Current_2_10km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out2_10km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out2_10km = sqrt(e*b*Pulse_Out2_10km_mA_final).*randn(size(Pulse_Out2_10km_mA_final));  % Shot noise current mA

Total_Noise_Out2_10km_mA = Thermal_Noise_Current_2_10km + Shot_Noise_Current_Out2_10km; 
Total_Noise_Out2_10km_mW = Total_Noise_Out2_10km_mA /R ; 

Pulse_Out2_10km_mW_final_noise = (Pulse_Out2_10km_mA_final + Total_Noise_Out2_10km_mA)/R;

SNR_Out2_10km = Pulse_Out2_10km_mW_final_noise ./ Total_Noise_Out2_10km_mW; 

%pulse (1,1,0)
Pulse_Out4_10km = Pulse_Out1_10km + Pulse_Out2_10km; 
Pulse_Out4_10km_mW_final_noise = Pulse_Out1_10km_mW_final_noise + Pulse_Out2_10km_mW_final_noise;

SNR_Out4_10km = Pulse_Out4_10km_mW_final_noise ./ (Total_Noise_Out1_10km_mW + Total_Noise_Out2_10km_mW);

%pulse(0,1,1)
Pulse_Out5_10km = Pulse_Out1_10km + Pulse_Out3_10km; %(0,1,1)
Pulse_Out5_10km_mW_final_noise = Pulse_Out1_10km_mW_final_noise + Pulse_Out3_10km_mW_final_noise;

SNR_Out5_10km = Pulse_Out5_10km_mW_final_noise ./ (Total_Noise_Out1_10km_mW + Total_Noise_Out3_10km_mW);

%pulse(1,0,1)
Pulse_Out6_10km = Pulse_Out2_10km + Pulse_Out3_10km; %(1,0,1)
Pulse_Out6_10km_mW_final_noise = Pulse_Out2_10km_mW_final_noise + Pulse_Out3_10km_mW_final_noise;

SNR_Out6_10km = Pulse_Out6_10km_mW_final_noise ./ (Total_Noise_Out2_10km_mW + Total_Noise_Out3_10km_mW);

%(1,1,1)
Pulse_Out7_10km = Pulse_Out1_10km + Pulse_Out2_10km + Pulse_Out3_10km; %(1,1,1)
Pulse_Out7_10km_mW_final_noise = Pulse_Out1_10km_mW_final_noise + Pulse_Out2_10km_mW_final_noise + Pulse_Out3_10km_mW_final_noise;

SNR_Out7_10km = Pulse_Out7_10km_mW_final_noise ./ (Total_Noise_Out1_10km_mW + Total_Noise_Out2_10km_mW + Total_Noise_Out3_10km_mW );

plot(t,Pulse_Out0_10km_mW_final_noise,'.',t, Pulse_Out1_10km_mW_final_noise,'.', t, Pulse_Out2_10km_mW_final_noise,'.', t, Pulse_Out3_10km_mW_final_noise,'.', t, Pulse_Out4_10km_mW_final_noise,'.', t, Pulse_Out5_10km_mW_final_noise,'.', t, Pulse_Out6_10km_mW_final_noise,'.',t, Pulse_Out7_10km_mW_final_noise,'.'); 
title('Eye diagram of 3 pulses with noise separated by 0,1 ns (bit rate 10Gb/s) of width 0,05ns for L=10km')
legend('(0,0,0)','(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')
xlabel('time')
ylabel('Power in mW')

% plot the SNR
%plot(t, SNR_Out1_10km,'.', t, SNR_Out2_10km, '.', t, SNR_Out3_10km,'.', t, SNR_Out4_10km,'.', t, SNR_Out5_10km,'.', t, SNR_Out6_10km,'.',t, SNR_Out7_10km,'.')
%legend('(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')
 
subplot(3,1,2);


%5bits

% compute each possibility
Pulse_Out0_5bits_10km = 0* P0 * A1 * (exp(-(t-t01).^2 / (2*Tau1^2))); %(0,0,0,0,0)

Pulse_Out4_5bits_10km = P0 * A1 * (exp(-(t-t01+2*tsep).^2 / (2*Tau1^2))); %(0,0,0,0,1)
Pulse_Out3_5bits_10km = P0 * A1 * (exp(-(t-t01+tsep).^2 / (2*Tau1^2))); %(0,0,0,1,0)
Pulse_Out1_5bits_10km = P0 * A1 * (exp(-(t-t01).^2 / (2*Tau1^2))); %(0,0,1,0,0)
Pulse_Out2_5bits_10km = P0 * A1 * (exp(-(t-t01-tsep).^2 / (2*Tau1^2))); %(0,1,0,0,0)
Pulse_Out5_5bits_10km = P0 * A1 * (exp(-(t-t01-2*tsep).^2 / (2*Tau1^2))); %(1,0,0,0,0)



Pulse_Out6_5bits_10km = Pulse_Out4_5bits_10km + Pulse_Out3_5bits_10km; %(0,0,0,1,1)
Pulse_Out7_5bits_10km = Pulse_Out4_5bits_10km + Pulse_Out1_5bits_10km; %(0,0,1,0,1)
Pulse_Out8_5bits_10km = Pulse_Out4_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,0,0,1)
Pulse_Out9_5bits_10km = Pulse_Out4_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,0,0,1)

Pulse_Out10_5bits_10km = Pulse_Out3_5bits_10km + Pulse_Out1_5bits_10km; %(0,0,1,1,0)
Pulse_Out11_5bits_10km = Pulse_Out3_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,0,1,0)
Pulse_Out12_5bits_10km = Pulse_Out3_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,0,1,0)

Pulse_Out13_5bits_10km = Pulse_Out1_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,1,0,0)
Pulse_Out14_5bits_10km = Pulse_Out1_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,1,0,0)

Pulse_Out15_5bits_10km = Pulse_Out2_5bits_10km + Pulse_Out5_5bits_10km; %(1,1,0,0,0)



Pulse_Out16_5bits_10km = Pulse_Out6_5bits_10km + Pulse_Out1_5bits_10km; % (0,0,1,1,1)
Pulse_Out17_5bits_10km = Pulse_Out6_5bits_10km + Pulse_Out2_5bits_10km; % (0,1,0,1,1)
Pulse_Out18_5bits_10km = Pulse_Out6_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,0,1,1)

Pulse_Out19_5bits_10km = Pulse_Out7_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,1,0,1)
Pulse_Out20_5bits_10km = Pulse_Out7_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,1,0,1)
Pulse_Out21_5bits_10km = Pulse_Out8_5bits_10km + Pulse_Out5_5bits_10km; %(1,1,0,0,1)

Pulse_Out22_5bits_10km = Pulse_Out10_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,1,1,0)
Pulse_Out23_5bits_10km = Pulse_Out10_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,1,1,0)

Pulse_Out24_5bits_10km = Pulse_Out11_5bits_10km + Pulse_Out5_5bits_10km; %(1,1,0,1,0)
Pulse_Out25_5bits_10km = Pulse_Out8_5bits_10km + Pulse_Out1_5bits_10km; %(1,1,1,0,0)



Pulse_Out26_5bits_10km = Pulse_Out16_5bits_10km + Pulse_Out2_5bits_10km; %(0,1,1,1,1)
Pulse_Out27_5bits_10km = Pulse_Out16_5bits_10km + Pulse_Out5_5bits_10km; %(1,0,1,1,1)

Pulse_Out28_5bits_10km = Pulse_Out17_5bits_10km + Pulse_Out5_5bits_10km; %(1,1,0,1,1)

Pulse_Out29_5bits_10km = Pulse_Out19_5bits_10km + Pulse_Out5_5bits_10km; %(1,1,1,0,1)

Pulse_Out30_5bits_10km = Pulse_Out23_5bits_10km + Pulse_Out2_5bits_10km; %(1,1,1,1,0)



Pulse_Out31_5bits_10km = Pulse_Out28_5bits_10km + Pulse_Out1_5bits_10km; %(1,1,1,1,1)


%plot(t,Pulse_Out30_5bits_10km)

%plot(t,Pulse_Out0_5bits_10km, t, Pulse_Out1_5bits_10km,t, Pulse_Out2_5bits_10km,t, Pulse_Out3_5bits_10km,t, Pulse_Out4_5bits_10km, t, Pulse_Out5_5bits_10km, t, Pulse_Out6_5bits_10km,t, Pulse_Out7_5bits_10km, t, Pulse_Out8_5bits_10km, t, Pulse_Out9_5bits_10km,t, Pulse_Out10_5bits_10km, t, Pulse_Out11_5bits_10km, t, Pulse_Out12_5bits_10km, t, Pulse_Out13_5bits_10km, t, Pulse_Out14_5bits_10km, t, Pulse_Out15_5bits_10km, t, Pulse_Out16_5bits_10km, t, Pulse_Out17_5bits_10km, t, Pulse_Out18_5bits_10km, t, Pulse_Out19_5bits_10km, t, Pulse_Out20_5bits_10km, t,Pulse_Out21_5bits_10km, t, Pulse_Out22_5bits_10km, t, Pulse_Out23_5bits_10km, t, Pulse_Out24_5bits_10km, t, Pulse_Out25_5bits_10km, t, Pulse_Out26_5bits_10km,t, Pulse_Out27_5bits_10km, t, Pulse_Out28_5bits_10km, t, Pulse_Out29_5bits_10km, t, Pulse_Out30_5bits_10km, t, Pulse_Out31_5bits_10km)
%title('Eye diagram of 5 pulses separated by X ns of width Xns for L=10km')
%xlabel('time')
%ylabel('power in mW')


%50km 3bits
%compute each possibility
%(0,0,0)
Pulse_Out0_50km = 0* P0 * A2 * (exp(-(t-t01).^2 / (2*Tau2^2))); %(0,0,0)
Pulse_Out0_50km_dBm = 10*log10(Pulse_Out0_50km);%dBm
Pulse_Out0_50km_dBm_afterAtt = Pulse_Out0_50km_dBm - Attenuation * l2;
Pulse_Out0_50km_mW_final = 10.^(Pulse_Out0_50km_dBm_afterAtt / 10);

Pulse_Out0_50km_mA_final = Pulse_Out0_50km_mW_final * R;

Thermal_Noise_Current_0_50km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out0_50km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out0_50km = sqrt(e*b*Pulse_Out0_50km_mA_final).*randn(size(Pulse_Out0_50km_mA_final));  % Shot noise current mA

Total_Noise_Out0_50km_mA = Thermal_Noise_Current_0_50km + Shot_Noise_Current_Out0_50km; 
Total_Noise_Out0_50km_mW = Total_Noise_Out0_50km_mA / R;

Pulse_Out0_50km_mW_final_noise = (Pulse_Out0_50km_mA_final + Total_Noise_Out0_50km_mA)/R;

SNR_Out0_50km = Pulse_Out0_50km_mW_final_noise ./ Total_Noise_Out0_50km_mW; 
%(0,0,1)
Pulse_Out3_50km = P0 * A2 * (exp(-(t-t01-tsep).^2 / (2*Tau2^2))); %(0,0,1)
Pulse_Out3_50km_dBm = 10*log10(Pulse_Out3_50km);%dBm
Pulse_Out3_dBm_50km_afterAtt = Pulse_Out3_50km_dBm - Attenuation * l2;
Pulse_Out3_50km_mW_final = 10.^(Pulse_Out3_dBm_50km_afterAtt / 10);

Pulse_Out3_50km_mA_final = Pulse_Out3_50km_mW_final * R;

Thermal_Noise_Current_3_50km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out3_50km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out3_50km = sqrt(e*b*Pulse_Out3_50km_mA_final).*randn(size(Pulse_Out3_50km_mA_final));  % Shot noise current mA

Total_Noise_Out3_50km_mA = Thermal_Noise_Current_3_50km + Shot_Noise_Current_Out3_50km; 
Total_Noise_Out3_50km_mW = Total_Noise_Out3_50km_mA / R;

Pulse_Out3_50km_mW_final_noise = (Pulse_Out3_50km_mA_final + Total_Noise_Out3_50km_mA)/R;

SNR_Out3_50km = Pulse_Out3_50km_mW_final_noise ./ Total_Noise_Out3_50km_mW; 

%(0,1,0)
Pulse_Out1_50km = P0 * A2 * (exp(-(t-t01).^2 / (2*Tau2^2))); %(0,1,0)
Pulse_Out1_50km_dBm = 10*log10(Pulse_Out1_50km);%dBm
Pulse_Out1_dBm_50km_afterAtt = Pulse_Out1_50km_dBm - Attenuation * l2;
Pulse_Out1_50km_mW_final = 10.^(Pulse_Out1_dBm_50km_afterAtt / 10);

Pulse_Out1_50km_mA_final = Pulse_Out1_50km_mW_final * R;

Thermal_Noise_Current_1_50km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out1_50km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out1_50km = sqrt(e*b*Pulse_Out1_50km_mA_final).*randn(size(Pulse_Out1_50km_mA_final));  % Shot noise current mA

Total_Noise_Out1_50km_mA = Thermal_Noise_Current_1_50km + Shot_Noise_Current_Out1_50km; 
Total_Noise_Out1_50km_mW = Total_Noise_Out1_50km_mA / R;

Pulse_Out1_50km_mW_final_noise = (Pulse_Out1_50km_mA_final + Total_Noise_Out1_50km_mA)/R;

SNR_Out1_50km = Pulse_Out1_50km_mW_final_noise ./ Total_Noise_Out1_50km_mW ; 

%(1,0,0)
Pulse_Out2_50km = P0 * A2 * (exp(-(t-t01+tsep).^2 / (2*Tau2^2))); %(1,0,0)
Pulse_Out2_50km_dBm = 10*log10(Pulse_Out2_50km);%dBm
Pulse_Out2_dBm_50km_afterAtt = Pulse_Out2_50km_dBm - Attenuation * l2;
Pulse_Out2_50km_mW_final = 10.^(Pulse_Out2_dBm_50km_afterAtt / 10);

Pulse_Out2_50km_mA_final = Pulse_Out2_50km_mW_final * R;

Thermal_Noise_Current_2_50km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out2_50km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out2_50km = sqrt(e*b*Pulse_Out2_50km_mA_final).*randn(size(Pulse_Out2_50km_mA_final));  % Shot noise current mA

Total_Noise_Out2_50km_mA = Thermal_Noise_Current_2_50km + Shot_Noise_Current_Out2_50km; 
Total_Noise_Out2_50km_mW = Total_Noise_Out2_50km_mA ./ R; 

Pulse_Out2_50km_mW_final_noise = (Pulse_Out2_50km_mA_final + Total_Noise_Out2_50km_mA)/R;

SNR_Out2_50km = Pulse_Out2_50km_mW_final_noise ./ Total_Noise_Out2_50km_mW;
%(1,1,0)
Pulse_Out4_50km = Pulse_Out1_50km + Pulse_Out2_50km; %(1,1,0)
Pulse_Out4_50km_mW_final_noise = Pulse_Out1_50km_mW_final_noise + Pulse_Out2_50km_mW_final_noise;

SNR_Out4_50km = Pulse_Out4_50km_mW_final_noise ./ (Total_Noise_Out1_50km_mW + Total_Noise_Out2_50km_mW) ; 
%(0,1,1)
Pulse_Out5_50km = Pulse_Out1_50km + Pulse_Out3_50km; %(0,1,1)
Pulse_Out5_50km_mW_final_noise = Pulse_Out1_50km_mW_final_noise + Pulse_Out3_50km_mW_final_noise;

SNR_Out5_50km = Pulse_Out5_50km_mW_final_noise ./ (Total_Noise_Out1_50km_mW + Total_Noise_Out3_50km_mW) ; 
%(1,0,1)
Pulse_Out6_50km = Pulse_Out2_50km + Pulse_Out3_50km; %(1,0,1)
Pulse_Out6_50km_mW_final_noise = Pulse_Out2_50km_mW_final_noise + Pulse_Out3_50km_mW_final_noise;

SNR_Out6_50km = Pulse_Out6_50km_mW_final_noise ./ (Total_Noise_Out1_50km_mW + Total_Noise_Out3_50km_mW) ; 

%(1,1,1)
Pulse_Out7_50km = Pulse_Out1_50km + Pulse_Out2_50km + Pulse_Out3_50km; %(1,1,1)
Pulse_Out7_50km_mW_final_noise = Pulse_Out1_50km_mW_final_noise + Pulse_Out2_50km_mW_final_noise + Pulse_Out3_50km_mW_final_noise;

SNR_Out7_50km = Pulse_Out7_50km_mW_final_noise ./ (Total_Noise_Out1_50km_mW + Total_Noise_Out2_50km_mW + Total_Noise_Out3_50km_mW) ; 

plot(t,Pulse_Out0_50km_mW_final_noise,'.',t,Pulse_Out1_50km_mW_final_noise,'.',t,Pulse_Out2_50km_mW_final_noise,'.',t,Pulse_Out3_50km_mW_final_noise,'.',t,Pulse_Out4_50km_mW_final_noise,'.',t,Pulse_Out5_50km_mW_final_noise,'.',t,Pulse_Out6_50km_mW_final_noise,'.',t,Pulse_Out7_50km_mW_final_noise,'.')

title('Eye diagram of 3 pulses with noise separated by 0,1 ns (bit rate 10Gb/s) of width 0,05ns for L=50km')
legend('(0,0,0)','(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')
xlabel('time')
ylabel('Power in mW')

%plot SNR
%plot(t, SNR_Out1_50km, '.', t, SNR_Out2_50km,'.', t, SNR_Out3_50km,'.', t, SNR_Out4_50km,'.', t, SNR_Out5_50km,'.', t, SNR_Out6_50km,'.',t, SNR_Out7_50km,'.')
%legend('(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')
subplot(3,1,3);
%50 km 5 bits
%compute each possibility
Pulse_Out0_5bits_50km = 0* P0 * A2 * (exp(-(t-t01).^2 / (2*Tau2^2))); %(0,0,0,0,0)

Pulse_Out4_5bits_50km = P0 * A2 * (exp(-(t-t01+2*tsep).^2 / (2*Tau2^2))); %(0,0,0,0,1)
Pulse_Out3_5bits_50km = P0 * A2 * (exp(-(t-t01+tsep).^2 / (2*Tau2^2))); %(0,0,0,1,0)
Pulse_Out1_5bits_50km = P0 * A2 * (exp(-(t-t01).^2 / (2*Tau2^2))); %(0,0,1,0,0)
Pulse_Out2_5bits_50km = P0 * A2 * (exp(-(t-t01-tsep).^2 / (2*Tau2^2))); %(0,1,0,0,0)
Pulse_Out5_5bits_50km = P0 * A2 * (exp(-(t-t01-2*tsep).^2 / (2*Tau2^2))); %(1,0,0,0,0)



Pulse_Out6_5bits_50km = Pulse_Out4_5bits_50km + Pulse_Out3_5bits_50km; %(0,0,0,1,1)
Pulse_Out7_5bits_50km = Pulse_Out4_5bits_50km + Pulse_Out1_5bits_50km; %(0,0,1,0,1)
Pulse_Out8_5bits_50km = Pulse_Out4_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,0,0,1)
Pulse_Out9_5bits_50km = Pulse_Out4_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,0,0,1)

Pulse_Out10_5bits_50km = Pulse_Out3_5bits_50km + Pulse_Out1_5bits_50km; %(0,0,1,1,0)
Pulse_Out11_5bits_50km = Pulse_Out3_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,0,1,0)
Pulse_Out12_5bits_50km = Pulse_Out3_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,0,1,0)

Pulse_Out13_5bits_50km = Pulse_Out1_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,1,0,0)
Pulse_Out14_5bits_50km = Pulse_Out1_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,1,0,0)

Pulse_Out15_5bits_50km = Pulse_Out2_5bits_50km + Pulse_Out5_5bits_50km; %(1,1,0,0,0)



Pulse_Out16_5bits_50km = Pulse_Out6_5bits_50km + Pulse_Out1_5bits_50km; % (0,0,1,1,1)
Pulse_Out17_5bits_50km = Pulse_Out6_5bits_50km + Pulse_Out2_5bits_50km; % (0,1,0,1,1)
Pulse_Out18_5bits_50km = Pulse_Out6_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,0,1,1)

Pulse_Out19_5bits_50km = Pulse_Out7_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,1,0,1)
Pulse_Out20_5bits_50km = Pulse_Out7_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,1,0,1)
Pulse_Out21_5bits_50km = Pulse_Out8_5bits_50km + Pulse_Out5_5bits_50km; %(1,1,0,0,1)

Pulse_Out22_5bits_50km = Pulse_Out10_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,1,1,0)
Pulse_Out23_5bits_50km = Pulse_Out10_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,1,1,0)

Pulse_Out24_5bits_50km = Pulse_Out11_5bits_50km + Pulse_Out5_5bits_50km; %(1,1,0,1,0)
Pulse_Out25_5bits_50km = Pulse_Out8_5bits_50km + Pulse_Out1_5bits_50km; %(1,1,1,0,0)



Pulse_Out26_5bits_50km = Pulse_Out16_5bits_50km + Pulse_Out2_5bits_50km; %(0,1,1,1,1)
Pulse_Out27_5bits_50km = Pulse_Out16_5bits_50km + Pulse_Out5_5bits_50km; %(1,0,1,1,1)

Pulse_Out28_5bits_50km = Pulse_Out17_5bits_50km + Pulse_Out5_5bits_50km; %(1,1,0,1,1)

Pulse_Out29_5bits_50km = Pulse_Out19_5bits_50km + Pulse_Out5_5bits_50km; %(1,1,1,0,1)

Pulse_Out30_5bits_50km = Pulse_Out23_5bits_50km + Pulse_Out2_5bits_50km; %(1,1,1,1,0)




Pulse_Out31_5bits_50km = Pulse_Out28_5bits_50km + Pulse_Out1_5bits_50km; %(1,1,1,1,1)


%plot(t,Pulse_Out30_5bits_50km)

%plot(t,Pulse_Out0_5bits_50km, t, Pulse_Out1_5bits_50km,t, Pulse_Out2_5bits_50km,t, Pulse_Out3_5bits_50km,t, Pulse_Out4_5bits_50km, t, Pulse_Out5_5bits_50km, t, Pulse_Out6_5bits_50km,t, Pulse_Out7_5bits_50km, t, Pulse_Out8_5bits_50km, t, Pulse_Out9_5bits_50km,t, Pulse_Out10_5bits_50km, t, Pulse_Out11_5bits_50km, t, Pulse_Out12_5bits_50km, t, Pulse_Out13_5bits_50km, t, Pulse_Out14_5bits_50km, t, Pulse_Out15_5bits_50km, t, Pulse_Out16_5bits_50km, t, Pulse_Out17_5bits_50km, t, Pulse_Out18_5bits_50km, t, Pulse_Out19_5bits_50km, t, Pulse_Out20_5bits_50km, t,Pulse_Out21_5bits_50km, t, Pulse_Out22_5bits_50km, t, Pulse_Out23_5bits_50km, t, Pulse_Out24_5bits_50km, t, Pulse_Out25_5bits_50km, t, Pulse_Out26_5bits_50km,t, Pulse_Out27_5bits_50km, t, Pulse_Out28_5bits_50km, t, Pulse_Out29_5bits_50km, t, Pulse_Out30_5bits_50km, t, Pulse_Out31_5bits_50km)
%title('Eye diagram of 5 pulses separated by X ns of width Xns for L=50km')
%xlabel('time')
%ylabel('power in mW')

%100km
%(0,0,0)
Pulse_Out0_100km = 0* P0 * A3 * (exp(-(t-t01).^2 / (2*Tau3^2))); %(0,0,0)
Pulse_Out0_100km_dBm = 10*log10(Pulse_Out0_100km);%dBm
Pulse_Out0_dBm_100km_afterAtt = Pulse_Out0_100km_dBm - Attenuation * l3;
Pulse_Out0_100km_mW_final = 10.^(Pulse_Out0_dBm_100km_afterAtt / 10);

Pulse_Out0_100km_mW_final = Pulse_Out0_100km_mW_final * R;
%we have everything in mA, we can add noise now
Thermal_Noise_Current_0_100km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out0_100km_mW_final)); % thermal noise current mA

Shot_Noise_Current_Out0_100km = sqrt(e*b*Pulse_Out0_100km_mW_final).*randn(size(Pulse_Out0_100km_mW_final)); % Shot noise current mA

Total_Noise_Out0_100km_mA = Thermal_Noise_Current_0_100km + Shot_Noise_Current_Out0_100km; 
Total_Noise_Out0_100km_mW = Total_Noise_Out0_100km_mA / R;

Pulse_Out0_100km_mW_final_noise = (Pulse_Out0_100km_mW_final + Total_Noise_Out0_100km_mA)/R;

SNR_Out0_100km = Pulse_Out0_100km_mW_final_noise ./ Total_Noise_Out0_100km_mW; 
%(0,0,1)
Pulse_Out3_100km = P0 * A3 * (exp(-(t-t01-tsep).^2 / (2*Tau3^2))); %(0,0,1)
Pulse_Out3_100km_dBm = 10*log10(Pulse_Out3_100km);%dBm
Pulse_Out3_dBm_100km_afterAtt = Pulse_Out3_100km_dBm - Attenuation * l3;
Pulse_Out3_100km_mW_final = 10.^(Pulse_Out3_dBm_100km_afterAtt / 10);

Pulse_Out3_100km_mA_final = Pulse_Out3_100km_mW_final * R;
%we have everything in mA, we can add noise now
Thermal_Noise_Current_3_100km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out3_100km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out3_100km = sqrt(e*b*Pulse_Out3_100km_mA_final).*randn(size(Pulse_Out3_100km_mA_final));  % Shot noise current mA

Total_Noise_Out3_100km_mA = Thermal_Noise_Current_3_100km + Shot_Noise_Current_Out3_100km; 
Total_Noise_Out3_100km_mW = Total_Noise_Out3_100km_mA / R;

Pulse_Out3_100km_mW_final_noise = (Pulse_Out3_100km_mA_final + Total_Noise_Out3_100km_mA)/R;

SNR_Out3_100km = Pulse_Out3_100km_mW_final_noise ./ Total_Noise_Out3_100km_mW;

%(0,1,0)
Pulse_Out1_100km = P0 * A3 * (exp(-(t-t01).^2 / (2*Tau3^2))); %(0,1,0)
Pulse_Out1_100km_dBm = 10*log10(Pulse_Out1_100km);%dBm
Pulse_Out1_dBm_100km_afterAtt = Pulse_Out1_100km_dBm - Attenuation * l3;
Pulse_Out1_100km_mW_final = 10.^(Pulse_Out1_dBm_100km_afterAtt / 10);

Pulse_Out1_100km_mA_final = Pulse_Out1_100km_mW_final * R;
%we have everything in mA, we can add noise now
Thermal_Noise_Current_1_100km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out1_100km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out1_100km = sqrt(e*b*Pulse_Out1_100km_mA_final).*randn(size(Pulse_Out1_100km_mA_final)); % Shot noise current mA

Total_Noise_Out1_100km = Thermal_Noise_Current_1_100km + Shot_Noise_Current_Out1_100km; 
Total_Noise_Out1_100km_mW = Total_Noise_Out1_100km / R ; 

Pulse_Out1_100km_mW_final_noise = (Pulse_Out1_100km_mA_final + Total_Noise_Out1_100km)/R;

SNR_Out1_100km = Pulse_Out1_100km_mW_final_noise ./ Total_Noise_Out1_100km_mW;

%(1,0,0)
Pulse_Out2_100km = P0 * A3 * (exp(-(t-t01+tsep).^2 / (2*Tau3^2))); %(1,0,0)
Pulse_Out2_100km_dBm = 10*log10(Pulse_Out2_100km);%dBm
Pulse_Out2_dBm_100km_afterAtt = Pulse_Out2_100km_dBm - Attenuation * l3;
Pulse_Out2_100km_mW_final = 10.^(Pulse_Out2_dBm_100km_afterAtt / 10);

Pulse_Out2_100km_mA_final = Pulse_Out2_100km_mW_final * R;

Thermal_Noise_Current_2_100km = R * NEP * sqrt((b/2))*randn(size(Pulse_Out2_100km_mA_final)); % thermal noise current mA

Shot_Noise_Current_Out2_100km = sqrt(e*b*Pulse_Out2_100km_mA_final).*randn(size(Pulse_Out2_100km_mA_final));  % Shot noise current mA

Total_Noise_Out2_100km = Thermal_Noise_Current_2_100km + Shot_Noise_Current_Out2_100km; 
Total_Noise_Out2_100km_mW = Total_Noise_Out2_100km / R;

Pulse_Out2_100km_mW_final_noise = (Pulse_Out2_100km_mA_final + Total_Noise_Out2_100km)/R;

SNR_Out2_100km = Pulse_Out2_100km_mW_final_noise ./ Total_Noise_Out2_100km_mW;

%(1,1,0)
Pulse_Out4_100km_mW_final_noise = Pulse_Out1_100km_mW_final_noise + Pulse_Out2_100km_mW_final_noise; 
SNR_Out4_100km = Pulse_Out4_100km_mW_final_noise ./ (Total_Noise_Out1_100km_mW + Total_Noise_Out2_100km_mW) ;

%(0,1,1)
Pulse_Out5_100km_mW_final_noise = Pulse_Out1_100km_mW_final_noise + Pulse_Out3_100km_mW_final_noise; 
SNR_Out5_100km = Pulse_Out5_100km_mW_final_noise ./ (Total_Noise_Out1_100km_mW + Total_Noise_Out3_100km_mW) ;

%(1,0,1)
Pulse_Out6_100km_mW_final_noise = Pulse_Out2_100km_mW_final_noise + Pulse_Out3_100km_mW_final_noise; 
SNR_Out6_100km = Pulse_Out6_100km_mW_final_noise ./ (Total_Noise_Out2_100km_mW + Total_Noise_Out3_100km_mW) ;

%(1,1,1)
Pulse_Out7_100km_mW_final_noise = Pulse_Out1_100km_mW_final_noise + Pulse_Out2_100km_mW_final_noise + Pulse_Out3_100km_mW_final_noise; 
SNR_Out7_100km = Pulse_Out7_100km_mW_final_noise ./ (Total_Noise_Out1_100km_mW + Total_Noise_Out2_100km_mW + Total_Noise_Out3_100km_mW) ;

plot(t,Pulse_Out0_100km_mW_final_noise,'.', t, Pulse_Out1_100km_mW_final_noise,'.', t, Pulse_Out2_100km_mW_final_noise,'.', t, Pulse_Out3_100km_mW_final_noise,'.', t, Pulse_Out4_100km_mW_final_noise,'.', t, Pulse_Out5_100km_mW_final_noise,'.', t, Pulse_Out6_100km_mW_final_noise,'.', t, Pulse_Out7_100km_mW_final_noise,'.')


% plot the SNR 
%plot(t, SNR_Out1_100km,'.', t, SNR_Out2_100km,'.', t, SNR_Out3_100km,'.', t, SNR_Out4_100km,'.', t, SNR_Out5_100km,'.', t, SNR_Out6_100km,'.',t, SNR_Out7_100km,'.')
%legend('(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')

title('Eye diagram of 3 pulses with noise separated by 0,1 ns (bit rate 10Gb/s) of width 0,05ns for L=100km') 
legend('(0,0,0)','(0,1,0)','(1,0,0)', '(0,0,1)', '(1,1,0)', '(0,1,1)', '(1,0,1)', '(1,1,1)')
xlabel('time')
ylabel('power in mW')

%5bitS
%compute each possibility
Pulse_Out0_5bits_100km = 0* P0 * A3 * (exp(-(t-t01).^2 / (2*Tau3^2))); %(0,0,0,0,0)

Pulse_Out4_5bits_100km = P0 * A3 * (exp(-(t-t01+2*tsep).^2 / (2*Tau3^2))); %(0,0,0,0,1)
Pulse_Out3_5bits_100km = P0 * A3 * (exp(-(t-t01+tsep).^2 / (2*Tau3^2))); %(0,0,0,1,0)
Pulse_Out1_5bits_100km = P0 * A3 * (exp(-(t-t01).^2 / (2*Tau3^2))); %(0,0,1,0,0)
Pulse_Out2_5bits_100km = P0 * A3 * (exp(-(t-t01-tsep).^2 / (2*Tau3^2))); %(0,1,0,0,0)
Pulse_Out5_5bits_100km = P0 * A3 * (exp(-(t-t01-2*tsep).^2 / (2*Tau3^2))); %(1,0,0,0,0)



Pulse_Out6_5bits_100km = Pulse_Out4_5bits_100km + Pulse_Out3_5bits_100km; %(0,0,0,1,1)
Pulse_Out7_5bits_100km = Pulse_Out4_5bits_100km + Pulse_Out1_5bits_100km; %(0,0,1,0,1)
Pulse_Out8_5bits_100km = Pulse_Out4_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,0,0,1)
Pulse_Out9_5bits_100km = Pulse_Out4_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,0,0,1)

Pulse_Out10_5bits_100km = Pulse_Out3_5bits_100km + Pulse_Out1_5bits_100km; %(0,0,1,1,0)
Pulse_Out11_5bits_100km = Pulse_Out3_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,0,1,0)
Pulse_Out12_5bits_100km = Pulse_Out3_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,0,1,0)

Pulse_Out13_5bits_100km = Pulse_Out1_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,1,0,0)
Pulse_Out14_5bits_100km = Pulse_Out1_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,1,0,0)

Pulse_Out15_5bits_100km = Pulse_Out2_5bits_100km + Pulse_Out5_5bits_100km; %(1,1,0,0,0)



Pulse_Out16_5bits_100km = Pulse_Out6_5bits_100km + Pulse_Out1_5bits_100km; % (0,0,1,1,1)
Pulse_Out17_5bits_100km = Pulse_Out6_5bits_100km + Pulse_Out2_5bits_100km; % (0,1,0,1,1)
Pulse_Out18_5bits_100km = Pulse_Out6_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,0,1,1)

Pulse_Out19_5bits_100km = Pulse_Out7_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,1,0,1)
Pulse_Out20_5bits_100km = Pulse_Out7_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,1,0,1)
Pulse_Out21_5bits_100km = Pulse_Out8_5bits_100km + Pulse_Out5_5bits_100km; %(1,1,0,0,1)

Pulse_Out22_5bits_100km = Pulse_Out10_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,1,1,0)
Pulse_Out23_5bits_100km = Pulse_Out10_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,1,1,0)

Pulse_Out24_5bits_100km = Pulse_Out11_5bits_100km + Pulse_Out5_5bits_100km; %(1,1,0,1,0)
Pulse_Out25_5bits_100km = Pulse_Out8_5bits_100km + Pulse_Out1_5bits_100km; %(1,1,1,0,0)



Pulse_Out26_5bits_100km = Pulse_Out16_5bits_100km + Pulse_Out2_5bits_100km; %(0,1,1,1,1)
Pulse_Out27_5bits_100km = Pulse_Out16_5bits_100km + Pulse_Out5_5bits_100km; %(1,0,1,1,1)

Pulse_Out28_5bits_100km = Pulse_Out17_5bits_100km + Pulse_Out5_5bits_100km; %(1,1,0,1,1)

Pulse_Out29_5bits_100km = Pulse_Out19_5bits_100km + Pulse_Out5_5bits_100km; %(1,1,1,0,1)

Pulse_Out30_5bits_100km = Pulse_Out23_5bits_100km + Pulse_Out2_5bits_100km; %(1,1,1,1,0)




Pulse_Out31_5bits_100km = Pulse_Out28_5bits_100km + Pulse_Out1_5bits_100km; %(1,1,1,1,1)

    
%plot(t,Pulse_Out30_5bits_100km)

%%plot(t,Pulse_Out0_5bits_100km, t, Pulse_Out1_5bits_100km,t, Pulse_Out2_5bits_100km,t, Pulse_Out3_5bits_100km,t, Pulse_Out4_5bits_100km, t, Pulse_Out5_5bits_100km, t, Pulse_Out6_5bits_100km,t, Pulse_Out7_5bits_100km, t, Pulse_Out8_5bits_100km, t, Pulse_Out9_5bits_100km,t, Pulse_Out10_5bits_100km, t, Pulse_Out11_5bits_100km, t, Pulse_Out12_5bits_100km, t, Pulse_Out13_5bits_100km, t, Pulse_Out14_5bits_100km, t, Pulse_Out15_5bits_100km, t, Pulse_Out16_5bits_100km, t, Pulse_Out17_5bits_100km, t, Pulse_Out18_5bits_100km, t, Pulse_Out19_5bits_100km, t, Pulse_Out20_5bits_100km, t,Pulse_Out21_5bits_100km, t, Pulse_Out22_5bits_100km, t, Pulse_Out23_5bits_100km, t, Pulse_Out24_5bits_100km, t, Pulse_Out25_5bits_100km, t, Pulse_Out26_5bits_100km,t, Pulse_Out27_5bits_100km, t, Pulse_Out28_5bits_100km, t, Pulse_Out29_5bits_100km, t, Pulse_Out30_5bits_100km, t, Pulse_Out31_5bits_100km)
%title('Eye diagram of 5 pulses separated by X ns of width X ns for L=100km')
%xlabel('time')
%ylabel('power in mW')
