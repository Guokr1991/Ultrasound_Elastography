% Copyright, 2009, Hassan Rivaz.  All Rights reserved.
% Permission granted to download and use AM2D, Kal_LSQ and LSQ for 
% non-commercial research use only.  No further distribution or copying 
% permitted without express written consent of Hassan Rivaz. contact: rivaz@jhu.edu 
% Contact Taylor L. Jordan, Licensing Associate, 
% We might release the source codes in the future

% Please cite the following two papers if you are using this code or the
% accompanying RF-data 
% [1] 3. Rivaz, H., Boctor, E., Foroughi, P., Zellars, R., Fichtinger, G.,
% Hager, G., “Ultrasound Elastography: a Dynamic Programming Approach”, 
% IEEE Trans. Medical Imaging Oct. 2008, vol. 27 pp 1373-1377
% [2] Rivaz, H., Boctor, E., Choti, M., Hager, G., “Regularized Ultrasound
% Elastography”, IEEE Trans. Medical Imaging (submitted)

% IRF: the range of variation of axial displacement. set it to:
% [-a 0] if compression
% [0 b] if extension
% [-a b] if not sure

% IA: the range of variation of lateral displacement. set it to:
% [-a 0] if moved to the left
% [0 b] if moved to the right
% [-a b] if not sure

% midA: the seed RF-line where the displacement calculation starts from

clear all
% close all

load 'rf043.mat'
Im1 = RfDataFilt(:,:);
maxIm = max(Im1(:));
Im1 = Im1/maxIm;

load 'rf030.mat'
Im2 = RfDataFilt(:,:);
Im2 = Im2/maxIm;
clear signals_matrix

% ------------------------------------------------------- %
% ------------- set See_B_mode to see B-mode ------------ %
% ------------------------------------------------------- %
See_B_mode = 0;
if See_B_mode
    BMODE1 = log(abs(hilbert(Im1(40:end-40,10:end-10)))+.01);
    figure,  imagesc(BMODE1);colormap(gray), colorbar
    BMODE2 = log(abs(hilbert(Im2(40:end-40,10:end-10)))+.01);
    figure, imagesc(BMODE2);colormap(gray), colorbar
end
% ---------------------------------------- %
% ------------- DP Paerametes ------------ %
% ---------------------------------------- %
IRF = [-45 0];
IA = [0 0]; %Maximum allowed disparity in lateral D
alfa_DP = 0.15; % DP regularization weight
% ---------------------------------------- %
% ------------ 2D AM Paerametes ---------- %
% ---------------------------------------- %
midA = 200;
alfa = 5; % axial regularization
beta = 20; % lateral regularization
gamma = 0.005; % lateral regularization 
T = .2; % threshold for IRLS
a_t = 1; %attenuation, in dB/cm/MHz
f_0 = FileHeader.rfsd(1,1).TxFrequencyMhz/1e6; %ultrasound center freq. in MHz
f_s = 40; % sampling freq. in MHz
xx = calc_att (a_t, f_0, f_s); % to compensate for attenuation

[D1 D2 DPdisp] = AM2D(Im1, Im2, IRF, IA, midA, .15, alfa, beta, gamma, T, xx);
% the disp. of the first 40 and last 40 samples is not calculated in AM2D: 
% the disp. of the first 10 and last 10 A-lines is not calculated in AM2D: 
figure, imagesc(D1), colorbar, title('axial displacement'), colormap(hot)
figure, imagesc(D2), colorbar, title('lateral displacement'), colormap(hot)
% DPdisp is the displacement of the midA calculated by DP
% ---------------------------------------------------- %
% ------------ Calculating Strain from Disp ---------- %
% ---------------------------------------------------- %
wDIff = 93; % window length of the differentiation kernel
[strain1 strain2] = Kal_LSQ(D1(41:end-41,11:end-10),wDIff);

strain2 = -strain2((wDIff+1)/2:end-(wDIff-1)/2,:);
strain2(strain2>.03) = .03;
strain2(strain2<-.000) = -.000;

startA = 10; endA = size(Im1,2)-10;
startRF = 40+(wDIff+1)/2; endRF = size(Im1,1)-40-(wDIff+1)/2; 
xRat=(40.3-2.5)/256;
yRat=42.95/2233;
figure,  imagesc([0 xRat*(endA-startA)],yRat*[startRF endRF],strain2);colormap(gray), 
colormap(gray); colorbar, title('axial strain')
xlabel('width (mm)'); ylabel('depth (mm)')
