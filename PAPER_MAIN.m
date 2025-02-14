%%%% PAPER 2021 %%%%%%%%%%
clc % limpa comm win
clear % limpa mem
close all % fecha win

addpath('G:\Meu Drive\Projeto Lagartos\Funções utilizadas\DisciplinaMatlab\CircStat')
cd     ('G:\Meu Drive\Projeto Lagartos\ANALISES')
load 'Animais-estimulo-aversivo.mat'

%% PARAMETROS
DATA = AnimalGato0502;
srate = 500; % frequencia de amostragem em Hz
dt = 1/srate; % Passo em segundos
time_vector = dt:dt:length(DATA)/srate; % Vetor de tempo em segundos
LFP = DATA;

%TFD
frequency_vector = 1:0.1:50;
frequency_width = 1;
clear TFDcont
TFDcont = zeros(length(frequency_vector),length(LFP));
for j=1:length(frequency_vector)
    filtrado = eegfilt(LFP,srate,frequency_vector(j),...
        frequency_vector(j)+frequency_width);
    TFDcont(j,:)= abs(hilbert(filtrado ));
end

%FILTRA O SINAL
delta = eegfilt(LFP,srate,2,6);
theta = eegfilt(LFP,srate,7,12);
beta  = eegfilt(LFP,srate,15,25);
gamma = eegfilt(LFP,srate,36,40);
%OBTEM O ENVELOPE DE AMPLITUDE
deltaAmp = abs(hilbert(delta));
thetaAmp = abs(hilbert(theta));
betaAmp = abs(hilbert(beta));
gammaAmp = abs(hilbert(gamma));

Bands={'Delta','Theta','Beta','Gamma'};

%% plot
t1=0;
t2=t1+100;
close all
% figure('units','normalized','outerposition',[0 0 1 1])
fig1 = figure(1);clf
% set(gcf, 'Position', get(0, 'Screensize'))
set(gcf,'color','white')
subplot(4,2,[1 3])
plot(time_vector,LFP,'k','linewidth',1)
title('RAW LFP signal - CD','fontsize',16)
set(gca,'fontsize',12)
xlim([t1 t2])
ylim([-0.4 0.4])
x11 = gca;
ylabel('Voltage (mV)')
xlabel('Time (s)')
box on


subplot(4,2,[5 7])
imagesc(time_vector,frequency_vector+frequency_width/2,(TFDcont))
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([t1 t2])
ylim([0 40])
caxis([0 0.05])
set(gca,'fontsize',12)
box on
title('Time Frequency Decomposition','fontsize',16)
x12 = gca;


subplot(4,2,2)
var=0.2;
plot(time_vector,delta,'b','linewidth',1)
hold on
plot(time_vector,deltaAmp,'b','linewidth',1)
xlim([95 100])
ylim([-0.1 0.1])
set(gca,'YTick',[-0.6:var:0],'YTickLabel',flip(Bands))
box off
subplot(4,2,4)
plot(time_vector,theta-1*var,'c-','linewidth',1)
hold on
plot(time_vector,thetaAmp-1*var,'c-','linewidth',1)
xlim([95 100])
ylim([-0.25 -0.15])
set(gca,'YTick',[-0.6:var:0],'YTickLabel',flip(Bands))
box off
subplot(4,2,6)
plot(time_vector,beta-2*var,'g-','linewidth',1)
hold on
plot(time_vector,betaAmp-2*var,'g-','linewidth',1)
xlim([95 100])
ylim([-0.42 -0.38])
set(gca,'YTick',[-0.6:var:0],'YTickLabel',flip(Bands))
box off
subplot(4,2,8)
plot(time_vector,gamma-3*var,'m-','linewidth',1)
hold on
plot(time_vector,gammaAmp-3*var,'m-','linewidth',1)
title('Filtered signals','fontsize',16)
set(gca,'fontsize',12)
xlim([95 100])
ylim([-0.62 -0.58])
x13 = gca;
xlabel('Time (s)')
set(gca,'YTick',[-0.6:var:0],'YTickLabel',flip(Bands))
box off
%%
xlim([1 300])

%% INSPEÇÃO GERAL DO SINAL DE GRUPO
clc
clear
close all
addpath('G:\Meu Drive\Projeto Lagartos\Funções utilizadas\DisciplinaMatlab\CircStat')
cd     ('G:\Meu Drive\Projeto Lagartos\ANALISES')
load 'Animais-estimulo-aversivo.mat'
load 'Animais-controle.mat'
srate=500;

%% Tratamento do ruído dia 01
G0201 = AnimalGato0201([1:112*srate, ...
    128*srate:161.6*srate, 161.7*srate:227*srate, ...
    233*srate:514.5*srate, 514.6*srate:532*srate, ...
    552*srate:568*srate, 582*srate:end]);
G0301 = AnimalGato0301([1:3*srate, 7*srate:152*srate, 158*srate:292*srate, ...
    298*srate:313.5*srate, 313.7*srate:527.25*srate, ...
    528.25*srate:end]);
G0501 = AnimalGato0501([1:102.3*srate, 102.5*srate:112.3*srate, ...
    112.6*srate:130.9*srate, 131*srate:264.2*srate, ...
    264.25*srate:556.8*srate, 556.9*srate:end]);
G0701 = AnimalGato0701([1:28.6*srate, 28.8*srate:323*srate, 332*srate:579.9*srate, ...
    580*srate:end]);
G1101 = AnimalGato1101([1:382*srate, 393*srate:end]);

G0202 = AnimalGato0202([1:52*srate, 62*srate:93*srate, 95*srate:133*srate, ...
    152*srate:187*srate, 190*srate:197*srate, ...
    200*srate:201.8*srate, 202.4*srate, 202.5*srate:324.8*srate, ...
    325*srate:482*srate, 488*srate:507*srate, 518*srate:549*srate, ...
    556*srate:end]);
G0302 = AnimalGato0302([1:2.5*srate, 2.7*srate:4.85*srate, 4.95*srate:204.4*srate,...
    204.5*srate:430.7*srate, 430.8*srate:end]);
G0502 = AnimalGato0502([1:72.5*srate, 72.6*srate:184*srate, 184.1*srate:304.3*srate,...
    304.5*srate:479.65*srate, 479.75*srate:end]);
G0702 = AnimalGato0702([5*srate:13*srate, 20*srate:23*srate, 28*srate:247.8*srate, ...
    248*srate:553.6*srate, 553.7*srate:end]);
G1102 = AnimalGato1102([1:21*srate, 24*srate:97*srate, 103*srate:137*srate, 143*srate:192*srate, ...
    198*srate:217*srate, 227*srate:287*srate, 293*srate:338*srate, ...
    347*srate:353*srate, 362*srate:370*srate,378*srate:403*srate, ...
    407*srate:427*srate, 432*srate:447*srate, 452*srate:456*srate, ...
    463*srate:486*srate,495*srate:end]);

XAMBG = {G0201;G0301;G0501;G0701;G1101};
XEXPG = {G0202;G0302;G0502;G0702;G1102};


C0101 = AnimalControle0101([12*srate:end]);
C0401 = AnimalControle0401([1:453*srate, 457*srate:end]);
C0401 = C0401([1:163*srate, 183*srate:228*srate, 240*srate:474*srate, ...
    476*srate:end])
C0601 = AnimalControle0601([1:end]);
C0901 = AnimalControle0901([1:3.5*srate,4*srate:411*srate, 414*srate:end]);
C1001 = AnimalControle1001([1:104*srate, 108*srate:263*srate, 269*srate:end]);

C0102 = AnimalControle0102([12*srate:end]);
C0402 = AnimalControle0402([1:312*srate, 317*srate:528*srate, 535*srate:end]);
C0402 = C0402([5*srate:52*srate, 54*srate:286*srate, 289*srate:end])
C0602 = AnimalControle0602([1:end]);
C0902 = AnimalControle0902([12*srate:78.5*srate, 81.5*srate:146*srate, ...
    158*srate:201*srate, 207*srate:231*srate, ...
    235*srate:240*srate, 244*srate:255*srate, ...
    262*srate:268*srate, 277*srate:288*srate, ...
    312*srate:318*srate, 335*srate:345*srate, ...
    349*srate:468*srate, 472*srate:492*srate, ...
    495*srate:519*srate, 525*srate:543*srate, ...
    550*srate:570*srate, 574*srate:579*srate, ...
    583*srate:end]);
C1002 = AnimalControle1002([1:60*srate, 63*srate:130*srate, 133*srate:239*srate, ...
    246*srate:286*srate, 293*srate:358*srate, 363*srate:end]);
%sem ruído
XAMBC = {C0101;C0401;C0601;C0901;C1001};
XEXPC = {C0102;C0402;C0602;C0902;C1002};

%% Tratamento do ruído dia 02
G0201d = AnimalGato0201d2 ([1:45.1*srate, 45.3*srate:54.2*srate, ...
    54.6*srate:133.6*srate, 134.4*srate:205*srate,...
    210*srate:288*srate, 293*srate:362*srate,...
    364*srate:396.3*srate, 396.6*srate:435.5*srate,...
    436*srate:455*srate, 455.5*srate:488*srate,...
    488.5*srate:539*srate, 539.5*srate:591.5*srate,...
    593*srate:end]);
G0301d = AnimalGato0301d2([1:15*srate, 15.3*srate:58.6*srate,...
    58.9*srate:74.9*srate, 75.2*srate:94.8*srate,...
    95*srate:100.8*srate, 101*srate:132*srate,...
    138*srate:206*srate, 207*srate:211.7*srate, ...
    212.3*srate:225.5*srate, 226*srate:244.3*srate,...
    244.5*srate:298.8*srate, 299*srate:338.9*srate,...
    339*srate:359.5*srate, 360*srate:400.7*srate, ...
    401*srate:509.5*srate, 510*srate:541.2*srate, ...
    541.4*srate:end]);
G0501d = AnimalGato0501d2([1:27.6*srate, 27.8*srate:38*srate, ...
    38.2*srate:49.8*srate, 50*srate:56.7*srate, ...
    56.8*srate:272*srate, 282*srate:298.8*srate, ...
    298.9*srate:352.6*srate, 352.9*srate:369.2*srate, ...
    369.4*srate:456.1*srate, 456.2*srate:471.6*srate, ...
    471.7*srate:558.1*srate, 558.2*srate:end]);
G0701d = AnimalGato0701d2([1:214.9*srate, 215.1*srate:242.7*srate, ...
    242.8*srate:247.1*srate, 247.2*srate:320.1*srate,...
    320.2*srate:332.9*srate, 333*srate:364*srate, ...
    364.2*srate:395.7*srate, 395.8*srate:565.4*srate, ...
    565.5*srate:595.7*srate, 595.8*srate:end]);
G1101d = AnimalGato1101d2([1:96*srate, 96.8*srate:168.2*srate, ...
    169*srate:245*srate, 249*srate:276.8*srate, ...
    279*srate:301.6*srate, 301.8*srate:339*srate, ...
    353*srate:377*srate, 382*srate:397.5*srate, ...
    406.3*srate:413.5*srate, 416*srate:417*srate, ...
    419*srate:430*srate, 432*srate:446.3*srate, ...
    446.5*srate:498*srate, 505*srate:571*srate, 579*srate:end]);

G0202d = AnimalGato0202d2([1*srate:50.1*srate, 50.5*srate:63*srate, ...
    64*srate:65*srate, 65.5*srate:66*srate, ...
    66.5*srate:129.2*srate, 129.6*srate:130.5*srate, ...
    131*srate:150.5*srate, 151*srate:209.5*srate,...
    210*srate:220.5*srate, 221.5*srate:238.5*srate,...
    239*srate:309*srate, 310*srate:321.3*srate,...
    321.5*srate:478*srate, 478.5*srate:544.5*srate, ...
    545*srate:549.9*srate, 550*srate:end]);
G0302d = AnimalGato0302d2([1:65*srate, 65.5*srate:99.8*srate, ...
    100*srate:127*srate, 128*srate:371.9*srate,...
    372.1*srate:401.3*srate, 401.5*srate:463.2*srate,...
    463.4*srate:494.8*srate, 495*srate:end]);
G0502d = AnimalGato0502d2([1:135.1*srate, 135.2*srate:142.9*srate, ...
    143*srate:327.9*srate, 328*srate:426.4*srate, ...
    426.8*srate:431.4*srate, 431.6*srate:451.5*srate, ...
    451.6*srate:563.4*srate, 563.5*srate:end]);
G0702d = AnimalGato0702d2([1:9.05*srate, 9.15*srate:361.88*srate, ...
    361.92*srate:end]);
G1102d = AnimalGato1102d2([ 1:80.3*srate, 80.5*srate:113.2*srate, ...
    113.5*srate:160.1*srate, 160.3*srate:276.5*srate, ...
    286*srate:339.5*srate, 339.6*srate:378*srate, ...
    383*srate:394.2*srate, 394.6*srate:524*srate, ...
    528*srate:end]);

XAMBGd = {G0201d;G0301d;G0501d;G0701d;G1101d};
XEXPGd = {G0202d;G0302d;G0502d;G0702d;G1102d};

C0101d = AnimalControle0101d2([20*srate:274*srate, 280.4*srate:520*srate, ...
    527*srate:end]);
C0401d = AnimalControle0401d2([1:27*srate, 31*srate:165*srate, 184*srate:234*srate, 244*srate:483*srate, ...
    487*srate:end]);
C0601d = AnimalControle0601d2([1:end]);
C0901d = AnimalControle0901d2([1:514*srate,520*srate:end]);
C1001d = AnimalControle1001d2([1:end]);

C0102d = AnimalControle0102d2([20*srate:end]);
C0402d = AnimalControle0402d2([1:50*srate, 55*srate:285*srate, 289*srate:end]);
C0602d = AnimalControle0602d2([1:end]);
C0902d = AnimalControle0902d2([12*srate:543*srate, 547*srate:590*srate]);
C1002d = AnimalControle1002d2([1:end]);

XAMBCd = {C0101d;C0401d;C0601d;C0901d;C1001d};
XEXPCd = {C0102d;C0402d;C0602d;C0902d;C1002d};

% clear Animal* G*
%% COMPARAÇÃO DIA 01 - TREINO
clear PsdAMB PsdEXP SpectAMB SpectEXP T1 T2 TAMB TEXP
for a=1:5   %animals
    clear data data2 PSD P
    srate      = 500;
    idx        = 1/srate;
    
    data       = XAMBC{a}; %Habituation
    data2      = XEXPC{a}; %Exposition
    
    window     = 5*srate;
    overlap    = window*0.25;
    nfft       = 2^13;
    %Habituation
    [PSDAMB, F]   = pwelch(data(1:300*srate),window,overlap,nfft,srate);
    PsdAMB(a,:) = PSDAMB;
    [S F TAMB{a} PAMB{a}]  = spectrogram(data(1:300*srate),window,overlap,nfft,srate);
    [ACGAMB(a,:), lags] = xcorr(data,data,srate,'coeff');
    %Exposition
    [PSDEXP, F]   = pwelch(data2(1:300*srate),window,overlap,nfft,srate);
    PsdEXP(a,:) = PSDEXP;
    [S F TEXP{a} PEXP{a}]  = spectrogram(data2(1:300*srate),window,overlap,nfft,srate);
    [ACGEXP(a,:), lags] = xcorr(data2,data2,srate,'coeff');
end

%
clear Spect*
TAMBmax = size(TAMB{1},2); % GRUPO CONTROLE DIA 01 E 02- {2} / GRUPO EXPOSIÇÃO DIA 01- {1}, DIA 02 {5}
TEXPmax = size(TEXP{5},2); % GRUPO CONTROLE DIA 01 E 02- {4} / GRUPO EXPOSIÇÃO DIA 01- {5}, DIA 02 {5}
for a = 1:5
    SpectAMB(:,:,a) = PAMB{a}(:,1:TAMBmax);
    SpectEXP(:,:,a) = PEXP{a}(:,1:TEXPmax);
end

%%  TESTE T NOS DADOS NORMALIZADOS PELA MEDIA - TEMPO TOTAL
clc
Bands=[2,6; 7,12; 15,25; 36,40];%definir bandas de interesse
clear h p IC effect s
for b = 1:4
    clear idx Norm
    idx          = find(F>Bands(b,1) & F<Bands(b,2));
    Norm         = mean([mean(PsdAMB(:,idx),2), mean(PsdEXP(:,idx),2)],2);
    %     NPAMB{b} = PsdAMB(:,idx)
    %     NPEXP{b} = PsdEXP(:,idx)
    %     NP{b} = mean(mean(NPAMB{1},2),1)/mean(mean(NPAMB{b},2),1);
    %     NP2{b} = mean(mean(NPEXP{1},2),1)/mean(mean(NPEXP{b},2),1);
    NormAMB(:,b) = mean(PsdAMB(:,idx),2)./Norm;
    NormEXP(:,b) = mean(PsdEXP(:,idx),2)./Norm;
    [h(b) p(b) ci stats]  = ttest(NormEXP(:,b), NormAMB(:,b),'Tail','both');
    
    IC(:,:,b) = ci;
    s{b} = stats;
    effect(:,b) = computeCohen_d(NormEXP(:,b), NormAMB(:,b), 'paired');
end
disp([p(1) p(2) p(3) p(4)])
%% T-TEST - TREINO

FIG1=figure(1);clf
subplot(4,4,[1 4])
LFP_AMB = C0601;
LFP_EXP = C0602;
dt = 1/srate;
time_vector = dt:dt:length(LFP_AMB)/srate;
plot(time_vector,LFP_AMB,'k','linewidth',1)
hold on
plot(time_vector,LFP_EXP-1,'r','linewidth',1)
ylim([-1.8 0.5])
box off
xlabel 'Time (s)'
x1=38;
xlim([x1 x1+10])
set(gca,'ytick',[-1 0],'yticklabels',{'Exposition','Habituation'})

% %
delta_AMB = eegfilt(LFP_AMB,srate,2,6);
delta_EXP = eegfilt(LFP_EXP,srate,2,6);

theta_AMB = eegfilt(LFP_AMB,srate,7,12);
theta_EXP = eegfilt(LFP_EXP,srate,7,12);

beta_AMB = eegfilt(LFP_AMB,srate,15,25);
beta_EXP = eegfilt(LFP_EXP,srate,15,25);

gama_AMB = eegfilt(LFP_AMB,srate,36,40);
gama_EXP = eegfilt(LFP_EXP,srate,36,40);

subplot(4,4,[1 4])
plot(time_vector,delta_AMB-0.2,'k','linewidth',2)
plot(time_vector,delta_EXP-1.2,'r','linewidth',2)
plot(time_vector,theta_AMB-0.35,'k','linewidth',2)
plot(time_vector,theta_EXP-1.35,'r','linewidth',2)
plot(time_vector,beta_AMB-0.5,'k','linewidth',2)
plot(time_vector,beta_EXP-1.5,'r','linewidth',2)
plot(time_vector,gama_AMB-0.5,'k','linewidth',2)
plot(time_vector,gama_EXP-1.5,'r','linewidth',2)

% %
%Spectrograms
subplot(4,4,[5 6])
imagesc(TAMB{5},F,mean(SpectAMB,3))
axis xy
ylim([20 40])
caxis([0 0.00005])
title 'Habituation'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'


subplot(4,4,[7 8])
imagesc(TEXP{5},F,mean(SpectEXP,3))
axis xy
ylim([20 40])
caxis([0 0.00005])
title 'Exposition'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'

set(gcf,'color','w')

subplot(4,4,[9 10])
imagesc(TAMB{5},F,mean(SpectAMB,3))
axis xy
ylim([0 20])
caxis([0 0.0005])
title 'Habituation'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'


subplot(4,4,[11 12])
imagesc(TEXP{5},F,mean(SpectEXP,3))
axis xy
ylim([0 20])
caxis([0 0.0005])
title 'Exposition'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'

set(gcf,'color','w')

subplot(4,4,[13 14])
plot(F,log(mean(PsdAMB)),'k','linewidth',2)
hold on
plot(F,log(mean(PsdEXP)),'r','linewidth',2)
xlim([0 40])
plot(F,log(mean(PsdAMB)+std(PsdAMB)/sqrt(5)),'k--')
plot(F,log(mean(PsdAMB)-std(PsdAMB)/sqrt(5)),'k--')
plot(F,log(mean(PsdEXP)+std(PsdEXP)/sqrt(5)),'r--')
plot(F,log(mean(PsdEXP)-std(PsdEXP)/sqrt(5)),'r--')
ylim([-14 -7])
title 'Habituation x Exposition'
xlabel 'Frequency (Hz)'
ylabel 'Power (log)'
Legend={'Hab','Exp'};
legend(Legend,'box','off')
box off

subplot(4,4,[15 16])
bar([1:1:4],[mean(NormAMB)],0.2,'k')
hold on
bar([1.25:1:4.25],[mean(NormEXP)],0.2,'r')
errorbar([1:1:4],[mean(NormAMB)],[std(NormAMB)/sqrt(5)],'.k')
errorbar([1.25:1:4.25],[mean(NormEXP)],[std(NormEXP)/sqrt(5)],'.k')
xlabel 'Band Frequency (Hz)'
ylabel 'Norm Power'
ylim([0.6 1.4])
xlim([0.75 4.75])
box off
set(gca,'xtick',[1.125:1:4.125],'xticklabels',{'2-6','7-12','15-25','36-40'})
% set(gca,'xtick',[1:5],'xticklabels',{'Delta','Theta','Beta','Gamma'})
for b=1:4
    if p(b)<0.05
        text(b,1.3,'*')
    end
end

%% ACUMULATIVO - TREINO
clc
clear Psd* Spect* T* ACG*
c = 1;
for time = 60:60:300             % 60:60:300
    for a=1:5                    % animals
        clear data data2 PSD P
        srate      = 500;
        idx        = 1/srate;
        
        data       = XAMBG{a};  % Habituation
        data2      = XEXPG{a};  % Exposition
        
        window     = 5*srate;
        overlap    = window*0.25;
        nfft       = 2^13;
        
        %Habituation
        [PSD1, F]   = pwelch(data(1:300*srate),window,overlap,nfft,srate);
        Psd1(a,:)   = PSD1;
        [S F T1{a} P1d{a}]  = spectrogram(data(1:300*srate),window,overlap,nfft,srate);
        [ACGEXP(a,:), lags]    = xcorr(data,data,srate,'coeff');
        
        %Exposition
        [PSD2, F]   = pwelch(data2(1:time*srate),window,overlap,nfft,srate);
        Psd2(a,:)   = PSD2;
        Psd2D{c}     = Psd2;
        [S F T2{a} P2d{a}]  = spectrogram(data2(1:time*srate),window,overlap,nfft,srate);
        [ACGEXPd(a,:), lags]     = xcorr(data2,data2,srate,'coeff');
        
        
    end
    
    ACGEXPP{c} = ACGEXP;
    c = c+1;
    
end

clear Spect*
T1max = size(T1{5},2); % GRUPO CONTROLE DIA 01 E 02- {2} / GRUPO EXPOSIÇÃO DIA 01- {1}, DIA 02 {5}
T2max = size(T2{5},2); % GRUPO CONTROLE DIA 01 E 02- {4} / GRUPO EXPOSIÇÃO DIA 01- {5}, DIA 02 {5}
for a = 1:5
    Spect1(:,:,a) = P1d{a}(:,1:T1max);
    Spect2(:,:,a) = P2d{a}(:,1:T2max);
    
end

%%  TESTE T - TREINO - TEMPO ACUMULATIVO
clc
Bands=[2,6; 7,12; 15,25; 36,40];%definir bandas de interesse
clear h p IC effect EP*
c = 1;
clc
for time = 1:5  %acumulative time
    
    for b = 1:4 %bands    
        clear idx Norm
        idx          = find(F>Bands(b,1) & F<Bands(b,2));
        Norm         = mean([mean(Psd1(:,idx),2), mean(Psd2D{time}(:,idx),2)],2);
        Norm1(:,b) = mean(Psd1(:,idx),2)./Norm;
        Norm1D{time} = Norm1;
%         EPNorm1(time,b)  = std(Norm1D{time}(:,b))/sqrt(5)
        Norm2(:,b)= mean(Psd2D{time}(:,idx),2)./Norm;
        Norm2D{time} = Norm2;
%         EPNorm2D(time,b) = std(Norm2D{time}(:,b))/sqrt(5)
        [h(b) p{time}(b) ci stats]  = ttest(Norm1D{time}(:,b), Norm2D{time}(:,b),'Tail','both')
        IC{time}(:,:,b) = ci;
        s{time,b} = stats;
        effect{time}(:,b) = computeCohen_d(Norm1(:,b), Norm2D{time}(:,b), 'paired');
        
    end
    c = c+1;
    
end



%% COMPARAÇÃO DIA 02 - TESTE
clear PsdAMB PsdEXP SpectAMB SpectEXP T1 T2 TAMB TEXP
for a=1:5                   %animals
    clear data data2 PSD P
    srate      = 500;
    idx        = 1/srate;
    
    data       = XAMBCd{a};
    data2      = XEXPCd{a};
    
    window     = 5*srate;
    overlap    = window*0.25;
    nfft       = 2^13;
    [PSDAMB, F]   = pwelch(data(1:300*srate),window,overlap,nfft,srate);
    PsdAMB(a,:) = PSDAMB;%./sum(PSDAMB);
    [S F TAMB{a} PAMB{a}]  = spectrogram(data(1:300*srate),window,overlap,nfft,srate);
    [ACGAMB(a,:), lags] = xcorr(data,data,srate,'coeff');
    
    nfft       = 2^13;
    [PSDEXP, F]   = pwelch(data2(1:300*srate),window,overlap,nfft,srate);
    PsdEXP(a,:) = PSDEXP;%./sum(PSDEXP);
    [S F TEXP{a} PEXP{a}]  = spectrogram(data2(1:300*srate),window,overlap,nfft,srate);
    [ACGEXP(a,:), lags] = xcorr(data2,data2,srate,'coeff');
end

%
clear Spect*
TAMBmax = size(TAMB{1},2); % GRUPO CONTROLE DIA 01 E 02- {2} / GRUPO EXPOSIÇÃO DIA 01- {1}, DIA 02 {5}
TEXPmax = size(TEXP{5},2); % GRUPO CONTROLE DIA 01 E 02- {4} / GRUPO EXPOSIÇÃO DIA 01- {5}, DIA 02 {5}
for a = 1:5
    SpectAMB(:,:,a) = PAMB{a}(:,1:TAMBmax);
    SpectEXP(:,:,a) = PEXP{a}(:,1:TEXPmax);
end
%% T-TEST - TEST
%  TESTE T NOS DADOS NORMALIZADOS PELA MEDIA - TEMPO TOTAL
clc
Bands=[2,6; 7,12; 15,25; 36,40];%definir bandas de interesse
clear h p IC effect s Norm*
for b = 1:4 %bands
    clear idx Norm
    idx          = find(F>Bands(b,1) & F<Bands(b,2));
    Norm         = mean([mean(PsdAMB(:,idx),2), mean(PsdEXP(:,idx),2)],2);
    %     NPAMB{b} = PsdAMB(:,idx)
    %     NPEXP{b} = PsdEXP(:,idx)
    %     NP{b} = mean(mean(NPAMB{1},2),1)/mean(mean(NPAMB{b},2),1);
    %     NP2{b} = mean(mean(NPEXP{1},2),1)/mean(mean(NPEXP{b},2),1);
    NormAMB(:,b) = mean(PsdAMB(:,idx),2)./Norm;
    NormEXP(:,b) = mean(PsdEXP(:,idx),2)./Norm;
    [h(b) p(b) ci stats]  = ttest(NormAMB(:,b), NormEXP(:,b),'Tail','both')
    IC(:,:,b) = ci;
    s{b} = stats;
    effect(:,b) = computeCohen_d(NormAMB(:,b), NormEXP(:,b), 'paired');
end

%% T-TEST - TEST

FIG1=figure(1);clf
subplot(4,4,[1 4])
LFP_AMB = C0601d;
LFP_EXP = C0602d;
dt = 1/srate;
time_vector = dt:dt:length(LFP_AMB)/srate;
plot(time_vector,LFP_AMB,'k','linewidth',1)
hold on
plot(time_vector,LFP_EXP-1,'r','linewidth',1)
ylim([-1.8 0.5])
box off
xlabel 'Time (s)'
x1=38;
xlim([x1 x1+10])
set(gca,'ytick',[-1 0],'yticklabels',{'Exposition','Habituation'})

% %
delta_AMB = eegfilt(LFP_AMB,srate,2,6);
delta_EXP = eegfilt(LFP_EXP,srate,2,6);

theta_AMB = eegfilt(LFP_AMB,srate,7,12);
theta_EXP = eegfilt(LFP_EXP,srate,7,12);

beta_AMB = eegfilt(LFP_AMB,srate,15,25);
beta_EXP = eegfilt(LFP_EXP,srate,15,25);

gama_AMB = eegfilt(LFP_AMB,srate,36,40);
gama_EXP = eegfilt(LFP_EXP,srate,36,40);

subplot(4,4,[1 4])
plot(time_vector,delta_AMB-0.2,'k','linewidth',2)
plot(time_vector,delta_EXP-1.2,'r','linewidth',2)
plot(time_vector,theta_AMB-0.35,'k','linewidth',2)
plot(time_vector,theta_EXP-1.35,'r','linewidth',2)
plot(time_vector,beta_AMB-0.5,'k','linewidth',2)
plot(time_vector,beta_EXP-1.5,'r','linewidth',2)
plot(time_vector,gama_AMB-0.5,'k','linewidth',2)
plot(time_vector,gama_EXP-1.5,'r','linewidth',2)

% %
%Spectrograms
subplot(4,4,[5 6])
imagesc(TAMB{5},F,mean(SpectAMB,3))
axis xy
ylim([20 40])
caxis([0 0.00005])
title 'Habituation'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'


subplot(4,4,[7 8])
imagesc(TEXP{5},F,mean(SpectEXP,3))
axis xy
ylim([20 40])
caxis([0 0.00005])
title 'Exposition'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'

set(gcf,'color','w')

subplot(4,4,[9 10])
imagesc(TAMB{5},F,mean(SpectAMB,3))
axis xy
ylim([0 20])
caxis([0 0.0005])
title 'Habituation'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'


subplot(4,4,[11 12])
imagesc(TEXP{5},F,mean(SpectEXP,3))
axis xy
ylim([0 20])
caxis([0 0.0005])
title 'Exposition'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'

set(gcf,'color','w')

subplot(4,4,[13 14])
plot(F,log(mean(PsdAMB)),'k','linewidth',2)
hold on
plot(F,log(mean(PsdEXP)),'r','linewidth',2)
xlim([0 40])
plot(F,log(mean(PsdAMB)+std(PsdAMB)/sqrt(5)),'k--')
plot(F,log(mean(PsdAMB)-std(PsdAMB)/sqrt(5)),'k--')
plot(F,log(mean(PsdEXP)+std(PsdEXP)/sqrt(5)),'r--')
plot(F,log(mean(PsdEXP)-std(PsdEXP)/sqrt(5)),'r--')
ylim([-14 -7])
title 'Habituation x Exposition'
xlabel 'Frequency (Hz)'
ylabel 'Power (log)'
Legend={'Hab','Exp'};
legend(Legend,'box','off')
box off

subplot(4,4,[15 16])
bar([1:1:4],[mean(NormAMB)],0.2,'k')
hold on
bar([1.25:1:4.25],[mean(NormEXP)],0.2,'r')
errorbar([1:1:4],[mean(NormAMB)],[std(NormAMB)/sqrt(5)],'.k')
errorbar([1.25:1:4.25],[mean(NormEXP)],[std(NormEXP)/sqrt(5)],'.k')
xlabel 'Band Frequency (Hz)'
ylabel 'Norm Power'
ylim([0.6 1.4])
xlim([0.75 4.75])
box off
set(gca,'xtick',[1.125:1:4.125],'xticklabels',{'2-6','7-12','15-25','36-40'})
% set(gca,'xtick',[1:5],'xticklabels',{'Delta','Theta','Beta','Gamma'})
for b=1:4
    if p(b)<0.05
        text(b,1.3,'*')
    end
end
%% ACUMULATIVO - TEST
clc
clear Spect* T* ACG* P*
c = 1;
for time = 60:60:300             % 60:60:300
    for a=1:5                    % animals
        clear data data2 PSD* 
        srate      = 500;
        idx        = 1/srate;
        
        data       = XAMBGd{a}; %habituation
        data2      = XEXPGd{a}; %exposition
        
        window     = 5*srate;
        overlap    = window*0.25;
        nfft       = 2^13;
        
        %Habituation in total time
        [PSD1, F]   = pwelch(data(1:300*srate),window,overlap,nfft,srate);
        Psd1(a,:)   = PSD1;
        [S F T1{a} P1d{a}]  = spectrogram(data(1:300*srate),window,overlap,nfft,srate);
        [ACGEXP(a,:), lags]    = xcorr(data,data,srate,'coeff');
        %Exposition in acumulative time
        [PSD2, F]   = pwelch(data2(1:time*srate),window,overlap,nfft,srate);
        Psd2(a,:)   = PSD2;
        Psd2d{c}     = Psd2;
        [S F T2{c} P2d{c}(:,:,a)]  = spectrogram(data2(1:time*srate),window,overlap,nfft,srate);
        [ACGEXPd(a,:), lags]     = xcorr(data2,data2,srate,'coeff');
        
        
    end
    Psd2D(:,:,c)     = Psd2;

    
    ACGEXPP{c} = ACGEXP;
    c = c+1;
    
end

clear Spect*

for c = 1:5 %time
    T1max = size(T1{c},2); % GRUPO CONTROLE DIA 01 E 02- {2} / GRUPO EXPOSIÇÃO DIA 01- {1}, DIA 02 {5}
    T2max = size(T2{c},2); % GRUPO CONTROLE DIA 01 E 02- {4} / GRUPO EXPOSIÇÃO DIA 01- {5}, DIA 02 {5}
    for a = 1:5 %animal
        Spect1(:,:,a) = P1d{a}(:,1:T1max);
        Spect2{c}(:,:,a) = P2d{c}(:,1:T2max,a);
        
    end
end
%%  TESTE T - TESTE - TEMPO ACUMULATIVO
clc
Bands=[2,6; 7,12; 15,25; 36,40];%definir bandas de interesse
clear h p IC effect EP*
c = 1;
clc
for time = 1:5 %acumulative time -> 0 - 60' + 60'
    
    for b = 1:4 %bands
        clear idx Norm
        idx          = find(F>Bands(b,1) & F<Bands(b,2));
        Norm         = mean([mean(Psd1(:,idx),2), mean(Psd2d{time}(:,idx),2)],2);
        Norm1(:,b) = mean(Psd1(:,idx),2)./Norm;
        Norm1D{time} = Norm1;
        %         EPNorm1(time,b)  = std(Norm1D{time}(:,b))/sqrt(5)
        Norm2(:,b)= mean(Psd2d{time}(:,idx),2)./Norm;
        Norm2D{time} = Norm2;
        %         EPNorm2D(time,b) = std(Norm2D{time}(:,b))/sqrt(5)
        [h{time}(b) p{time}(b) ci stats]  = ttest(Norm1D{time}(:,b), Norm2D{time}(:,b), 'Alpha',0.01,'Tail','both')
        IC{time}(:,:,b) = ci;
        s{time,b} = stats;
        effect{time}(:,b) = computeCohen_d(Norm1(:,b), Norm2D{time}(:,b), 'paired');
        
    end
    c = c+1;
    
end

%%
% Acumulative Spectrograms
time = 0;
for c = 1:1:5
subplot(2,5,c)
imagesc(T1{c},F,mean(Spect1,3))
axis xy
ylim([0 20])
xlim([0 (60 + time)])
caxis([0 0.0005])
title 'Habituation'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'


subplot(2,5,c+5)
imagesc(T2{c},F,mean(Spect2{c},3))
axis xy
ylim([0 20])
caxis([0 0.0005])
title 'Exposition'
xlabel 'Time (s)'
ylabel 'Frequency (Hz)'

set(gcf,'color','w')
time = time + 60
end
%% Acumulative Bars
for c = 1:1:5
subplot(1,5,c)
bar([1:1:4],[mean(Norm1D{c})],0.2,'k')
hold on
bar([1.25:1:4.25],[mean(Norm2D{c})],0.2,'r')
errorbar([1:1:4],[mean(Norm1D{c})],[std(Norm1D{c})/sqrt(5)],'.k')
errorbar([1.25:1:4.25],[mean(Norm2D{c})],[std(Norm2D{c})/sqrt(5)],'.k')
xlabel 'Band Frequency (Hz)'
ylabel 'Norm Power'
ylim([0.5 1.5])
xlim([0.75 4.75])
box off
set(gca,'xtick',[1.125:1:4.125],'xticklabels',{'2-6','7-12','15-25','36-40'})
%set(gca,'xtick',[1:5],'xticklabels',{'Delta','Theta','Beta','Gamma'})
for b=1:4
    if p{c}(b)<0.05
        text(b,1.3,'*')
    end
end
end