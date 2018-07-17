% Propagation process equivalence in networks
% Finds the equivalent networks (see [1]) and compares
% for a two-tier system Vs one-tier system cellular network
% Calculates A and beta values for COST231-Hata model
% Heights in metres, distances in kilometres
% Author: H.P. Keeler, Inria Paris/ENS, 2013
%

% [1] B. Błaszczyszyn and H. Keeler
% 'Equivalence and comparison of heterogeneous cellular networks',
% accepted at WDN-CN2013, 2013,
% http://arxiv.org/abs/1306.0772
%


%Propagation loss check for two-tier system Vs one-tier system
%Calculates A and beta values for COST231-Hata model
%heights in metres, distances in kilometres

clear all; clc;
close all;
lambda1=1.8; lambda2=2.2;%base station density for two-tier system

diskRadius=15;

sigmaDb1=5;  %fading log-normal std parameter franges from 6 to 12


%path-loss function parameters - calculate using COST231-Hata model
%inputs
hata_hb1=20;    hata_hb2=100; %antenna heights (m)
%hata_hb3=(hata_hb1+hata_hb2)/2;       %regular average of antenna heights (m)
hata_hb3=(hata_hb1*lambda1+hata_hb2*lambda2)/(lambda1+lambda2);      %spatial average of antenna heights (m)

hata_fc1=1800;  hata_fc2=hata_fc1; hata_fc3=hata_fc1;   %carrier frequency (MHz)
hata_hm1=1;    hata_hm2=hata_hm1;  %user height
hata_hm3=hata_hm1;
largecity=1;
%model parameter calculations
if (largecity==1)
    hata_ahm1=3.2*(log10(11.75*hata_hm1))^2-4.97;
    hata_ahm2=3.2*(log10(11.75*hata_hm2))^2-4.97;
    hata_ahm3=3.2*(log10(11.75*hata_hm3))^2-4.97;
    hata_C1=3;      hata_C2=hata_C1;          %3 for metropolitan areas
    hata_C3=hata_C1;
else
    hata_ahm1=(1.1*log10(hata_fc1)-0.7)*hata_hm1-(1.56*log10(hata_fc1) -0.8);
    hata_ahm2=(1.1*log10(hata_fc2)-0.7)*hata_hm2-(1.56*log10(hata_fc2) -0.8);
    hata_ahm3=(1.1*log10(hata_fc3)-0.7)*hata_hm3-(1.56*log10(hata_fc3) -0.8);
    hata_C1=0;      hata_C2=hata_C1; hata_C3=hata_C1;         %3 for metropolitan areas
    
end

hata_A1=46.3+33.9*log10(hata_fc1)-13.82*log10(hata_hb1)-hata_ahm1;
hata_A2=46.3+33.9*log10(hata_fc2)-13.82*log10(hata_hb2)-hata_ahm2;
hata_B1=44.9-6.55*log10(hata_hb1);
hata_B2=44.9-6.55*log10(hata_hb2);

hata_A3=46.3+33.9*log10(hata_fc3)-13.82*log10(hata_hb3)-hata_ahm3;
hata_B3=44.9-6.55*log10(hata_hb3);

betaConst1=hata_B1/10
A1=10^((hata_A1+hata_C1)/10)
betaConst2=hata_B2/10
A2=10^((hata_A2+hata_C2)/10)
betaConst3=hata_B3/10
A3=10^((hata_A3+hata_C3)/10)

betaPrime=betaConst3;
%betaConst1=3.1; betaConst2=3.8; %path-loss exponent
%A1=1;A2=3;A3=2;

%switch for random fading
fadingSwitch=1; %0 for constant fading; 1 for random fading
%choose distribution type %1 to 3 (0= no shadowing)
distType=1;
ES1=1; ES2=1;
switch(distType)
    case 0 %no shadowing -dist 0
        ESTwoBeta1=1;
        ESTwoBeta2=1;
    case 1  %log normal parameters - dist 1
        %v1=1/10;                
        %sigma1=sqrt(log(v1/ES1^2+1));        
        sigma1=sigmaDb1/10*log(10);        
        v1=ES1^2*(exp(sigma1^2)-1);        
        mu1=log(ES1^2/sqrt(v1+ES1^2));        
        sigmaDb1=sigma1*10/log(10);
        v2=v1; mu2=mu1;sigma2=sigma1;        
        
        ESTwoBeta1=exp(2*(sigma1^2+betaConst1*mu1)/betaConst1^2); %log-normal
        ESTwoBeta2=exp(2*(sigma2^2+betaConst2*mu2)/betaConst2^2); %log-normal
    case 2 %weibull paramter - dist 2
        kW1=1;        kW2=1; %kw=1 corresponds to exponential
        lambdaW1=ES1/gamma(1/kW1+1); lambdaW2=ES2/gamma(1/kW2+1);
        ESTwoBeta1=lambdaW1^(2/betaConst1)*gamma((2/betaConst1)/kW1+1); %Weibull
        ESTwoBeta2=lambdaW2^(2/betaConst1)*gamma((2/betaConst2)/kW2+1); %Weibull
    case 3 %Nakagami paramters - dist 3
        omega1=3; omega2=3;
        mN1=5; mN2=5;
        ESTwoBeta1=gamma((1/betaConst1)+mN1)/gamma(mN1)*(omega1/mN1)^(1/betaConst1); %Nakagami
        ESTwoBeta2=gamma((1/betaConst2)+mN2)/gamma(mN2)*(omega2/mN2)^(1/betaConst2); %Nakagami
end

%model constant - incorporates model parameters
a1=lambda1*pi*ESTwoBeta1/A1^(2/betaConst1);
a2=lambda2*pi*ESTwoBeta2/A2^(2/betaConst2);

a1NoFade=a1/ESTwoBeta1;
a2NoFade=a2/ESTwoBeta2;

%calculate phi
%single tier - distinct betas and A
rmin=.5; rmax=15;
rValues=linspace(rmin,rmax,20);
phi2tier=(a1*rValues.^(2*betaPrime/betaConst1-2)*(betaPrime/betaConst1) ...
    +a2*rValues.^(2*betaPrime/betaConst2-2)*(betaPrime/betaConst2))/pi;
phi2tier_rescaled=phi2tier*A3^(2/betaConst3);

phi2tierNoFade=(a1NoFade*rValues.^(2*betaPrime/betaConst1-2)*(betaPrime/betaConst1) ...
    +a2NoFade*rValues.^(2*betaPrime/betaConst2-2)*(betaPrime/betaConst2))/pi;
phi2tier_rescaledNoFade=phi2tierNoFade*A3^(2/betaConst3);


%single tier - equal beta and A
A1=A3; A2=A3;
betaConst1=betaConst3;
betaConst2=betaConst3;
%betaPrime=(betaConst1+betaConst2)/2; %average of two %other besta Values
a1=lambda1*pi*ESTwoBeta1/A1^(2/betaConst1);
a2=lambda2*pi*ESTwoBeta2/A2^(2/betaConst2);
a1NoFade=a1/ESTwoBeta1;
a2NoFade=a2/ESTwoBeta2;

phi1tier=(a1*rValues.^(2*betaPrime/betaConst1-2)*(betaPrime/betaConst1) ...
    +a2*rValues.^(2*betaPrime/betaConst2-2)*(betaPrime/betaConst2))/pi;
phi1tier_rescaled=phi1tier*A1^(2/betaConst1);

phi1tierNoFade=(a1NoFade*rValues.^(2*betaPrime/betaConst1-2)*(betaPrime/betaConst1) ...
    +a2NoFade*rValues.^(2*betaPrime/betaConst2-2)*(betaPrime/betaConst2))/pi;
phi1tier_rescaledNoFade=phi1tierNoFade*A1^(2/betaConst1);


%plotting
figure;
set(gcf,'DefaultLineLineWidth',2);
plot(rValues,phi2tier_rescaled,rValues,phi1tier_rescaled,'--',...
    rValues,phi2tier_rescaledNoFade,'x',rValues,phi1tier_rescaledNoFade,'o');
xlabel('r (km)','fontsize',16); ylabel('\phi(r)E(A^{2/\beta})','fontsize',16); grid;
legend('Two-tier with fading','One-tier with fading','Two-tier without fading','One-tier without fading');
ymax=max(max(phi1tier_rescaled,phi2tier_rescaled));
axis([rmin rmax 0 1.5*ymax]);
set(gca, 'FontSize', 14);

