% Calculates A and beta values for COST231-Hata model
% heights in metres, distances in kilometres
%
% This work was used for [1,3] to show the equivalence in different types of 
% cellular network models based on Poisson point processes
%
% Author: H.P. Keeler, Inria Paris/ENS, 2013
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
%  shadowing ', presented to ISIT, 2013
% http://arxiv.org/abs/1301.6491
% [2] B. Błaszczyszyn and H. Keeler
% 'Equivalence and comparison of heterogeneous cellular networks',
%  presented at WDN-CN2013, 2013,
% http://arxiv.org/abs/1306.0772
% [3] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% random heterogeneous cellular networks' , submitted to TOIT, 2013013


function [betaValues, KValues]=funHataParameters(height1,height2,lambda1,lambda2)



%path-loss function parameters - calculate using COST231-Hata model
%inputs
hata_hb1=height1; hata_hb2=height2; %antenna heights (m)
%hata_hb1=20;    hata_hb2=100; 
%hata_hb3=(hata_hb1+hata_hb2)/2;       %antenna heights (m)
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

betaConst1=hata_B1/10;
A1=10^((hata_A1+hata_C1)/10);
%caclulates K - propagation constant used in [1] and [3]
K1=A1^(1/betaConst1);
betaConst2=hata_B2/10;
A2=10^((hata_A2+hata_C2)/10);
K2=A2^(1/betaConst2);
betaConst3=hata_B3/10;
A3=10^((hata_A3+hata_C3)/10);
K3=A3^(1/betaConst3);


betaValues=[betaConst1,betaConst2,betaConst3];
KValues=[K1,K2,K3];
