%Cellular Beam
clc;

Smm  = 300;  % Centroid Distance Between Holes
Domm = 210; % Diameter of Holes
emm  = Smm - Domm;  % Boundary Distance Between Holes
dmm  = 270;
twmm = 9.5;
bfmm = 210;
tfmm = 9.5;

FyNmm2 = 355;

S  = Smm/25.4;  % Centroid Distance Between Holes
Do = Domm/25.4; % Diameter of Holes
e  = emm/25.4;  % Boundary Distance Between Holes
d  = dmm/25.4;
tw = twmm/25.4;
bf = bfmm/25.4;
tf = tfmm/25.4;

YukKNm = 5.88;
YukKipft = YukKNm*0.06854;

SapMesh = 20;

%Fy = FyNmm2/6.895;

Fy = 50;
Fu = 65;


loss = Do/2 - sqrt((Do/2)^2 - ((S-Do)/2)^2);
dg = d + Do/2 - loss;


E = 29000;
G = 11200;


Check1 = S/Do;
Check2 = dg/Do;
Check3 = e;



[KesmeMoment,Makaslar] = xlsread('MomentKesme.xlsx','A1:D10000');
KesmeMoment(:) = abs(KesmeMoment(:));

Bsp = e;
Artan = S*mod((39.3700787*max(KesmeMoment(:,1))-2*e)/S,1);

if e > Artan/2
    Bsp = e + Artan/2;
elseif e < Artan/2
    Bsp = Artan/2;
end

NumOfHoles = round(39.3700787*max(KesmeMoment(:,1))/S-mod(39.3700787*max(KesmeMoment(:,1))/S,1));

HoleLoc(1,1) = Bsp;

%Net Kesitlerin Lokasyonlarýnýn Bulunmasý
for i=2:NumOfHoles 
    HoleLoc(i,1) = HoleLoc(i-1,1)+S;
end

%Brüt Kesitlerin Lokasyonlarýnýn Bulunmasý

GrossLoc(1,1) = Bsp/2;
GrossLoc((length(HoleLoc(:,1))+1),1) = (Bsp/2 + HoleLoc(length(HoleLoc(:,1)),1));

for i=2:NumOfHoles
    GrossLoc(i,1) = (HoleLoc(i-1,1)+HoleLoc(i,1))/2;
end

GrossLocm = GrossLoc*0.0254;
HoleLocm = HoleLoc*0.0254;

%for i=1:SapMesh-1
%Mati(:,i) = find(KesmeMoment(:,1)==KesmeMoment(i,1));

%HoleLocmU(Mati(:,i),1) = ((KesmeMoment(i+1,2) - KesmeMoment(i,2))/(KesmeMoment(i+1,1) - KesmeMoment(i,1)))*(HoleLocm(i,1)-KesmeMoment(i,1)) + KesmeMoment(i,2);
%end


if 1.08 < S/Do && S/Do < 1.5 
    fprintf('Check Geometric Limits 1.08 < S/Do = %.2f < 1.5 is OKAY \n',Check1)
else
    fprintf('Check Geometric Limits 1.08 > S/Do = %.2f < 1.5 is NOT PROVIDED \n',Check1)
end


if 1.25 < dg/Do && dg/Do < 1.75
    fprintf('Check Geometric Limits 1.25 < dg/Do = %.2f < 1.75 is OKAY \n',Check2)
else
    fprintf('Check Geometric Limits 1.25 > dg/Do = %.2f < 1.75 is NOT PROVIDED \n',Check2)
end


if e > 24.4
    fprintf('Check Geometric Limits e = %.2f > 24.4 is OKAY \n',Check3)
else
    fprintf('Check Geometric Limits e = %.2f < 24.4 is NOT PROVIDED \n',Check3)
end

fprintf('\n')


dtnet = (dg- Do)/2;
y = sqrt((0.5*Do)^2 - (0.225*Do)^2);
dtcrit = (Do/2) - y + dtnet;

AteeNet = dtnet*tw + tf*bf;
AteeCrit = dtcrit*tw + tf*bf;
ANet = 2*AteeNet;
ACrit = 2*AteeCrit;
TeeGovdeNet = dtnet - tf;
TeeGovdeCrit = dtcrit - tf;

ycnet = (bf*tf*(TeeGovdeNet+tf/2) + tw*TeeGovdeNet*(TeeGovdeNet/2))/(tw*TeeGovdeNet + bf*tf);

IxteeNet = tw*(TeeGovdeNet-ycnet)^3/12 + (TeeGovdeNet-ycnet)*tw*((TeeGovdeNet-ycnet)/2)^2 + tw*ycnet^3/12 + ycnet*tw*(ycnet/2)^2 + tf*bf*(TeeGovdeNet + tf/2 -ycnet)^2 + (tf^3)*bf/12;

yccrit = (bf*tf*(TeeGovdeCrit+tf/2) + tw*TeeGovdeCrit*(TeeGovdeCrit/2))/(tw*TeeGovdeCrit + bf*tf);

IxteeCrit = tw*(TeeGovdeCrit-yccrit)^3/12 + (TeeGovdeCrit-yccrit)*tw*((TeeGovdeCrit-yccrit)/2)^2 + tw*yccrit^3/12 + yccrit*tw*(yccrit/2)^2 + tf*bf*(TeeGovdeCrit + tf/2 -yccrit)^2 + (tf^3)*bf/12;

SxTeeNetTop = IxteeNet/(dtnet - ycnet);
SxTeeNetBot = IxteeNet/(ycnet);

SxTeeCritTop = IxteeCrit/(dtcrit - yccrit);
SxTeeCritBot = IxteeCrit/(yccrit);

if bf*tf > (dtnet-tf)*tw;
ypNet = (bf*tf-dtnet*tw+tf*tw)/(2*bf);
xPnet = dtnet - tf + ypNet;
ZpxNet = (bf*(tf-ypNet)^2)/2 + (dtnet-tf)*tw*(dtnet-tf+ypNet)/2 + bf*(ypNet^2)/2;
else
ypNet = (dtnet*tw-tf*tw-bf*tw)/(2*tw);
xPnet = dtnet - tf - ypNet;
ZpxNet = (tw*(dtnet-tf-ypNet)^2)/2 + (tw*ypNet^2)/2 + bf*tf*(tf+ypNet);
end

if bf*tf > (dtcrit-tf)*tw
ypCrit = (bf*tf-dtcrit*tw+tf*tw)/(2*bf);
xPcrit = dtcrit - tf + ypCrit;
ZpxNet = (bf*(tf-ypCrit)^2)/2 + (dtcrit-tf)*tw*(dtcrit-tf+ypCrit)/2 + bf*(ypCrit^2)/2;
else
ypCrit = (dtcrit*tw-tf*tw-bf*tw)/(2*tw);
xPcrit = dtcrit - tf - ypCrit;
ZpxNet = (tw*(dtcrit-tf-ypCrit)^2)/2 + (tw*ypCrit^2)/2 + bf*tf*(tf+ypCrit);
end

IyteeNet = (tf*bf^3)/12 + ((dtnet-tf)*tw^3)/12;
IyteeCrit = (tf*bf^3)/12 + ((dtcrit-tf)*tw^3)/12;

rxteeNet = sqrt(IxteeNet/AteeNet);
ryteeNet = sqrt(IyteeNet/AteeNet);

rxteeCrit = sqrt(IxteeCrit/AteeCrit);
ryteeCrit = sqrt(IyteeCrit/AteeCrit);

JteeNet = (bf*tf^3 + (dtnet-tf/2)*tw^3)/3;
JteeCrit = (bf*tf^3 + (dtcrit-tf/2)*tw^3)/3;

yUssu = dg/2;
xUssu = dg/2;

deffectnet = dg - 2*(dtnet - ycnet);
deffectcrit = dg - 2*(dtcrit - yccrit);

IxNet = 2*IxteeNet + 2*AteeNet*(deffectnet/2)^2;
IxCrit = 2*IxteeCrit + 2*AteeCrit*(deffectcrit/2)^2;
Iy = 1.18;

SxNet = IxNet/(dg/2);
SxCrit = IxCrit/(dg/2);

ZxNet = 2*AteeNet*(deffectnet/2);
ZxCrit = 2*AteeCrit*(deffectcrit/2);

% Gross Seciton Properties

Agross = AteeNet + Do*tw;

Ixgross = IxNet + (tw*Do^3)/12;

Sxgross = Ixgross/(dg/2);

Qcrit = (dtcrit - yccrit - tf/2)*tf*(bf^2)/2;

ecrit = (dtcrit-tf/2)*(Qcrit*((bf^2)/2)*tf)/IyteeCrit;

%yoCr = ecrit - yccrit;

yoCr = yccrit - tf/2;

roCr = sqrt(yoCr^2 + (IxteeCrit + IyteeCrit)/AteeCrit);


% Narinlik Kontrolleri

NarinDegil = 0;
if (bf/2)/tf > 0.56*sqrt(E/Fy)
 fprintf('(bf/2)/tf > 0.56*sqrt(E/Fy) ---> %.2f > %.2f Kesit Narin\n',(bf/2)/tf,0.56*sqrt(E/Fy))
else
 fprintf('(bf/2)/tf < 0.56*sqrt(E/Fy) ---> %.2f < %.2f Kesit Narin Deðil\n',(bf/2)/tf,0.56*sqrt(E/Fy)) 
 NarinDegil = 1;
end

Kompakt = 0;
if dtnet/tw > 0.84*sqrt(E/Fy)
 fprintf('d/tw > 0.84*sqrt(E/Fy) ---> %.2f > %.2f Kesit Kompakt Deðil\n',dtnet/tw,0.84*sqrt(E/Fy))
else
 fprintf('d/tw < 0.84*sqrt(E/Fy) ---> %.2f < %.2f Kesit Kompakt\n',dtnet/tw,0.84*sqrt(E/Fy))
 Kompakt = 1;
end



if Kompakt == 1 && NarinDegil == 1
for i=1:length(KesmeMoment(:,1))

Vr = KesmeMoment(i,2)*0.2248;
Mr = KesmeMoment(i,3)*8.8507;

if i < length(KesmeMoment(:,1))
Mrnext = KesmeMoment(i+1,3)*8.8507;
end

Pr = Mr/deffectcrit;

VrkN = KesmeMoment(i,2);
PrkN = Pr*4.4482;
MrkN = KesmeMoment(i,3);

Mvr = (Vr/2)*(Do/4);

xBukrulmaBoyu = 0.65*(Do/2)/rxteeCrit;
yBukrulmaBoyu = (Do/2)/ryteeCrit;

% Basýnc Hesabý

BurkulmaBoyu = max(yBukrulmaBoyu,xBukrulmaBoyu);

Fe1 = ((pi^2)*E)/(BurkulmaBoyu)^2;

if BurkulmaBoyu <= 4.71*sqrt(E/Fy) || (Fy/Fe1) <= 2.25
    Fcr1 = (0.658^(Fy/Fe1))*Fy;
else
    Fcr1 = 0.887*Fe1;
end

Pn1 = Fcr1*AteeCrit;
Pn1kN = Pn1*4.4482;

% Eðilmeli Burulmalý Burkulma Sýnýr Durumu

Fey = ((pi^2)*E)/(yBukrulmaBoyu)^2;
Fez = ((pi^2)*E/(Do/2)^2 + G*JteeCrit)/(AteeCrit*(roCr^2));


H = 1 - (yoCr^2/roCr^2);


Fe2 = ((Fey + Fez)/(2*H))*(1-sqrt(1-((4*Fey*Fez*H)/(Fey+Fez)^2)));

if yBukrulmaBoyu <= 4.71*sqrt(E/Fy) || (Fy/Fe2) <= 2.25
    Fcr2 = (0.658^(Fy/Fe2))*Fy;
else
    Fcr2 = 0.887*Fe2;
end

Pn2 = Fcr2*AteeCrit;
Pn2kN = Pn2*4.4482;

if Pr < max(Pn1,Pn2)
    %Kontroller{i,1} = fprintf('Pr < max(Pn1,Pn2) = %.2f < %.2f %s için %.2fm de basýnç dayanýmý koþulu saðlandý kesit güvenli \n',Pr,max(Pn1,Pn2),Makaslar{i,1},KesmeMoment(i,1));
else
    Kontroller{i,1} = fprintf('Pr > max(Pn1,Pn2) = %.2f > %.2f %s için %.2fm de basýnç dayanýmý koþulu SAÐLANAMADI kesit GÜVENLÝ DEÐÝL \n',Pr,max(Pn1,Pn2),Makaslar{i,1},KesmeMoment(i,1));
end



% T Kesitte Eðilme Kontrolü

My = Fy*SxTeeCritBot;
Mp = My;
MpkN = Mp*0.113;

if Mvr < Mp
    %fprintf('Mvr < Mp = %.2f < %.2f %s için %.2fm de Eðilme Akmasý koþulu saðlandý kesit güvenli \n',Mvr,Mp,Makaslar{i,1},KesmeMoment(i,1));
else
    fprintf('Mvr > Mp = %.2f > %.2f %s için %.2fm de Eðilme Akmasý koþulu SAÐLANAMADI kesit GÜVENLÝ DEÐÝL \n',Mvr,Mp,Makaslar{i,1},KesmeMoment(i,1));
end


% Gövdede Lokal Burkulma

if dtcrit/tw <= 0.84*sqrt(E/Fy)
    FcrLok = Fy;
elseif 0.84*sqrt(E/Fy) <= dtcrit/tw && dtcrit/tw  <= sqrt(E/Fy)
    FcrLok = (1.43-0.515*(dtcrit/tw)*sqrt(Fy/E))*Fy;
elseif dtcrit/tw >= sqrt(E/Fy)
    FcrLok = 1.52*E/(dtcrit/tw)^2;
end

Mu = 0.9*FcrLok*SxTeeCritBot;

if Mvr < Mu
    %fprintf('Mvr < Mu = %.2f < %.2f %s için %.2f  de Lokal Burkulma Kontrolü Güvenli \n',Mvr,Mu,Makaslar{i,1},KesmeMoment(i,1));
else
    fprintf('Mvr > Mu = %.2f > %.2f %s için %.2f  de Lokal Burkulma Kontrolü GÜVENLÝ DEÐÝL \n',Mvr,Mu,Makaslar{i,1},KesmeMoment(i,1));
end

% Bileþik Etkiler Kontrolü


if Pr/max(Pn2,Pn1) > 0.2
    Ratio = Pr/max(Pn2,Pn1) + Mvr/min(Mu,Mp);
elseif Pr/max(Pn2,Pn1) < 0.2
    Ratio = Pr/max(Pn2,Pn1)*2 + Mvr/min(Mu,Mp);
end


if Ratio < 1
    %fprintf('Ratio > 1 = %.2f < %.2f %s için %.2fm de Basýnç + Eðilme Kontrolü Güvenli \n',Ratio,1,Makaslar{i,1},KesmeMoment(i,1));
elseif Ratio > 1
    fprintf('Ratio < 1 = %.2f > %.2f %s için %.2fm de Basýnç + Eðilme Kontrolü GÜVENLÝ DEÐÝL \n',Ratio,1,Makaslar{i,1},KesmeMoment(i,1));
end

% Boþluklar Arasý Brüt Kesitte Lokal Burkulma Kontrolü

if i < length(KesmeMoment(:,1))
Vrh = (Mrnext - Mr)/deffectcrit;
Mrh = 0.9*(Do/2)*Vrh;

Me = (tw*(S-Do+0.564*Do^2)*Fy)/6;
C1 = 5.097 + 0.1464*(Do/tw) - 0.00174*(Do/tw)^2;
C2 = 1.441 + 0.0625*(Do/tw) - 0.000683*(Do/tw)^2;
C3 = 3.645 + 0.0853*(Do/tw) - 0.00108*(Do/tw)^2;

MallowBoluMe = (C1*(S/Do) - C2*(S/Do)^2 - C3);
Mallow = abs(0.9*MallowBoluMe*Me);

if Mallow > Mrh
    %fprintf('Mallow > Mrh = %.2f > %.2f %s için %.2fm de Gövdede Lokal Burkulma Kontrolü Güvenli \n',Mallow,Mrh,Makaslar{i,1},KesmeMoment(i,1));
elseif Mallow < Mrh
    fprintf('Mallow < Mrh = %.2f < %.2f %s için %.2fm de Gövdede Lokal Burkulma GÜVENLÝ DEÐÝL \n',Mallow,Mrh,Makaslar{i,1},KesmeMoment(i,1));
end

% Boþluklar Arasý Brüt Kesitte Yatay Kesme Kontrolü

Vnh = 0.6*Fy*e*tw;

if Vnh > Vrh
    %fprintf('Vnh > Vrh = %.2f > %.2f %s için %.2fm' de Gövdede Yatay Kesme Kontrolü Güvenli \n',Vnh,Vrh,Makaslar{i,1},KesmeMoment(i,1));
elseif Vnh < Vrh
    fprintf('Vnh < Vrh = %.2f < %.2f %s için %.2fm de Gövdede Yatay Kesme GÜVENLÝ DEÐÝL \n',Vnh,Vrh,Makaslar{i,1},KesmeMoment(i,1));
end

end


% Net Kesitte Dikey Kesme Kontrolü

hBolutw = dtnet/tw;

if hBolutw <= 1.10*sqrt(1.2*E/Fy)
    Cv1 = 1;
elseif hBolutw > 1.10*sqrt(1.2*E/Fy)
    Cv1 = (1.10*sqrt(1.2*E/Fy))/(dtnet/tw);
end


if hBolutw <= 1.10*sqrt(1.2*E/Fy)
    Cv2 = 1;
elseif hBolutw > 1.10*sqrt(1.2*E/Fy) && hBolutw < 1.37*sqrt(1.2*E/Fy)
    Cv2 = (1.10*sqrt(1.2*E/Fy))/(dtnet/tw);
elseif hBolutw > 1.37*sqrt(1.2*E/Fy)
    Cv2 = (1.51*1.2*E)/(((dtnet/tw)^2)*Fy);
end

if hBolutw <= 2.24*sqrt(E/Fy)
phiv = 1;
elseif hBolutw > 2.24*sqrt(E/Fy)
phiv = 0.9;
end

Vnvnet = 0.6*Fy*dtnet*tw*Cv2*2;

if phiv*Vnvnet > Vr
    %fprintf('phiv*Vnvnet > Vr = %.2f > %.2f %s için %.2fm' de Net T Kesitinde Dikey Kesme Kontrolü Güvenli \n',phiv*Vnvnet,Vr,Makaslar{i,1},KesmeMoment(i,1));
elseif phiv*Vnvnet < Vr
    fprintf('phiv*Vnvnet < Vr = %.2f < %.2f %s için %.2fm de Net T Kesitinde Dikey Kesme GÜVENLÝ DEÐÝL \n',phiv*Vnvnet,Vr,Makaslar{i,1},KesmeMoment(i,1));
end

% Brüt Kesitte Dikey Kesme Kontrolü

VnvBrut = 0.6*Fy*dg*tw*Cv1;

if VnvBrut*phiv > Vr
    %fprintf('phiv*VnvBrut > Vr = %.2f > %.2f %s için %.2fm' de Brüt I Kesitinde Dikey Kesme Kontrolü Güvenli \n',VnvBrut*phiv,Vr,Makaslar{i,1},KesmeMoment(i,1));
elseif VnvBrut*phiv < Vr
    fprintf('phiv*VnvBrut < Vr = %.2f < %.2f %s için %.2fm de Brüt I Kesitinde Dikey Kesme GÜVENLÝ DEÐÝL \n',VnvBrut*phiv,Vr,Makaslar{i,1},KesmeMoment(i,1));
end

deplasman = (YukKipft*(1/12)*(39.37*max(KesmeMoment(:,1)))^4)/(384*E*IxNet*0.9);

if deplasman/max(KesmeMoment(:,1)) < max(KesmeMoment(:,1))/180
    %fprintf('delta/L < L/180 %.2f < %.2f %s için Deplasman Kontrolü Güvenli \n',deplasman/max(KesmeMoment(:,1)),max(KesmeMoment(:,1))/180,Makaslar{i,1});
elseif deplasman/max(KesmeMoment(:,1)) >= max(KesmeMoment(:,1))
    fprintf('delta/L > L/180 %.2f > %.2f %s için Deplasman Kontrolü GÜVENLÝ DEÐÝL \n',deplasman/max(KesmeMoment(:,1)),max(KesmeMoment(:,1))/180,Makaslar{i,1});
end


end
end

%Baþlýk Kaynak Uzunluðu Hesabý

StationQty = 39;
ConditionQty = length(KesmeMoment(:,1))/StationQty;

[KesmeMoment,Makaslar] = xlsread('MomentKesme.xlsx','A1:D10000');

for j=0:ConditionQty
for i=1:StationQty
    
    if i ~= StationQty && j ~= ConditionQty
        
    LocIlk = KesmeMoment(j*StationQty+i,1);
    LocSon = KesmeMoment(j*StationQty+i+1,1);
    
    VIlk = KesmeMoment(j*StationQty+i,2);
    Vson = KesmeMoment(j*StationQty+i+1,2);
    
    Vort = ((Vson + VIlk)/2);
    VLoc = (LocSon-LocIlk);
    
    ToKayma(i,j+1) = twmm*VLoc*1000*(Vort*1000*(bfmm*tfmm*(dg*25.4-tfmm/2))/((Ixgross*25.4^4)*twmm));
    
    
    end
    end
end

for i=1:length(ToKayma(1,:))
ToKayma(StationQty,i) = sum(ToKayma(:,i));
end 

VmaxFl = max(abs(ToKayma(StationQty,1:ConditionQty)));

LKaynak = VmaxFl/(3*0.6*420*2);

