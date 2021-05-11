%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformadores trifasicos     %
% Autor: R. S. Salgado           %
% Data: 14/04/2021               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

format shortEng;

% constantes
j = sqrt(-1);
radeg = 180.00/pi;
a = -(1/2) + j*(sqrt(3)/2);
Af = [ 1;  a^2;  a; ];

% Valores nominais do transformador trif�sico *')

%transformador trif�sico
Stn = 150e3;
Vtnat = 2400;              % conex�o Delta
Vtnbt = 240*sqrt(3);       % conex�o Y
Itnat = Stn/(sqrt(3)*Vtnat);
Itnbt = Stn/(sqrt(3)*Vtnbt);

%Medidas de potencia ativa - m�t. dos 2 watt-metros, transformador trif�sico
%Ensaio de CA (lado de BT) *')
W1 = -474.72;
W2 = 1503.70;
Pw3F = W1 + W2;       %   Perda 3F no n�cleo
Qw3F = sqrt(3)*(W2-W1);
Sw3F = Pw3F + j*Qw3F;
Ymbt = conj(Sw3F/3)/(Vtnbt/sqrt(3))^2;    % em cada fase da conex�o Y

%Ensaio de CC (lado de AT) 

W1 = 2128.60;
W2 = 721.44;
Vcc = 60;
Pw3F = W1 + W2;
Qw3F = sqrt(3)*(W2-W1);
Sw3F = Pw3F + j*Qw3F;
Zeqat = conj(Sw3F/3)/abs(Itnat/sqrt(3))^2;
Zeqbt = Zeqat*((Vtnbt/sqrt(3))/Vtnat)^2;

%Ensaio CA: trafo 3F, Vnom no lado de BT, lado de AT em CA
%trafo 3F: fase "a"
Vna = Vtnbt/sqrt(3);
Ic0 = Vna/((Zeqbt/2)+(1/Ymbt));
Ic0m = abs(Ic0);
Ic0a = radeg*angle(Ic0);
S0 = Vna*conj(Ic0);
S0m = abs(S0);
S0a = radeg*angle(S0);

%trafo 3F: fase "a"
Vfase_bt = (Vtnbt/sqrt(3))*Af;
Vfase_bt_modulo = abs(Vfase_bt);
Vfase_bt_angulo = radeg*angle(Vfase_bt);

%Trafo 3F: Fasores tensao de linha no lado de BT
Vlinha_bt = Vfase_bt*sqrt(3)*(cos(30/radeg) + j*sin(30/radeg));
Vlinha_bt_modulo = abs(Vlinha_bt);
Vlinha_bt_angulo = radeg*angle(Vlinha_bt);

%Trafo 3F: Fasores corrente de linha no lado de BT
Ilinha_bt = Ic0*Af;
Ilinha_bt_modulo = abs(Ilinha_bt);
Ilinha_bt_angulo = radeg*angle(Ilinha_bt);

%Trafo 3F: Potencia complexa no lado de BT
S03F = sum(Vfase_bt.*conj(Ilinha_bt));

%Trafo 3F: medicao de potencia ativa via metodo dos 2 wattimetros no lado de BT
%Comprovacao atrav�s da solu��o de um circuito trif�sico
Vab = Vlinha_bt(1);
Vabm = abs(Vab);
Vaba = radeg*angle(Vab);
Ia = Ilinha_bt(1);
Iam = abs(Ia);
Iaa = radeg*angle(Ia);

Vcb = -Vlinha_bt(2);
Vcbm = abs(Vcb);
Vcba = radeg*angle(Vcb);
Ic = Ilinha_bt(3);
Icm = abs(Ic);
Ica = radeg*angle(Ic);
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia));
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic));

%Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)
P03F = W1 + W2;
Q03F = sqrt(3)*(W2-W1);
S03F = P03F + j*Q03F;
Zmbt = abs(Vtnbt/sqrt(3))^2/conj(S03F/3);
Ymbt = 1/Zmbt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ensaio CC: trafo 3F, Inom no lado de AT, lado de BT em CC
%Trafo 3F: Fasores tens�o de linha no lado de AT
VLat = Vcc*Af*(cos(30/radeg) + j*sin(30/radeg));
VLatm = abs(VLat);
VLata = radeg*angle(VLat);

%trafo 3F: tens�es fase-neutro
Vfat = (Vcc/sqrt(3))*Af;
Vfatm = abs(Vfat);
Vfata = radeg*angle(Vfat);

%Trafo 3F: Fasores corrente de linha no lado de AT
ILat = Itnat*(cos(angle(Zeqat)) + j*sin(angle(Zeqat)))*Af;
ILatm = abs(ILat);
ILata = radeg*angle(ILat);

%Trafo 3F: Pot�ncia complexa no lado de BT
S3F = sum(Vfat.*conj(ILat));

%Trafo 3F: medi��o via m�todo dos 2 watt�metros no lado de BT
Vab = VLat(1);
Vabm = abs(Vab);
Vaba = radeg*angle(Vab);
Ia = ILat(1);
Iam = abs(Ia);
Iaa = radeg*angle(Ia);

Vcb = -VLat(2);
Vcbm = abs(Vcb);
Vcba = radeg*angle(Vcb);
Ic = ILat(3);
Icm = abs(Ic);
Ica = radeg*angle(Ic);
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia));
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic));

%Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)
Pc3F = W1 + W2;
Qc3F = sqrt(3)*abs((W2-W1));
Sc3F = Pc3F + j*Qc3F;
Sc3Fm = abs(Sc3F);
Sc3Fa = radeg*angle(Sc3F);

%Trafo 3F: Imped�ncia dos enrolamentos (via pot�ncia complexa total)
Zeqat = Vcc^2 / conj(Sc3F/3);
Zeqbt = (Zeqat*(Vtnbt/(Vtnat))^2)/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************************+++++++++**  ')
%100% da carga nominal no lado de BT, fator de pot�ncia 0,85 atrasado *')
atraso = 1 % -1 para adiantado, 1 para atrasado
carga_nominal = .85;
Scr = carga_nominal*Stn/3;
tensao_carga = Vtnbt/sqrt(3);
fator_pot_carga = 0.85;
Pcr = Scr*fator_pot_carga;
Qcr = atraso*Scr*sin(acos(fator_pot_carga));
Icr = (Scr/tensao_carga)*(cos(acos(fator_pot_carga))-atraso*j*sin(acos(fator_pot_carga)));
Icm = abs(Icr);
Ica = (180/pi)*angle(Icr);
Scr = tensao_carga*conj(Icr);
Scr3F = 3*Scr;

%mped�ncia equivalente a carga nominal no lado de BT *')
Zcr = abs(tensao_carga)^2/conj(Scr);
Zcrm = abs(Zcr);
Zcra = radeg*angle(Scr);

%tens�o e corrente no n�cleo do trafo
Vn = Icr*(Zeqbt/2) + tensao_carga;
Vnm = abs(Vn);
Vna = radeg*angle(Vn);

%entrada do trafo
I1 = Vn*Ymbt + Icr;
I1m = abs(I1);
I1a = (180/pi)*angle(I1);

V1 = I1*(Zeqbt/2) + Vn;
V1m = abs(V1);
V1a = radeg*angle(V1);

S1 = V1*conj(I1)
cosfig = cos(atan(imag(S1)/real(S1)));
S13F = 3*S1;
S13Fm = abs(S13F);
S13Fa = radeg*angle(S13F);

%perdas no transformador monof�sico *')
St = abs(I1)^2*(Zeqbt/2) + abs(Icr)^2*(Zeqbt/2) + (conj(Ymbt)*abs(Vn)^2);
St3F = 3*St;

%rendimento do trafo *')
n = (real(Scr3F)/real(S13F))*100;

%tens�o de circuito aberto com V1 constante na entrada do autotrafo
V10 = V1*((1/Ymbt)/((1/Ymbt)+(Zeqbt/2)));
%regula��o
Reg = 100*(abs(V10)-abs(tensao_carga))/abs(tensao_carga);

%***********************************************************

%trafo 3F: fase "a"
Vfase_bt = V1*Af;
Vfase_bt_modulo = abs(Vfase_bt);
Vfase_bt_angulo = radeg*angle(Vfase_bt);

%Trafo 3F: Fasores tens�o de linha no lado de BT
Vlinha_bt = Vfase_bt*sqrt(3)*(cos(30/radeg) + j*sin(30/radeg));
Vlinha_bt_modulo = abs(Vlinha_bt);
Vlinha_bt_angulo = radeg*angle(Vlinha_bt);

%Trafo 3F: Fasores corrente de linha no lado de BT
Ilinha_bt = I1*Af;
Ilinha_bt_modulo = abs(Ilinha_bt);
Ilinha_bt_angulo = radeg*angle(Ilinha_bt);

%Trafo 3F: Pot�ncia complexa no lado de BT
S13F = sum(Vfase_bt.*conj(Ilinha_bt));

%Trafo 3F: medi��o via m�todo dos 2 watt�metros no lado de BT
Vab = Vlinha_bt(1);
Vabm = abs(Vab);
Vaba = radeg*angle(Vab);
Ia = Ilinha_bt(1);
Iam = abs(Ia);
Iaa = radeg*angle(Ia);

Vcb = -Vlinha_bt(2);
Vcbm = abs(Vcb);
Vcba = radeg*angle(Vcb);
Ic = Ilinha_bt(3);
Icm = abs(Ic);
Ica = radeg*angle(Ic);
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia));
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic));

%Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)
Pw3F = W1 + W2;
Qw3F = sqrt(3)*(W2-W1);
Sw3F = Pw3F + j*Qw3F;

%perdas no transformador trif�sico *')
St3F = 3*abs(I1)^2*(Zeqbt/2) + 3*abs(Icr)^2*(Zeqbt/2) + 3*(conj(Ymbt)*abs(Vn)^2);
