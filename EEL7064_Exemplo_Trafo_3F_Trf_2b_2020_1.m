%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformadores trif�sicos     %
% Autor: R. S. Salgado           %
% Data: 14/04/2021               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

format shortEng

% constantes
j = sqrt(-1);
radeg = 180.00/pi;
a = -(1/2) + j*(sqrt(3)/2);
Af = [ 1;  a^2;  a; ];

disp('   *  Valores nominais do transformador trif�sico *')

disp('   *  transformador trif�sico  *')
Stn = 150e3
Vtnat = 2400              % conex�o Delta
Vtnbt = 240*sqrt(3)       % conex�o Y
Itnat = Stn/(sqrt(3)*Vtnat)
Itnbt = Stn/(sqrt(3)*Vtnbt)

disp('   *  Medidas de pot�ncia ativa - m�t. dos 2 watt�metros, transformador trif�sico  *')
disp('   *  Ensaio de CA (lado de BT) *')
W1 = -474.72
W2 = 1503.70
Pw3F = W1 + W2;       %   Perda 3F no n�cleo
Qw3F = sqrt(3)*(W2-W1);
Sw3F = Pw3F + j*Qw3F
Ymbt = conj(Sw3F/3)/(Vtnbt/sqrt(3))^2    % em cada fase da conex�o Y

disp('   *  Ensaio de CC (lado de AT)   *')

W1 = 2128.60
W2 = 721.44
Vcc = 60
Pw3F = W1 + W2;
Qw3F = sqrt(3)*(W2-W1);
Sw3F = Pw3F + j*Qw3F
Zeqat = conj(Sw3F/3)/abs(Itnat/sqrt(3))^2
Zeqbt = Zeqat*((Vtnbt/sqrt(3))/Vtnat)^2

disp('   *  Ensaio CA: trafo 3F, Vnom no lado de BT, lado de AT em CA  *')
disp('   *  trafo 3F: fase "a"  *')
Vna = Vtnbt/sqrt(3)
Ic0 = Vna/((Zeqbt/2)+(1/Ymbt))
Ic0m = abs(Ic0)
Ic0a = radeg*angle(Ic0)
S0 = Vna*conj(Ic0)
S0m = abs(S0)
S0a = radeg*angle(S0)

disp('   *  trafo 3F: fase "a"  *')
Vfbt = (Vtnbt/sqrt(3))*Af;
Vfbtm = abs(Vfbt)
Vfbta = radeg*angle(Vfbt)

disp('   *  Trafo 3F: Fasores tens�o de linha no lado de BT  *')
VLbt = Vfbt*sqrt(3)*(cos(30/radeg) + j*sin(30/radeg));
Vlinha_bt_modulo = abs(VLbt)
Vlinha_bt_angulo = radeg*angle(VLbt)

disp('   *  Trafo 3F: Fasores corrente de linha no lado de BT  *')
ILbt = Ic0*Af;
ILbtm = abs(ILbt)
ILbta = radeg*angle(ILbt)

disp('   *  Trafo 3F: Pot�ncia complexa no lado de BT  *')
S03F = sum(Vfbt.*conj(ILbt))

disp('   *  Trafo 3F: medi��o de pot�ncia ativa via m�todo dos 2 watt�metros no lado de BT  *')
disp('   *  Comprova��o atrav�s da solu��o de um circuito trif�sico  *')
Vab = VLbt(1)
Vabm = abs(Vab)
Vaba = radeg*angle(Vab)
Ia = ILbt(1)
Iam = abs(Ia)
Iaa = radeg*angle(Ia)

Vcb = -VLbt(2)
Vcbm = abs(Vcb)
Vcba = radeg*angle(Vcb)
Ic = ILbt(3)
Icm = abs(Ic)
Ica = radeg*angle(Ic)
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia))
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic))

disp('   *  Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)  *')
P03F = W1 + W2;
Q03F = sqrt(3)*(W2-W1);
S03F = P03F + j*Q03F
Zmbt = abs(Vtnbt/sqrt(3))^2/conj(S03F/3)
Ymbt = 1/Zmbt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('   *  Ensaio CC: trafo 3F, Inom no lado de AT, lado de BT em CC  *')
disp('   *  Trafo 3F: Fasores tens�o de linha no lado de AT  *')
VLat = Vcc*Af*(cos(30/radeg) + j*sin(30/radeg));
VLatm = abs(VLat)
VLata = radeg*angle(VLat)

disp('   *  trafo 3F: tens�es fase-neutro  *')
Vfat = (Vcc/sqrt(3))*Af;
Vfatm = abs(Vfat)
Vfata = radeg*angle(Vfat)

disp('   *  Trafo 3F: Fasores corrente de linha no lado de AT  *')
ILat = Itnat*(cos(angle(Zeqat)) + j*sin(angle(Zeqat)))*Af;
ILatm = abs(ILat)
ILata = radeg*angle(ILat)

disp('   *  Trafo 3F: Pot�ncia complexa no lado de BT  *')
S3F = sum(Vfat.*conj(ILat))

disp('   *  Trafo 3F: medi��o via m�todo dos 2 watt�metros no lado de BT  *')
Vab = VLat(1);
Vabm = abs(Vab)
Vaba = radeg*angle(Vab)
Ia = ILat(1)
Iam = abs(Ia)
Iaa = radeg*angle(Ia)

Vcb = -VLat(2);
Vcbm = abs(Vcb)
Vcba = radeg*angle(Vcb)
Ic = ILat(3)
Icm = abs(Ic)
Ica = radeg*angle(Ic)
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia))
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic))

disp('   *  Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)  *')
Pc3F = W1 + W2;
Qc3F = sqrt(3)*abs((W2-W1));
Sc3F = Pc3F + j*Qc3F
Sc3Fm = abs(Sc3F)
Sc3Fa = radeg*angle(Sc3F)

disp('   *  Trafo 3F: Imped�ncia dos enrolamentos (via pot�ncia complexa total)  *')
Zeqat = Vcc^2 / conj(Sc3F/3)
Zeqbt = (Zeqat*(Vtnbt/(Vtnat))^2)/3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   *********************************************************************************+++++++++**  ')
disp('   *  100% da carga nominal no lado de BT, fator de pot�ncia 0,85 adiantado *')
Scr = 1*Stn/3
Vcr = Vtnbt/sqrt(3)
fpcr = 0.85      % adiantado
Pcr = Scr*fpcr;
Qcr = -Scr*sin(acos(fpcr)); % Aqui mudamos o sinal para representar o FP adiantado
Icr = (Scr/Vcr)*(cos(acos(fpcr))+j*sin(acos(fpcr)));
Icm = abs(Icr)
Ica = (180/pi)*angle(Icr)
Scr = Vcr*conj(Icr)
Scr3F = 3*Scr

disp('   * Imped�ncia equivalente a carga nominal no lado de BT *')
Zcr = abs(Vcr)^2/conj(Scr)
Zcrm = abs(Zcr)
Zcra = radeg*angle(Scr)

disp('   *  tens�o e corrente no n�cleo do trafo  *')
Vn = Icr*(Zeqbt/2) + Vcr;
Vnm = abs(Vn)
Vna = radeg*angle(Vn)

disp('   *  entrada do trafo  *')
I1 = Vn*Ymbt + Icr;
I1m = abs(I1)
I1a = (180/pi)*angle(I1)

V1 = I1*(Zeqbt/2) + Vn;
V1m = abs(V1)
V1a = radeg*angle(V1)

S1 = V1*conj(I1)
cosfig = cos(atan(imag(S1)/real(S1)))
S13F = 3*S1
S13Fm = abs(S13F)
S13Fa = radeg*angle(S13F)

disp('   *  perdas no transformador monof�sico *')
St = abs(I1)^2*(Zeqbt/2) + abs(Icr)^2*(Zeqbt/2) + (conj(Ymbt)*abs(Vn)^2)
St3F = 3*St

disp('   *  rendimento do trafo *')
n = (real(Scr3F)/real(S13F))*100

disp('   *  tens�o de circuito aberto com V1 constante na entrada do autotrafo  *')
V10 = V1*((1/Ymbt)/((1/Ymbt)+(Zeqbt/2)))
disp('   *  regula��o  *')
Reg = 100*(abs(V10)-abs(Vcr))/abs(Vcr)

disp('***********************************************************')

disp('   *  trafo 3F: fase "a"  *')
Vfbt = V1*Af;
Vfbtm = abs(Vfbt)
Vfbta = radeg*angle(Vfbt)

disp('   *  Trafo 3F: Fasores tens�o de linha no lado de BT  *')
VLbt = Vfbt*sqrt(3)*(cos(30/radeg) + j*sin(30/radeg));
Vlinha_bt_modulo = abs(VLbt)
Vlinha_bt_angulo = radeg*angle(VLbt)

disp('   *  Trafo 3F: Fasores corrente de linha no lado de BT  *')
ILbt = I1*Af;
ILbtm = abs(ILbt)
ILbta = radeg*angle(ILbt)

disp('   *  Trafo 3F: Pot�ncia complexa no lado de BT  *')
S13F = sum(Vfbt.*conj(ILbt))

disp('   *  Trafo 3F: medi��o via m�todo dos 2 watt�metros no lado de BT  *')
Vab = VLbt(1)
Vabm = abs(Vab)
Vaba = radeg*angle(Vab)
Ia = ILbt(1)
Iam = abs(Ia)
Iaa = radeg*angle(Ia)

Vcb = -VLbt(2)
Vcbm = abs(Vcb)
Vcba = radeg*angle(Vcb)
Ic = ILbt(3)
Icm = abs(Ic)
Ica = radeg*angle(Ic)
W1 = abs(Vab)*abs(Ia)*cos(angle(Vab)-angle(Ia))
W2 = abs(Vcb)*abs(Ic)*cos(angle(Vcb)-angle(Ic))

disp('   *  Trafo 3F: pot�ncias ativa e reatva totais (m�t. dos 2 watt�metros)  *')
Pw3F = W1 + W2
Qw3F = sqrt(3)*(W2-W1)
Sw3F = Pw3F + j*Qw3F

disp('   *  perdas no transformador trif�sico *')
St3F = 3*abs(I1)^2*(Zeqbt/2) + 3*abs(Icr)^2*(Zeqbt/2) + 3*(conj(Ymbt)*abs(Vn)^2)
