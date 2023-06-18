close all;clear all; clc;

DC = 0.3443;                %Operational profile duty  cycle  (based on case study)
CR = 730;                   %Cycling rate (cycles per year)
T_AE = 14;                  %Ambient temperature, nonoperating (in degrees C)
OperationYear = 10;         %years
%% CIRCUIT BREAKER
B10 = 6000;
lambda_CB = 0.1 * (2/24)/6000 *10^6              %failure [FPMH]

%% INVERTER
% input parameter for every power electronics
f_sw = 1000000;     %switching frequency [Hz]
U_DC = 650;         %DC voltage at the DC side [V DC]
P_B = 4*80;       %ouput power of the battery packs  [kwh]
U_RMS = 230;        %AC terminal voltage 
n_I = 4;            %number of the IGBTs
n_D = 10;           %number of the diodes
n_C = 4;            %number of the capacitors
n_fan = 1;          %number of the fan

% input parameter for IGBT
U_CE0 = 1.43;       %voltage drop on the IGBT
r_CE = 0;           %resistance of IGBT
cosphi = 1;         %angle between voltahe and current
E_on = 2.72e-3;     %energy loss of IGBT in the on-state
E_off = 1.88e-3;    %energy loss of IGBT in the off-state
UI_ref = 650;       % reference commutation voltage of IGBT
II_ref = 260;       % reference commutation current of IGBT
U_CEapplied = 350;
U_CErated = 650;    %[V]
R_TH = 0.33;        %thermal resistance from the heat sink to the ambient environment (from data sheet) [K/W]
R_CH = 3.4;         %thermal resistance from the junction to the interface case of an IGBT/diode

% input parameter for diodes
U_d0 =1.92;         %voltage drop on the Diode 
r_T =0;             %resistance of the diode
E_rec = 1.923e-6;   %reverse recovery energy loss of the diode
UD_ref =650;        %diode's reference commutation voltage
ID_ref =215;        %diode's reference commutation current
U_applied_reverse = 350;
U_rated_reverse = 650;

%Capacitor parameters to enter
C =33*10^-3;        %capacitance [microfarads]
S_A = U_DC/UD_ref;  %Stress ratio, the applied voltage stress divided by the rated voltage

%power losses IGBT & diodes
m = (sqrt(2)*U_RMS)/(sqrt(3)*U_DC/2);   %modulation index
U_L = 0.612*U_DC*m;                     %the AC output line-to-line voltage of the power electronic converter module (related to the DC bus voltage and modulation index m)
I = (sqrt(2)*P_B)/(sqrt(3)*U_L);        %peak phase current at the output interface under the ith operating condition

%for IGBTs & diodes
PI_cond_I = 0.5*(U_CE0*I/pi + r_CE*(I^2)/4) + m* cosphi*(U_CE0*I/8 + r_CE*(I^2)/(3*pi));    %+ during discharge and - during charge
PI_sw_I = (1/pi)*f_sw*(E_on+E_off)*(U_DC*I/(UI_ref*II_ref));
P_IGBT = PI_cond_I+PI_sw_I;                                                                 %IGBT power loss

PD_cond_D = 0.5*(U_d0*I/pi + r_T*(I^2)/4) - m* cosphi*(U_d0*I/8 + r_T*(I^2)/(3*pi));
PD_rec_D = (1/pi)*f_sw*E_rec*(U_DC*I/(UD_ref*ID_ref));
P_D = PD_cond_D + PD_rec_D;                                                                 %diode power loss

T_R = (n_I * P_IGBT+ n_D *P_D) * R_CH;              %ambient relative temperature rise in the interface
T_AO = 40;
T_vj = T_AO + T_R;                                  %junction temperature

%failure rate prediction parameters
%IGBT
lambda_OB_I = 0.000235;             %operation-related failure rate
lambda_EB_I = 0.0001657;            %nonoperation-related (environmental) failure rate 
lambda_TCB_I =  0.00016;            %temperature cycling base failure rate
lambda_IND_I = 0.008899;            %Failure rate, electrical overstress
lambda_SJB_I = 0.0015;              %solder joint failure rate contribution

beta_I = 0.281;                     %growth constant -IGBT
Y_I = 2023 ;                        %manufacture year of a power electronic component -IGBT
t_d =  1993;                        %the year of manufacture of parts on which data was collected
pi_G_I = exp(-beta_I*(Y_I-t_d));    %growth failure rate multiplier -IGBT
k = 8.617*10^-5;
DC_1op_I = 0.23;                    %specific constant-IGBT
DC_1nonop_I = 0.77;                 %specific constant-IGBT
Ea_op_I =  0.2;                     %Activation energy, operating-IGBT
Ea_nonop_I = 0.3;                   %Activation energy, nonoperating-IGBT
pi_TO_I  = exp((-Ea_op_I/k)*((1/(T_AO+T_R+273)-(1/298))));         %operation-related temperature acceleration factor , depends on different Ea
V_s_I = U_CEapplied/U_CErated;                                     %from the datasheet-IGBT
pi_S_I = 0.21*exp(0.31*V_s_I);                                     %Failure rate multiplier, stress-IGBT
pi_SJDT_I =  ((T_AO+T_R-T_AE)/44)^2.26;                            %Failure rate multiplier, solder joint delta temperature
pi_TE_I = exp((-Ea_nonop_I/k)*((1/(T_AE+273)-(1/298))));           %non-operation-related temperature acceleration factor , depends on different Ea
CR_1_I = 754.38;                    %specific constant for a certain power electronic component-IGBT
pi_CR_I = CR/CR_1_I;                %temperature cycling rate acceleration factor-IGBT
pi_DCO_I = DC /DC_1op_I;            %Failure rate multiplier for duty cycle
pi_DCN_I = (1-DC )/DC_1nonop_I;     %non-operation-related duty cycle factor

% DIODES
lambda_OB_D = 0.001603;     %operation-related failure rate -diode
lambda_EB_D = 0.002748;     %nonoperation-related (environmental) failure rate -diode
lambda_TCB_D =  0.02603;    %temperature cycling base failure rate-diode
lambda_IND_D = 0.01158;     %Failure rate, electrical overstress-diode
lambda_SJB_D =0.00021;      %solder joint failure rate contribution-diode

beta_D = 0.223;                     %growth constant -diode
Y_D = 2019 ;                        %manufacture year of a power electronic component-diode
pi_G_D = exp(-beta_D*(Y_D-t_d));    %growth failure rate multiplier
Ea_op_D =  0.3;                     %Activation energy, operating-diode
Ea_nonop_D = 0.4;                   %Activation energy, nonoperating-diode
pi_TO_D  = exp((-Ea_op_D/k)*((1/(T_AO+T_R+273)-(1/298))));         %operation-related temperature acceleration factor , depends on different Ea
%U_s_D = 0.29;                      %default Us-diode
U_s_D = U_applied_reverse/U_rated_reverse;
pi_S_D = 0.21*exp(0.31*U_s_D);                  %Failure rate multiplier, stress-diode
pi_SJDT_D =  ((T_AO+T_R-T_AE)/44)^2.26;         %Failure rate multiplier, solder joint delta temperature-diode
pi_TE_D = exp((-Ea_nonop_D/k)*((1/(T_AE+273)-(1/298))));          %non-operation-related temperature acceleration factor , depends on different Ea
CR_1_D = 736.84;                                %specific constant for a certain power electronic component
pi_CR_D = CR/CR_1_D;            %temperature cycling rate acceleration factor-diode
DC_1op_D = 0.23;                %specific constant-diode
DC_1nonop_D = 0.77;             %specific constant-diode
pi_DCO_D = DC /DC_1op_D;        %Failure rate multiplier for duty cycle-diode
pi_DCN_D = (1-DC )/DC_1nonop_D; %non-operation-related duty cycle factor-diode

%for both IGBT and diodes
DT_1 = 80;                          %component specific constant-IGBT&diodes are the same
pi_DT = ((T_vj - T_AE)/DT_1)^2;     %delta temperature cycling acceleration factor - for both IGBT and diodes

%Failure rate models for the capacitor
lambda_OB_C = 0.000175;         %operation-related failure rate-capacitor
lambda_EB_C = 0.000049;         %nonoperation-related (environmental) failure rate -capacitor
lambda_TCB_C = 0.000032;        %temperature cycling base failure rate-capacitor
lambda_IND_C = 0.000816;        %Failure rate, electrical overstress-capacitor
lambda_SJB_C =0.00095;          %solder joint failure rate contribution

beta_C = 0.229;                 %growth constant -capacitor
Y_C = 2023 ;                    %manufacture year of a power electronic component-capacitor
pi_G_C = exp(-beta_C*(Y_C-t_d));%growth failure rate multiplier
C1 = 7.6;                       %constant
S1 = 0.6;
CE = 0.23;                      %constant 
pi_C_C = (C/C1)^CE;             %capacitance failure rate multiplier
DC_1op_C = 0.17;                %specific constant-capacitor
DC_1nonop_C = 0.83;             %specific constant-capacitor
Ea_op_C =  0.2;
Ea_nonop_C = 0.3;
pi_TO_C  = exp((-Ea_op_C/k)*((1/(T_AO+273)-(1/298))));         %operation-related temperature acceleration factor , depends on different Ea
n_Cap = 17;                                     %constant, function of capacitor
pi_S_C = (S_A/S1)^n_Cap;                        %Failure rate multiplier, stress
pi_SJDT_C =  ((T_AO - T_AE)/44)^2.26;           %Failure rate multiplier, solder joint delta temperature
pi_TE_C = exp((-Ea_nonop_C/k)*((1/(T_AE+273)-(1/298))));          %non-operation-related temperature acceleration factor , depends on different Ea
CR_1_C = 1140.35;                               %specific constant for a certain power electronic component
pi_CR_C = CR/CR_1_C;                            %temperature cycling rate acceleration factor
DT_1_C = 21;                                    %component specific constant
pi_DT_C = ((T_AO - T_AE)/DT_1_C)^2;             %delta temperature cycling acceleration factor
pi_DCO_C = DC /DC_1op_C;                        %Failure rate multiplier for duty cycle
pi_DCN_C = (1-DC )/DC_1nonop_C;                 %non-operation-related duty cycle factor

%failure rates & MTTF - each component - INVERTER
lambda_I = (pi_G_I*(lambda_OB_I * pi_DCO_I * pi_TO_I * pi_S_I + lambda_EB_I * pi_DCN_I * pi_TE_I + lambda_TCB_I * pi_CR_I * pi_DT) + lambda_SJB_I*pi_SJDT_I + lambda_IND_I)             %IGBT failure rate [FPMH]
lambda_D = (pi_G_D*(lambda_OB_D * pi_DCO_D * pi_TO_D * pi_S_D + lambda_EB_D * pi_DCN_D * pi_TE_D + lambda_TCB_D * pi_CR_D * pi_DT) + lambda_SJB_D*pi_SJDT_D + lambda_IND_D)             %Diodes failure rate [FPMH]
lambda_C = (pi_G_C*pi_C_C*(lambda_OB_C * pi_DCO_C * pi_TO_C * pi_S_C + lambda_EB_C * pi_DCN_C * pi_TE_C + lambda_TCB_C * pi_CR_C * pi_DT_C) + lambda_SJB_C*pi_SJDT_C + lambda_IND_C)    %Capacitor failure rate [FPMH]

%fan - known from the datasheet
MTTF_fan =  490000/(DC*365*24);         %years
lambda_fan = (1/MTTF_fan )* 10^6/8760   %FPMH

%failure rates & MTTF - inverter
lambda_inv =(n_D*lambda_D + n_C*lambda_C + n_I*lambda_I );
lambda_invertersystem = 3*lambda_inv
%MTTF_invertersystem= (1/lambda_invertersystem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOFTWARE
KSLOC = 1000;               %Lines of Source Code (in thousands) - assumption
FD = 0.5;                   %Fault Density - SEI's CMM level :5
FL = 2;                     %Fault Latency 
FA = 1;                     %Fault Activation
AS = 0.5;                   %Average Severity 
ts= 48;                     %Time to Stabilization (months)
DSL = 0.01;                 %Defect Stabilization Level 
F0 = KSLOC * FD;            %initial defect density 
kGR = log(1/DSL)/ts;        %growth rate
ti = 12*OperationYear;      %time (in months) after deployment 
F_ti = F0*exp(-kGR*ti);
F_tiMinus1= F0 * exp(-kGR*(ti-1));                     %Number of faults remaining at time, ti-1 , assumption
lambda_SW = (((F_tiMinus1-F_ti)/730)*(DC*FL*FA*AS))    %Software failure rate prediction (in failures per million calendar hours)

%% system level model - INVERTER
PI_P =  0.243;      %Parts Quality process factor 
SS_ESS = 0;         %Screening strength of the screen(s) applied, if any
t = 10;             %Time in years. 
PI_IM = ((t^-0.62)/1.77)*(1-SS_ESS);   %Infant mortality factor 
deltaT = T_AO - T_AE;
G = 10;             %magnitude of random vibration while the system is operating (G_RMS)
PI_E = (0.855*(0.8*(1-exp(-0.065*(deltaT+0.6)^0.6)))+0.2*(1-exp(-0.046*G^1.71)))/0.205;     %Environmental factor
PI_D = 0.094 ;      %Design process factor
Gij = [0, 1, 1, 1, 1, 1, 1];
Wij = [8,8,6,6,4,10,5];
alpha = sum(Gij.*Wij)/sum(Wij);
PI_G = (1.12*(t+2)-alpha)/(2^-alpha);
PI_M = 0.142;       %Manufacturing process factor 
PI_S = 0.036 ;      %Systems Management process factor
PI_I = 0.141 ;      %Induced process factor
PI_N =  0.237;      %No-defect process factor 
PI_W = 0.106;       %Wearout process factor
lambda_P_diode_inv = lambda_D *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W)
lambda_P_IGBT_inv = lambda_I *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W)
lambda_P_cap_inv = lambda_C *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W)
lambda_P_invertersystem = lambda_invertersystem *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W) + 3*n_fan*lambda_fan
MTTF_invSystem = 1/lambda_P_invertersystem                  %*10^6 h

%% DCDC CONVERTER
% input parameter for every power electronics
f_sw_dcdc = 100000;         %switching frequency [Hz]grees C) 
n_I_dcdc = 12;              %number of the IGBTs
n_D_dcdc = 12;              %number of the diodes
n_C_dcdc = 6;               %number of the capacitors

% input parameter for IGBT - DCDC
U_CE0_dcdc = 1.67;          %voltage drop on the IGBT
r_CE_dcdc = 1.5;            %resistance of IGBT
E_on_dcdc = 2.72e-3;        %energy loss of IGBT in the on-state
E_off_dcdc = 1.88e-3;       %energy loss of IGBT in the off-state
UI_ref_dcdc = 650;          % reference commutation voltage of IGBT
II_ref_dcdc = 100;          % reference commutation current of IGBT
U_CEapplied_dcdc = 600;
U_CErated_dcdc = 950;       %[V]
R_TH_dcdc = 0.895;          %thermal resistance from the heat sink to the ambient environment (from data sheet) [K/W]
R_CH_dcdc = 0.66;           %thermal resistance from the junction to the interface case of an IGBT/diode

% input parameter for diodes - DCDC
U_d0_dcdc =1.92;            %voltage drop on the Diode 
r_T_dcdc =1.5;              %resistance of the diode
E_rec_dcdc = 1.1e-3;        %reverse recovery energy loss of the diode
UD_ref_dcdc =650;           %diode's reference commutation voltage
ID_ref_dcdc =100;           %diode's reference commutation current
U_applied_reverse_dcdc = 600;
U_rated_reverse_dcdc = 950;

%Capacitor parameters to enter
C_dcdc =47*10^-3;               %capacitance [microfarads]
S_A_dcdc = U_DC/UD_ref_dcdc;    %Stress ratio, the applied voltage stress divided by the rated voltage

%% power losses IGBT & diodes
PI_cond_I_dcdc = 0.5*(U_CE0_dcdc*I/pi + r_CE_dcdc*(I^2)/4) + m* cosphi*(U_CE0_dcdc*I/8 + r_CE_dcdc*(I^2)/(3*pi));       %+ during discharge and - during charge
PI_sw_I_dcdc = (1/pi)*f_sw_dcdc*(E_on_dcdc+E_off_dcdc)*(U_DC*I/(UI_ref_dcdc*II_ref_dcdc));
P_IGBT_dcdc = PI_cond_I_dcdc+PI_sw_I_dcdc;                                                                              %IGBT power loss

PD_cond_D_dcdc = 0.5*(U_d0_dcdc*I/pi + r_T_dcdc*(I^2)/4) - m* cosphi*(U_d0_dcdc*I/8 + r_T_dcdc*(I^2)/(3*pi));
PD_rec_D_dcdc = (1/pi)*f_sw_dcdc*E_rec_dcdc*(U_DC*I/(UD_ref_dcdc*ID_ref_dcdc));
P_D_dcdc = PD_cond_D_dcdc + PD_rec_D_dcdc;                                                                              %diode power loss

T_R_dcdc = (n_I_dcdc * P_IGBT_dcdc+ n_D_dcdc *P_D_dcdc) * R_CH_dcdc;       %ambient relative temperature rise in the interface
T_vj_dcdc = T_AO + T_R_dcdc;                                               %junction temperature

%%failure rate prediction parameters
%IGBT
pi_TO_I_dcdc  = exp((-Ea_op_I/k)*((1/(T_AO+T_R_dcdc+273)-(1/298))));       %operation-related temperature acceleration factor , depends on different Ea
V_s_I_dcdc = U_CEapplied_dcdc/U_CErated_dcdc;                              %from the datasheet-IGBT
pi_S_I_dcdc = 0.21*exp(0.31*V_s_I_dcdc);                                   %Failure rate multiplier, stress-IGBT
pi_SJDT_I_dcdc =  ((T_AO+T_R_dcdc-T_AE)/44)^2.26;                          %Failure rate multiplier, solder joint delta temperature
pi_TE_I = exp((-Ea_nonop_I/k)*((1/(T_AE+273)-(1/298))));                   %non-operation-related temperature acceleration factor , depends on different Ea

% DIODES
pi_TO_D_dcdc  = exp((-Ea_op_D/k)*((1/(T_AO+T_R_dcdc+273)-(1/298))));       %operation-related temperature acceleration factor , depends on different Ea
U_s_D_dcdc = U_applied_reverse_dcdc/U_rated_reverse_dcdc;
pi_S_D_dcdc = 0.21*exp(0.31*U_s_D_dcdc);                                   %Failure rate multiplier, stress-diode
pi_SJDT_D_dcdc =  ((T_AO+T_R_dcdc-T_AE)/44)^2.26;                          %Failure rate multiplier, solder joint delta temperature-diode

%for both IGBT and diodes
pi_DT_dcdc = ((T_vj_dcdc - T_AE)/DT_1)^2;                                  %delta temperature cycling acceleration factor - for both IGBT and diodes

%%Failure rate models for the capacitor
pi_C_C_dcdc = (C_dcdc/C1)^CE;                                              %capacitance failure rate multiplier
pi_TO_C_dcdc  = exp((-Ea_op_C/k)*((1/(T_AO+273)-(1/298))));                %operation-related temperature acceleration factor , depends on different Ea
n_Cap = 17;                                                                %constant, function of capacitor
pi_S_C_dcdc = (S_A_dcdc/S1)^n_Cap;                                         %Failure rate multiplier, stress
pi_SJDT_C_dcdc =  ((T_AO - T_AE)/44)^2.26;                                 %Failure rate multiplier, solder joint delta temperature
pi_DT_C_dcdc = ((T_AO - T_AE)/DT_1_C)^2;                                   %delta temperature cycling acceleration factor

%%failure rates & MTTF - each component DCDC CONVERTER
lambda_I_dcdc = (pi_G_I*(lambda_OB_I * pi_DCO_I * pi_TO_I_dcdc * pi_S_I_dcdc + lambda_EB_I * pi_DCN_I * pi_TE_I + lambda_TCB_I * pi_CR_I * pi_DT_dcdc) + lambda_SJB_I*pi_SJDT_I_dcdc + lambda_IND_I);  %IGBT failure rate
lambda_D_dcdc = (pi_G_D*(lambda_OB_D * pi_DCO_D * pi_TO_D_dcdc * pi_S_D_dcdc + lambda_EB_D * pi_DCN_D * pi_TE_D + lambda_TCB_D * pi_CR_D * pi_DT_dcdc) + lambda_SJB_D*pi_SJDT_D_dcdc + lambda_IND_D); %Diodes failure rate
lambda_C_dcdc = (pi_G_C*pi_C_C_dcdc*(lambda_OB_C * pi_DCO_C * pi_TO_C_dcdc * pi_S_C_dcdc + lambda_EB_C * pi_DCN_C * pi_TE_C + lambda_TCB_C * pi_CR_C * pi_DT_C_dcdc) + lambda_SJB_C*pi_SJDT_C_dcdc + lambda_IND_C);   %Capacitor failure rate

lambda_P_I_dcdc = lambda_I_dcdc *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W);
lambda_P_D_dcdc = lambda_D_dcdc *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W);
lambda_P_C_dcdc = lambda_C_dcdc*(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W);

%%failure rates & MTTF - DCDC
lambda_DCDC =(n_D_dcdc*lambda_D_dcdc + n_C_dcdc*lambda_C_dcdc + n_I_dcdc*lambda_I_dcdc );
lambda_P_DCDC = lambda_DCDC *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W)

%% TRANSFORMER
%Transformer parameters to enter
T_AO_Tr = T_AO;                     %default from the book
Y_Tr = 2023;                        %manufacture year of a power electronic component

lambda_OB_Tr =  0.0001214;          %operation-related failure rate
lambda_EB_Tr = 0.0001499;           %nonoperation-related (environmental) failure rate 
lambda_TCB_Tr =  0.0000518;         %temperature cycling base failure rate
lambda_IND_Tr = 0.0000386;

beta_Tr = 0;                        %growth constant 
pi_G_Tr = exp(-beta_Tr*(Y_Tr-t_d)); %growth failure rate multiplier
Ea_op_Tr =  0.24;
Ea_nonop_Tr = 0.24;
T_R_Tr = 175;                       %temperature rise
pi_TO_Tr  = exp((-Ea_op_Tr/k)*((1/(T_AO_Tr+T_R_Tr+273)-(1/298))));      %operation-related temperature acceleration factor , depends on different Ea
pi_TE_Tr = exp((-Ea_nonop_Tr/k)*((1/(T_AE+273)-(1/298))));              %non-operation-related temperature acceleration factor , depends on different Ea
CR_1_Tr =  312;                                                         %specific constant for a certain power electronic component

pi_CR_Tr = CR/CR_1_Tr;                                          %temperature cycling rate acceleration factor
DT_1_Tr = 21.94;                                                %component specific constant
pi_DT_Tr = ((T_AO_Tr + T_R_Tr - T_AE)/DT_1_Tr)^2;               %delta temperature cycling acceleration factor
DC_1op_Tr = 0.38;                                               %specific constant
DC_1nonop_Tr = 0.62;                                            %specific constant
pi_DCO_Tr = DC /DC_1op_Tr;                                      %Failure rate multiplier for duty cycle
pi_DCN_Tr = (1-DC )/DC_1nonop_Tr;                               %non-operation-related duty cycle factor
lambda_T = (pi_G_Tr*(lambda_OB_Tr * pi_DCO_Tr * pi_TO_Tr + lambda_EB_Tr * pi_DCN_Tr * pi_TE_Tr + lambda_TCB_Tr * pi_CR_Tr * pi_DT_Tr) + lambda_IND_Tr);

%% system level model

lambda_P_T = lambda_T *(PI_P*PI_IM*PI_E+PI_D*PI_G + PI_M*PI_IM*PI_E*PI_E*PI_G + PI_S*PI_G + PI_I + PI_N +PI_W)

%% CONTROL BOARD
lambdaControlBoard = 5.284;                 %FPMH

%% TMS - assigned value
MTTF_fanTMS =  50000/(DC*365*24);           %years
lambda_TMS = (1/MTTF_fanTMS  )* 10^6/8760   %FPMH

%% DCPM
B10_Cont = 1000000;
lambda_Cont = 0.1 * (0.1429/24)/6000 *10^6	%FPMH
lambda_DCPM = 2* lambda_Cont;

%% BATTERY PACKS
%code for battery packs reliability is provided separately
lambda_batterySystem_FPMH = 15.9949
%% BESS Reliability
lambda_BESS= lambda_batterySystem_FPMH + lambda_P_invertersystem + lambda_CB + lambda_P_T + lambdaControlBoard +lambda_P_DCDC + lambda_TMS + lambda_SW + lambda_DCPM
BESS_MTTF = 1/(lambda_BESS)         %FPMH
BESS_MTTF_year = BESS_MTTF*10^6 / (8760*DC)
%% Lambda [f/year]
lambda_batterySystem_year = lambda_batterySystem_FPMH /10^6 *DC*8760
lambda_DCPM_year = lambda_DCPM/10^6 *DC*8760
lambda_P_invertersystem_year =lambda_P_invertersystem/10^6 *DC*8760
lambda_P_T_year = lambda_P_T/10^6 *DC*8760
lambda_TMS_year  = lambda_TMS/10^6 *DC*8760
lambda_EMS_year = (lambdaControlBoard+lambda_SW)/10^6 *DC*8760
lambda_CB_year = lambda_CB/10^6 *DC*8760
lambda_P_DCDC_year = lambda_P_DCDC/10^6 *DC*8760
lambda_BESS_year = lambda_batterySystem_year+lambda_DCPM_year+lambda_P_invertersystem_year+lambda_P_T_year+lambda_TMS_year+lambda_EMS_year+lambda_CB_year+lambda_P_DCDC_year
MTTF_BESS_year_2 = 1/lambda_BESS_year

%% plot the MTTF
time = 0:0.001:0.15;                  %operating time [hrs]

R_BESS = zeros(length(time),1);
for i=1:length(time)
    R_BESS(i) = exp(-lambda_BESS *time(i));
end

figure(1)
plot(time, R_BESS, 'LineWidth',2);
h=xlabel('Operation Time $-t$ [million hrs] ');
set(h, 'interpreter', 'latex');
h=ylabel('Reliability $ - R_{BESS}(t) $');
set(h, 'interpreter', 'latex');
xline(BESS_MTTF, '--r', 'MTTF_{BESS}','LineWidth',1.0);
set(findall(gcf,'type','axes'),'fontsize',18)
set(findall(gcf,'type','text'),'fontSize',18) 
grid on;
