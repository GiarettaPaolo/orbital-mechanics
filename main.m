% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRODUZIONE ALL'ANALISI DI MISSIONI SPAZIALI%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Final project%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Due date: 30 / 11 / 2022
% Presentation date: 5 / 12 / 2022
% 
% Paolo Giaretta
% Giovanni Licitra
% Francesco Monai
% ASSIGNED DATA (INITIAL AND FINAL ORBITS)
%
% r1 = [2254.3254; -8092.3126; -4199.8027];
% v1 = [5.6120; 2.4220; -1.7020];
% 
% a2 = 16410.0000;
% e2 = 0.2678;
% i2 = 0.5612;
% OM2 =0.4075;
% om2 = 1.0700;
% th2 = 1.3420;
%
%
%
close all
addpath(genpath("functions"))
addpath(genpath("transfers"))
%% PLOT INITIAL AND FINAL ORBITS
Orbit_plots

%% 1. STANDARD TRASNFER
disp("\nSTandard Transfer")
standard_transfer
% 
%% 2. ALTERNATIVE TRANSFER N1
disp("\nFirst Transfer")
alternative_transfer_1
% 
%% 3. ALTERNATIVE TRANSFER N2
disp("\nSecond Transfer")
alternative_transfer_2
% 
%% 4. ALTERNATIVE TRANSFER N3
disp("\nThird Transfer")
alternative_transfer_3
% 
%% 5. ALTERNATIVE TRANSFER N4
disp("\nFourth Transfer")
alternative_transfer_4
% 
%% 6. ALTERNATIVE TRANSFER N5
disp("\nFifth Transfer")
alternative_transfer_5
% 
%% 7. ALTERNATIVE TRANSFER N6
disp("\nSisxth Transfer")
alternative_transfer_6
%
%
% extra_transfers_(1-5) are supplementary transfers or optimization codes
% that have been excluded from the report or deemed unnecessary.












