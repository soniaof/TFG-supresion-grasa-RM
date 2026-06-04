function T2=t2_values_1_5
% Provides T2 value structure, all in millisecond [msec]
% T2 values taken from 
% https://mri-q.com/why-is-t1--t2.html

    T2.fat=70;                     % grasa
    T2.liver=40;                   % higado
    T2.muscle=50;                  % musculo
    T2.wm=75;                       % sustancia blanca
    T2.gm=90;                      % sustancia gris
    T2.csf=2000;                   % liquido cefalorraquideo (LCR)

return