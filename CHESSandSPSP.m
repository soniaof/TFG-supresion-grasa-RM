
% Comparación CHESS vs SPSP usando simulación de Bloch para ambos
% haber ejecutado test_fat_sat.m y test_2D_spec_espiral.m antes

%% 

%  BLOQUE 0 - Parámetros comunes

freqband = 15;
zoneband = 25;
flip_deg = 90;

fatfreq = -(Sys.fat * 1e-6 * Sys.gamma * Sys.field);

fatval = find(abs(frange - fatfreq) < freqband);
watval = find(abs(frange - 0)       < freqband);
roival = find((abs(frange - fatfreq) < zoneband) | (abs(frange - 0) < zoneband));

%fatval y watval son los índices dentro del vector frange que corresponden a frecuencias
% de grasa y agua respectivamente. Por ejemplo si frange va de -300 a 100 Hz,
% fatval señala las posiciones donde frange está cerca de -223 Hz.

dfat   = [1:10, 35:75, 99:118]; %columnas del espacio X donde hay grasa y agua. 
%Son números entre 1 y 128 que indican posiciones espaciales.
dwater = [11:34, 76:98, 119:128];

% Target y ROI (Nf x Nx)
dtarget = zeros(Nf, Nx);
dtarget(fatval, dfat) = flip_deg * pi/180; 
%matriz de 41 filas (frecuencias) por 128 columnas (posiciones X),
% donde se pone 90°×π/180 en las celdas que corresponden a frecuencia de grasa y
% posición de grasa, y cero en todo lo demás.

roi = logical(zeros(Nf, Nx));
roi(roival, :) = true;

fprintf('Frecuencia de grasa: %.1f Hz\n', fatfreq);

%% 

% BLOQUE 1 - Parámetros de simulación Bloch

% Gradientes del espiral de Arquímedes convertidos a G/cm
% (igual que en test_2D_spec_espiral.m sección "Cambiar dtarget/ROI y simular con Bloch")
gxsim = gx(:);   % [G/cm]
gysim = gy(:);   % [G/cm]
gzsim = zeros(size(gxsim));

% Tejido: grasa (igual en CHESS y SPSP)
T1sim = T1.fat / 1000;   % [s]
T2sim = T2.fat / 1000;   % [s]
sens  = ones([1 length(x)]); %sensibilidad del coil de RF, que aquí se pone a 1
% en todos los puntos porque no estamos modelando inhomogeneidades de antena.
%% 

%  BLOQUE 2 - Magnetización CHESS con Bloch
%  (ya la tenemos de test_fat_sat.m, solo transponemos)

% Mxyspec viene del bucle multiespectral de test_fat_sat.m
% tamaño (Nx x Nf) -> transponer a (Nf x Nx) para que coincida con dtarget
Mxy_chess = Mxyspec.';   % (Nf x Nx)

fprintf('Variación espacial CHESS (debe ser ~0): %.6f\n', ...
    std(abs(Mxy_chess(fatval(1),:))));

%% 

%  BLOQUE 3 - Pulso SPSP óptimo para flip=90°
%  (recalculamos con el mismo flip_deg que CHESS)

% Recalcular bopt para flip=90° (el mismo ángulo que usamos en CHESS)
% Usamos la misma fórmula WLS que en test_2D_spec.m pero sin diag() gigante
dtarget_spsp = zeros(Nf, Nx);
dtarget_spsp(fatval, dfat) = flip_deg * pi/180;

roi_spsp = logical(zeros(Nf, Nx));
roi_spsp(roival, :) = true;

% Optimización eficiente (sin crear matriz diag gigante)
Aw   = A .* roi_spsp(:);
AtA  = Aw' * A;
Atd  = A'  * (roi_spsp(:) .* dtarget_spsp(:));
bopt_spsp = (AtA + beta * eye(Nt)) \ Atd;
%Esto calcula el pulso RF óptimo para SPSP. La fórmula matemática completa es:
%bopt = inv(A'·W·A + β·I) · A'·W·d
%Donde A es la matriz del sistema que relaciona pulso con magnetización,
% W es la ROI como peso, d es el target y β es la regularización que controla
% cuánta energía puede tener el pulso.
fprintf('Pulso SPSP calculado. Norma: %.4f\n', norm(bopt_spsp));

%% 

%  BLOQUE 4 - Magnetización SPSP con Bloch

% ANTES usábamos la aproximación lineal:
%   Mxy_spsp = mopt = reshape(A*bopt, [Nf Nx])
%
% AHORA usamos el simulador de Bloch, igual que CHESS:
%   para cada frecuencia del rango, simulamos blochCim con el pulso bopt_spsp

fprintf('Simulando SPSP con Bloch (%d frecuencias)...\n', Nf);

Mxyspec_spsp = zeros(Nx, Nf);   % (Nx x Nf), igual que Mxyspec de CHESS

for ii = 1:Nf
    freq_ii = frange(ii) * ones(size(x));   % todos los puntos X con la misma frecuencia
    
    % blochCim: simula la magnetización con el pulso bopt_spsp
    % igual que se hace en CHESS pero con el pulso SPSP
    [Mx, My, ~] = blochCim(bopt_spsp, [gxsim, gysim, gzsim], ...
        Sys.dt, T1sim, T2sim, freq_ii, x, 0, sens);
    
    Mxyspec_spsp(:, ii) = Mx + 1i * My;
end

%El bucle recorre las 41 frecuencias de frange. Para cada una,
% llama a blochCim que resuelve las ecuaciones diferenciales de Bloch físicamente.
% Se le pasa el pulso bopt_spsp, los gradientes, el paso temporal, T1, T2,
% la frecuencia de esa iteración y las posiciones X. Te devuelve Mx y My que son
% las componentes de la magnetización transversal. Se combinan en la magnetización
% compleja Mxy = Mx + i·My.

% Transponer a (Nf x Nx) para que coincida con dtarget
Mxy_spsp_bloch = Mxyspec_spsp.';   % (Nf x Nx)

fprintf('Simulación Bloch SPSP completada.\n');

%% 

%  BLOQUE 5 - Calcular NRMSE de ambas técnicas
%  (ahora ambas usan Bloch)

func_nrmse = @(m, d, r) norm(abs(d(r)) - abs(m(r))) / norm(abs(d(r)));

nrmse_chess      = func_nrmse(Mxy_chess,      dtarget, roi);
nrmse_spsp_bloch = func_nrmse(Mxy_spsp_bloch, dtarget, roi);

% También calculamos el NRMSE de la aproximación STA para comparar
% con el Bloch y ver cuánto error introduce la aproximación
Mxy_spsp_sta = mopt;   % aproximación A*b que ya teníamos
nrmse_spsp_sta = func_nrmse(Mxy_spsp_sta, dtarget, roi);

fprintf('\n=== RESULTADOS ===\n');
fprintf('NRMSE CHESS:            %.4f (%.1f%%)\n', nrmse_chess,      nrmse_chess*100);
fprintf('NRMSE SPSP Bloch:       %.4f (%.1f%%)\n', nrmse_spsp_bloch, nrmse_spsp_bloch*100);
fprintf('NRMSE SPSP aprox. STA:  %.4f (%.1f%%)\n', nrmse_spsp_sta,   nrmse_spsp_sta*100);
fprintf('\nDiferencia Bloch vs STA: %.4f\n', abs(nrmse_spsp_bloch - nrmse_spsp_sta));

%% 

%  BLOQUE 6 - Figura comparativa principal
%  Target | CHESS (Bloch) | SPSP (Bloch)

figure('Name','Comparación CHESS vs SPSP - Bloch', 'Position', [50 50 1400 700]);
sgtitle('Comparación CHESS vs SPSP - Simulación Bloch', 'FontSize', 14);

clim_val = [0 flip_deg*pi/180];

% --- Fila 1: Magnetización completa ---
subplot(2,3,1);
imagesc(x, frange, dtarget, clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('TARGET');
hold on;
yline(fatfreq, 'r--', 'Grasa', 'LineWidth', 1.5);
yline(0, 'b--', 'Agua', 'LineWidth', 1.5);
hold off;

subplot(2,3,2);
imagesc(x, frange, abs(Mxy_chess), clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('CHESS - Bloch');
hold on;
yline(fatfreq, 'r--', 'LineWidth', 1.5);
yline(0, 'b--', 'LineWidth', 1.5);
hold off;

subplot(2,3,3);
imagesc(x, frange, abs(Mxy_spsp_bloch), clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('SPSP - Bloch');
hold on;
yline(fatfreq, 'r--', 'LineWidth', 1.5);
yline(0, 'b--', 'LineWidth', 1.5);
hold off;

% --- Fila 2: Solo ROI con NRMSE ---
subplot(2,3,4);
imagesc(x, frange, dtarget .* roi, clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('TARGET en ROI');

subplot(2,3,5);
imagesc(x, frange, abs(Mxy_chess) .* roi, clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('CHESS en ROI\nNRMSE = %.3f', nrmse_chess));

subplot(2,3,6);
imagesc(x, frange, abs(Mxy_spsp_bloch) .* roi, clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('SPSP Bloch en ROI\nNRMSE = %.3f', nrmse_spsp_bloch));

colormap('parula');

%% 

%  BLOQUE 7 - Comparar Bloch vs aproximación STA en SPSP
%  (figura extra que muestra cuánto error introduce la aproximación)

figure('Name','SPSP: Bloch vs aproximación STA', 'Position', [50 50 1200 400]);
sgtitle('SPSP: Simulación Bloch vs Aproximación STA (A*b)');

subplot(1,3,1);
imagesc(x, frange, abs(Mxy_spsp_sta), clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('SPSP aproximación STA\nNRMSE = %.3f', nrmse_spsp_sta));

subplot(1,3,2);
imagesc(x, frange, abs(Mxy_spsp_bloch), clim_val);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('SPSP simulación Bloch\nNRMSE = %.3f', nrmse_spsp_bloch));

subplot(1,3,3);
imagesc(x, frange, abs(abs(Mxy_spsp_bloch) - abs(Mxy_spsp_sta)));
xlabel('X [cm]'); ylabel('Frec. [Hz]');
colorbar;
title('Diferencia |Bloch - STA|');
% Si la diferencia es pequeña, la aproximación STA es buena
% Si es grande, el ángulo de giro es demasiado alto para STA

colormap('parula');

%% 

%  BLOQUE 8 - Perfil espectral en punto central

idx_centro = 50;   % columna 50 está dentro de dfat=[35:75]

figure('Name','Perfil espectral central - Bloch');
plot(frange, abs(Mxy_chess(:, idx_centro)),      'b-',  'LineWidth', 2, ...
    'DisplayName', 'CHESS (Bloch)');
hold on;
plot(frange, abs(Mxy_spsp_bloch(:, idx_centro)), 'r-',  'LineWidth', 2, ...
    'DisplayName', 'SPSP (Bloch)');
plot(frange, abs(Mxy_spsp_sta(:, idx_centro)),   'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'SPSP (STA aprox)');
xline(fatfreq, 'k:', 'Grasa', 'LineWidth', 1.5);
xline(0, 'g:',  'Agua',  'LineWidth', 1.5);
hold off;
xlabel('Frecuencia [Hz]'); ylabel('|M_{xy}|');
title(sprintf('Perfil espectral en x = %.1f cm', x(idx_centro)));
legend; grid on; ylim([0 2]);

%CHESS suprime la grasa espectralmente bien (pico alto en -220 Hz) pero es ciego
% en el espacio X, SPSP también suprime la grasa pero con menor amplitud porque
% tiene que cumplir restricciones espaciales además de espectrales
%La aproximación STA subestima mucho la magnetización real para ángulos grandes
% como 90°, por eso es importante usar Bloch para evaluación final
%% 

%  BLOQUE 9 - Resumen

fprintf('\n========================================\n');
fprintf('        RESUMEN COMPARATIVO\n');
fprintf('========================================\n');
fprintf('%-25s %-10s %-12s %-12s\n', 'Métrica', 'CHESS', 'SPSP Bloch', 'SPSP STA');
fprintf('%-25s %-10.4f %-12.4f %-12.4f\n', 'NRMSE', ...
    nrmse_chess, nrmse_spsp_bloch, nrmse_spsp_sta);
fprintf('%-25s %-10s %-12s %-12s\n', 'Selectividad X',    'NO',  'SÍ', 'SÍ');
fprintf('%-25s %-10s %-12s %-12s\n', 'Selectividad freq', 'SÍ',  'SÍ', 'SÍ');
fprintf('%-25s %-10s %-12s %-12s\n', 'Simulación',        'Bloch','Bloch','Aprox');
fprintf('========================================\n');