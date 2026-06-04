% chess_phantom.m
%
% ANTES haber ejecutado:
%   test_fat_sat.m -> genera: Mxyspec, frange, x, Sys, T1, T2

%% BLOQUE 0 - Cargar phantom y definir ROIs
load('/Users/sonita/Desktop/4º/TFG/Datos TFG/phantomdigital.mat');

roi_grasa    = logical(tejidos(:,:,3));   % capa 3 = anillo = GRASA
roi_interior = logical(tejidos(:,:,2));   % capa 2 = interior = AGUA
roi_fondo    = logical(tejidos(:,:,1));   % capa 1 = fondo/aire

figure('Name','ROIs del phantom - CHESS');
mapa_roi = zeros(128,128);
mapa_roi(roi_interior) = 1;
mapa_roi(roi_grasa)    = 2;
imagesc(mapa_roi); axis image;
colormap(gca, [0.1 0.1 0.1; 0.2 0.4 0.8; 0.9 0.2 0.2]);
colorbar;
title('ROIs: negro=fondo | azul=agua | rojo=grasa');
xlabel('X [px]'); ylabel('Y [px]');

fprintf('Píxeles grasa (capa 3):    %d\n', sum(roi_grasa(:)));
fprintf('Píxeles interior (capa 2): %d\n', sum(roi_interior(:)));

%% BLOQUE 1 - Parámetros
fatfreq = -(Sys.gamma * Sys.field * Sys.fat * 1e-6);
freqband = 15;
fatval   = find(abs(frange - fatfreq) < freqband);
watval   = find(abs(frange - 0)       < freqband);

fprintf('\nFrecuencia de grasa: %.1f Hz\n', fatfreq);
%% BLOQUE 1 - Parámetros
fatfreq_pos = Sys.gamma * Sys.field * Sys.fat * 1e-6;  % +223.5 Hz
fatfreq     = -abs(fatfreq_pos);                        % fuerza a -223.5 Hz 
freqband    = 15;

fatval = find(abs(frange - fatfreq) < freqband);
watval = find(abs(frange - 0)       < freqband);

fprintf('fatfreq = %.4f Hz\n', fatfreq);
fprintf('fatval: %d puntos\n', numel(fatval));

%% BLOQUE 2 - Respuesta espectral de CHESS
% Mxyspec viene de test_fat_sat.m, tamaño (Nx x Nf) = (128 x 41)
% Transponer a (Nf x Nx) = (41 x 128)
Mxy_chess = Mxyspec.';

% CHESS no tiene gradientes espaciales -> respuesta igual en todos los X
% Cogemos la columna central como perfil representativo
[~, idx_centro] = min(abs(x));
perfil_chess = abs(Mxy_chess(:, idx_centro));   % (Nf x 1)

% Cuánto excita CHESS en cada banda
sup_chess_grasa = mean(perfil_chess(fatval));
sup_chess_agua  = mean(perfil_chess(watval));

fprintf('CHESS excitación banda grasa: %.4f\n', sup_chess_grasa);
fprintf('CHESS excitación banda agua:  %.4f (ideal ~0)\n', sup_chess_agua);

% Verificar que CHESS es uniforme en X (todas las columnas iguales)
fprintf('Variación espacial CHESS (debe ser ~0): %.6f\n', ...
    std(abs(Mxy_chess(fatval(1), :))));

%% BLOQUE 3 - Aplicar CHESS al phantom
% Grasa:    señal × (1 - sup_grasa) -> suprimida
% Interior: señal × (1 - sup_agua)  -> casi sin cambio
% Fondo:    señal = 0
factor_chess = ones(128,128);
factor_chess(roi_grasa)    = 1 - sup_chess_grasa;
factor_chess(roi_interior) = 1 - sup_chess_agua;
factor_chess(roi_fondo)    = 0;

imagen_chess = factor_chess .* phantom;

%% BLOQUE 4 - Figura resultado CHESS
figure('Name','Resultado CHESS en phantom', 'Position', [50 50 1200 420]);

subplot(1,3,1);
imagesc(phantom); colormap(gca,'gray'); axis image; colorbar;
title('Sin supresión');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,2);
imagesc(imagen_chess); colormap(gca,'gray'); axis image; colorbar;
title(sprintf('Con CHESS\nSup. grasa: %.1f%%', ...
    100*(1 - mean(imagen_chess(roi_grasa)) / mean(phantom(roi_grasa)))));
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,3);
imagesc(phantom - imagen_chess); colormap(gca,'hot'); axis image; colorbar;
title('Señal suprimida por CHESS');
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
hold off;
% Debe brillar solo en el anillo rojo

%% BLOQUE 5 - Perfil espectral de CHESS
figure('Name','Perfil espectral CHESS');
plot(frange, perfil_chess, 'b-', 'LineWidth', 2);
hold on;
xline(fatfreq, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia grasa');
xline(0,       'g--', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia agua');
hold off;
xlabel('Frecuencia [Hz]'); ylabel('|M_{xy}|');
title('Perfil espectral CHESS (punto central)');
legend; grid on;
% Debe tener pico alto en ~-220 Hz y bajo en 0 Hz

%% BLOQUE 6 - Cuantificación
s_orig_grasa  = mean(phantom(roi_grasa));
s_chess_grasa = mean(imagen_chess(roi_grasa));
s_orig_agua   = mean(phantom(roi_interior));
s_chess_agua  = mean(imagen_chess(roi_interior));

fprintf('\n=== Resultado CHESS ===\n');
fprintf('%-20s %-12s %-15s\n', 'ROI', 'Señal orig', 'Señal CHESS');
fprintf('%-20s %-12.4f %-15.4f  ->  sup=%.1f%%\n', ...
    'Grasa (capa 3)', s_orig_grasa, s_chess_grasa, ...
    100*(1-s_chess_grasa/s_orig_grasa));
fprintf('%-20s %-12.4f %-15.4f  ->  cambio=%.1f%%\n', ...
    'Agua (capa 2)', s_orig_agua, s_chess_agua, ...
    100*(1-s_chess_agua/s_orig_agua));
fprintf('\n-> Ideal: supresión alta en grasa + cambio ~0%% en agua\n');