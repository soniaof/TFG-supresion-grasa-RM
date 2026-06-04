% chess_phantom_B0.m
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

%% BLOQUE 0.5 - Mapa de B0 sintético (1.5T)
% Se simula una variación local de frecuencia en el rango aproximado
% de -100 a 100 Hz. La componente positiva representa una variación
% suave del campo dentro del tejido, mientras que la perturbación negativa
% simula el efecto de una región con gas intestinal cercana al borde izquierdo.

[Xph, Yph] = meshgrid(1:128, 1:128);

% --- Tendencia global suave positiva ---
% Variación de fondo del campo, con valores positivos de hasta ~100 Hz
cx_pos = 95;
cy_pos = 64;
sigma_pos_x = 25;
sigma_pos_y = 40;

B0_global = 100 * exp(-((Xph-cx_pos).^2/(2*sigma_pos_x^2) + ...
                        (Yph-cy_pos).^2/(2*sigma_pos_y^2)));

% --- Perturbación negativa izquierda: simula gas intestinal ---
% Desplazamiento local negativo de hasta aproximadamente -100 Hz
% Se coloca cerca del borde izquierdo para afectar al anillo de grasa cercano.
cx_gas = 32;
cy_gas = 60;
sigma_gas = 18;

dist_gas = sqrt((Xph-cx_gas).^2 + (Yph-cy_gas).^2);
B0_gas = -100 * exp(-dist_gas.^2 / (2*sigma_gas^2));

% --- Mapa total B0 ---
B0map = B0_global + B0_gas;

% --- Enmascarar: solo tejido (grasa + interior), fuera del phantom = 0 ---
B0map = B0map .* double(roi_grasa | roi_interior);

% --- Limitar rango por seguridad ---
B0map(B0map > 100)  = 100;
B0map(B0map < -100) = -100;

% --- Visualizar mapa de B0 ---
figure('Name','Mapa de B0 sintético 1.5T - CHESS');
imagesc(B0map); colormap(gca,'jet'); axis image;
clim([-100 100]);
h = colorbar; ylabel(h, 'B_0 Field Map (Hz)');
title('Mapa de B0 sintético 1.5T [Hz] - CHESS');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
hold off;

%% BLOQUE 1 - Parámetros espectrales
fatfreq_pos = Sys.gamma * Sys.field * Sys.fat * 1e-6;
fatfreq     = -abs(fatfreq_pos);   % ~-223.5 Hz a 1.5T
freqband    = 15;
fatval = find(abs(frange - fatfreq) < freqband);
watval = find(abs(frange - 0)       < freqband);

fprintf('\nFrecuencia de grasa nominal: %.1f Hz\n', fatfreq);
fprintf('fatval: %d puntos\n', numel(fatval));

%% BLOQUE 2 - Respuesta espectral de CHESS
% Mxyspec viene de test_fat_sat.m, tamaño (Nx x Nf) = (128 x 41)
% Transponer a (Nf x Nx) = (41 x 128)
Mxy_chess = Mxyspec.';

% CHESS no tiene gradientes espaciales -> perfil igual en todas las X
[~, idx_centro] = min(abs(x));
perfil_chess = abs(Mxy_chess(:, idx_centro));   % (Nf x 1)

sup_chess_grasa = mean(perfil_chess(fatval));
sup_chess_agua  = mean(perfil_chess(watval));

fprintf('CHESS excitación banda grasa: %.4f\n', sup_chess_grasa);
fprintf('CHESS excitación banda agua:  %.4f (ideal ~0)\n', sup_chess_agua);

% Visualizar perfil espectral
figure('Name','Perfil espectral CHESS');
plot(frange, perfil_chess, 'b-', 'LineWidth', 2);
hold on;
xline(fatfreq, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia grasa nominal');
xline(0,       'g--', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia agua');
hold off;
xlabel('Frecuencia [Hz]'); ylabel('|M_{xy}|');
title('Perfil espectral CHESS (punto central)');
legend; grid on;

%% BLOQUE 3 - Aplicar CHESS al phantom con mapa de B0
% CHESS usa el mismo perfil espectral para todos los píxeles, ya que no tiene
% selectividad espacial. Con B0 inhomogéneo, cada píxel tiene su frecuencia
% local desplazada:
%
%   f_grasa_local(x,y) = fatfreq + B0map(y,x)
%   f_agua_local(x,y)  = 0       + B0map(y,x)

factor_chess = ones(128,128);

for ix = 1:128
    for iy = 1:128

        % Frecuencias locales según B0 en este píxel
        f_grasa_local = fatfreq + B0map(iy, ix);
        f_agua_local  = 0       + B0map(iy, ix);

        % Índice espectral más cercano en frange
        [~, if_grasa] = min(abs(frange - f_grasa_local));
        [~, if_agua]  = min(abs(frange - f_agua_local));

        % Supresión local: mismo perfil para todos los X
        sup_grasa_local = perfil_chess(if_grasa);
        sup_agua_local  = perfil_chess(if_agua);

        if roi_grasa(iy, ix)
            factor_chess(iy, ix) = 1 - sup_grasa_local;
        elseif roi_interior(iy, ix)
            factor_chess(iy, ix) = 1 - sup_agua_local;
        else
            factor_chess(iy, ix) = 0; % fondo
        end
    end
end

% Evitar valores fuera de rango por seguridad numérica
factor_chess = max(0, min(1, factor_chess));

imagen_chess = factor_chess .* phantom;

%% BLOQUE 4 - Figura resultado CHESS con B0
figure('Name','Resultado CHESS con B0', 'Position', [50 50 1200 420]);

subplot(1,3,1);
imagesc(phantom); colormap(gca,'gray'); axis image; colorbar;
title('Sin supresión');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,2);
imagesc(imagen_chess); colormap(gca,'gray'); axis image; colorbar;
title(sprintf('Con CHESS + B0\nSup. grasa: %.1f%%', ...
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

%% BLOQUE 5 - Mapa del factor de supresión
figure('Name','Mapa factor supresión CHESS con B0');
imagesc(factor_chess); colormap(gca,'parula'); axis image; colorbar;
title('Factor de supresión CHESS por píxel (1=sin cambio, 0=suprimido)');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

%% BLOQUE 6 - B0map sobre phantom
figure('Name','B0map sobre phantom - CHESS');
imagesc(B0map .* double(roi_grasa | roi_interior));
colormap(gca,'jet'); axis image;
clim([-100 100]);
h = colorbar; ylabel(h,'B_0 [Hz]');
title('B0 local [Hz] dentro del phantom');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
contour(dist_gas <= sigma_gas, 1, 'y--', 'LineWidth', 2);
hold off;

%% BLOQUE 7 - Comparación: supresión en zona afectada por B0 vs zona normal
% Para evitar medias vacías, se define la zona afectada como la grasa situada
% en las regiones con B0 más negativo. La zona normal se define como grasa
% con |B0| bajo.

b0_grasa_vals = B0map(roi_grasa);

% Zona gas: píxeles de grasa con B0 negativo importante
umbral_gas = prctile(b0_grasa_vals, 20);
mascara_grasa_gas = roi_grasa & (B0map <= umbral_gas);

% Zona normal: píxeles de grasa con menor |B0|
umbral_normal = prctile(abs(b0_grasa_vals), 30);
mascara_grasa_normal = roi_grasa & (abs(B0map) <= umbral_normal);

% Si alguna máscara queda vacía, se avisa y se evita NaN
if ~any(mascara_grasa_gas(:))
    warning('La máscara de grasa en zona de gas está vacía.');
    sup_zona_gas = NaN;
else
    sup_zona_gas = 1 - mean(factor_chess(mascara_grasa_gas));
end

if ~any(mascara_grasa_normal(:))
    warning('La máscara de grasa en zona normal está vacía.');
    sup_zona_normal = NaN;
else
    sup_zona_normal = 1 - mean(factor_chess(mascara_grasa_normal));
end

fprintf('\n=== Efecto del B0 en la supresión de grasa (CHESS) ===\n');
fprintf('Supresión grasa zona afectada por B0: %.1f%%\n', sup_zona_gas*100);
fprintf('Supresión grasa zona normal:          %.1f%%\n', sup_zona_normal*100);
fprintf('-> Diferencia: %.1f pp\n', (sup_zona_normal - sup_zona_gas)*100);

%% BLOQUE 8 - Cuantificación global
s_orig_grasa  = mean(phantom(roi_grasa));
s_chess_grasa = mean(imagen_chess(roi_grasa));
s_orig_agua   = mean(phantom(roi_interior));
s_chess_agua  = mean(imagen_chess(roi_interior));

fprintf('\n=== Resultado global CHESS con B0 (1.5T) ===\n');
fprintf('%-20s %-12s %-15s\n', 'ROI', 'Señal orig', 'Señal CHESS');
fprintf('%-20s %-12.4f %-15.4f  ->  sup=%.1f%%\n', ...
    'Grasa (capa 3)', s_orig_grasa, s_chess_grasa, ...
    100*(1-s_chess_grasa/s_orig_grasa));
fprintf('%-20s %-12.4f %-15.4f  ->  cambio=%.1f%%\n', ...
    'Agua (capa 2)', s_orig_agua, s_chess_agua, ...
    100*(1-s_chess_agua/s_orig_agua));
fprintf('\nNota: con B0 inhomogéneo, cada píxel presenta una frecuencia local distinta.\n');
fprintf('CHESS está centrado en la frecuencia nominal de grasa y no se adapta\n');
fprintf('a los desplazamientos locales producidos por B0.\n');

%% BLOQUE 9 - Comparación directa CHESS vs SPSP 
% Solo funciona si se ha ejecutado antes spsp_phantom_B0.m
if exist('factor_spsp', 'var')
    figure('Name','Comparación CHESS vs SPSP con B0', ...
        'Position', [50 50 1200 420]);

    subplot(1,3,1);
    imagesc(B0map .* double(roi_grasa | roi_interior));
    colormap(gca,'jet'); axis image; clim([-100 100]); colorbar;
    title('Mapa B0 [Hz]');
    xlabel('X [px]'); ylabel('Y [px]');

    subplot(1,3,2);
    imagesc(factor_chess); colormap(gca,'parula'); axis image;
    clim([0 1]); colorbar;
    title(sprintf('Factor supresión CHESS\nzona B0: %.1f%% | normal: %.1f%%', ...
        sup_zona_gas*100, sup_zona_normal*100));
    xlabel('X [px]'); ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
    hold off;

    subplot(1,3,3);
    imagesc(factor_spsp); colormap(gca,'parula'); axis image;
    clim([0 1]); colorbar;
    title('Factor supresión SPSP');
    xlabel('X [px]'); ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
    hold off;

    sgtitle('Comparación CHESS vs SPSP con B0 inhomogéneo');
end

%% 
%% BLOQUE FINAL - Comparación CHESS vs SPSP con B0

if exist('factor_spsp', 'var') && exist('factor_chess', 'var') && ...
   exist('imagen_spsp', 'var') && exist('imagen_chess', 'var')

    % --- Cálculo de porcentajes globales ---
    s_orig_grasa  = mean(phantom(roi_grasa));
    s_chess_grasa = mean(imagen_chess(roi_grasa));
    s_spsp_grasa  = mean(imagen_spsp(roi_grasa));

    s_orig_agua   = mean(phantom(roi_interior));
    s_chess_agua  = mean(imagen_chess(roi_interior));
    s_spsp_agua   = mean(imagen_spsp(roi_interior));

    sup_chess_grasa = 100 * (1 - s_chess_grasa / s_orig_grasa);
    sup_spsp_grasa  = 100 * (1 - s_spsp_grasa  / s_orig_grasa);

    cambio_chess_agua = 100 * (1 - s_chess_agua / s_orig_agua);
    cambio_spsp_agua  = 100 * (1 - s_spsp_agua  / s_orig_agua);

    % --- Figura comparativa final ---
    figure('Name','Comparación final CHESS vs SPSP con B0', ...
           'Position', [50 50 1400 650]);

    % Mapa B0
    subplot(2,3,1);
    imagesc(B0map .* double(roi_grasa | roi_interior));
    colormap(gca,'jet');
    axis image;
    clim([-100 100]);
    colorbar;
    title('Mapa B_0 [Hz]');
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
    hold off;

    % Imagen original
    subplot(2,3,4);
    imagesc(phantom);
    colormap(gca,'gray');
    axis image;
    colorbar;
    title('Imagen original');
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
    hold off;

    % CHESS con B0
    subplot(2,3,2);
    imagesc(imagen_chess);
    colormap(gca,'gray');
    axis image;
    colorbar;
    title(sprintf('CHESS + B_0\nSup. grasa = %.1f%% | cambio agua = %.1f%%', ...
        sup_chess_grasa, cambio_chess_agua));
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
    hold off;

    % Señal suprimida por CHESS
    subplot(2,3,5);
    imagesc(phantom - imagen_chess);
    colormap(gca,'hot');
    axis image;
    colorbar;
    title('Señal suprimida por CHESS');
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
    hold off;

    % SPSP con B0
    subplot(2,3,3);
    imagesc(imagen_spsp);
    colormap(gca,'gray');
    axis image;
    colorbar;
    title(sprintf('SPSP + B_0\nSup. grasa = %.1f%% | cambio agua = %.1f%%', ...
        sup_spsp_grasa, cambio_spsp_agua));
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
    hold off;

    % Señal suprimida por SPSP
    subplot(2,3,6);
    imagesc(phantom - imagen_spsp);
    colormap(gca,'hot');
    axis image;
    colorbar;
    title('Señal suprimida por SPSP');
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
    hold off;

    sgtitle('Comparación de CHESS y SPSP con inhomogeneidad de B_0');

    % Guardar figura para el TFG
    exportgraphics(gcf, 'comparacion_CHESS_SPSP_B0.png', 'Resolution', 300);

    % Mostrar resultados 
    fprintf('\n=== Comparación final CHESS vs SPSP con B0 ===\n');
    fprintf('CHESS: supresión grasa = %.1f%% | cambio agua = %.1f%%\n', ...
        sup_chess_grasa, cambio_chess_agua);
    fprintf('SPSP:  supresión grasa = %.1f%% | cambio agua = %.1f%%\n', ...
        sup_spsp_grasa, cambio_spsp_agua);

else
    warning('Faltan variables. Ejecuta antes los scripts de CHESS y SPSP con B0 en el mismo workspace.');
end

%% %% FIGURA FINAL COMPACTA - Comparación B0, CHESS y SPSP

if exist('B0map','var') && exist('imagen_chess','var') && exist('imagen_spsp','var')

    % --- Cálculo de métricas globales ---
    s_orig_grasa  = mean(phantom(roi_grasa));
    s_chess_grasa = mean(imagen_chess(roi_grasa));
    s_spsp_grasa  = mean(imagen_spsp(roi_grasa));

    sup_chess_grasa = 100 * (1 - s_chess_grasa / s_orig_grasa);
    sup_spsp_grasa  = 100 * (1 - s_spsp_grasa  / s_orig_grasa);

    % Si se quiere incluir también el cambio en agua:
    s_orig_agua   = mean(phantom(roi_interior));
    s_chess_agua  = mean(imagen_chess(roi_interior));
    s_spsp_agua   = mean(imagen_spsp(roi_interior));

    cambio_chess_agua = 100 * (1 - s_chess_agua / s_orig_agua);
    cambio_spsp_agua  = 100 * (1 - s_spsp_agua  / s_orig_agua);

    % --- Figura ---
    figure('Name','Comparación final B0 - CHESS - SPSP', ...
           'Position', [100 100 1350 420]);

    % Panel 1: mapa B0
    subplot(1,3,1);
    imagesc(B0map .* double(roi_grasa | roi_interior));
    colormap(gca,'jet');
    axis image;
    clim([-100 100]);
    colorbar;
    title('Mapa de B_0 [Hz]');
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.2);
    contour(dist_gas <= sigma_gas, 1, 'y--', 'LineWidth', 1.8);
    hold off;

    % Panel 2: CHESS
    subplot(1,3,2);
    imagesc(imagen_chess);
    colormap(gca,'gray');
    axis image;
    colorbar;
    title(sprintf('CHESS + B_0\nSup. grasa = %.1f%%', sup_chess_grasa));
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.2);
    hold off;

    % Panel 3: SPSP
    subplot(1,3,3);
    imagesc(imagen_spsp);
    colormap(gca,'gray');
    axis image;
    colorbar;
    title(sprintf('SPSP + B_0\nSup. grasa = %.1f%%', sup_spsp_grasa));
    xlabel('X [px]');
    ylabel('Y [px]');
    hold on;
    contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.2);
    hold off;

    sgtitle('Comparación de CHESS y SPSP en presencia de inhomogeneidades de B_0');

    exportgraphics(gcf, 'comparacion_B0_CHESS_SPSP_compacta.png', 'Resolution', 300);

    fprintf('\nFigura compacta guardada como: comparacion_B0_CHESS_SPSP_compacta.png\n');
    fprintf('CHESS: supresión grasa = %.1f%% | cambio agua = %.1f%%\n', ...
        sup_chess_grasa, cambio_chess_agua);
    fprintf('SPSP:  supresión grasa = %.1f%% | cambio agua = %.1f%%\n', ...
        sup_spsp_grasa, cambio_spsp_agua);

else
    warning('Faltan variables. Ejecuta antes los scripts de CHESS y SPSP con B0.');
end
