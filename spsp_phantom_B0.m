% spsp_phantom_B0.m
%
% ANTES haber ejecutado:
%   test_2D_spec_espiral.m -> genera: mopt, frange, x, Sys, T1, T2

%% BLOQUE 0 - Cargar phantom y definir ROIs
load('/Users/sonita/Desktop/4º/TFG/Datos TFG/phantomdigital.mat');

roi_grasa    = logical(tejidos(:,:,3));   % capa 3 = anillo = GRASA
roi_interior = logical(tejidos(:,:,2));   % capa 2 = interior = AGUA
roi_fondo    = logical(tejidos(:,:,1));   % capa 1 = fondo/aire

figure('Name','ROIs del phantom - SPSP');
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
%
%este bloque es el mismo que en chess_phantom_B0.m para que
% CHESS y SPSP se comparen bajo las mismas condiciones de B0.

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

% --- Enmascarar: fuera del phantom = 0 ---
B0map = B0map .* double(roi_grasa | roi_interior);

% --- Limitar rango por seguridad ---
B0map(B0map > 100)  = 100;
B0map(B0map < -100) = -100;

% --- Visualizar mapa de B0 ---
figure('Name','Mapa de B0 sintético 1.5T - SPSP');
imagesc(B0map); colormap(gca,'jet'); axis image;
clim([-100 100]);
h = colorbar; ylabel(h, 'B_0 Field Map (Hz)');
title('Mapa de B0 sintético 1.5T [Hz] - SPSP');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
hold off;

%% BLOQUE 1 - Parámetros espectrales
fatfreq  = -(Sys.fat * 1e-6 * Sys.gamma * Sys.field);   % ~-223.5 Hz a 1.5T
freqband = 15;
fatval   = find(abs(frange - fatfreq) < freqband);
watval   = find(abs(frange - 0)       < freqband);

fprintf('\nFrecuencia de grasa nominal: %.1f Hz\n', fatfreq);

%% BLOQUE 2 - Respuesta espectral del pulso SPSP
% mopt es (Nf x 128): perfil espectral para cada columna X del phantom

sup_spsp_grasa_por_x = mean(abs(mopt(fatval, :)), 1);   % (1 x 128)
sup_spsp_agua_por_x  = mean(abs(mopt(watval, :)), 1);   % (1 x 128)

figure('Name','Supresión SPSP a lo largo de X');
plot(1:128, sup_spsp_grasa_por_x, 'r-', 'LineWidth', 2, 'DisplayName', 'Banda grasa');
hold on;
plot(1:128, sup_spsp_agua_por_x,  'b-', 'LineWidth', 2, 'DisplayName', 'Banda agua');
hold off;
xlabel('X [px]');
ylabel('Supresión media');
title('Supresión SPSP según posición X');
legend; grid on;

%% BLOQUE 3 - Aplicar SPSP al phantom con mapa de B0
% Cada píxel tiene su propia frecuencia local:
%
%   f_grasa_local(x,y) = fatfreq + B0map(y,x)
%   f_agua_local(x,y)  = 0       + B0map(y,x)
%
% Si B0 desplaza la frecuencia fuera de la banda efectiva del pulso SPSP,
% la supresión empeora en esa región.

factor_spsp = ones(128,128);

for ix = 1:128
    for iy = 1:128

        % Frecuencias locales según B0 en este píxel
        f_grasa_local = fatfreq + B0map(iy, ix);
        f_agua_local  = 0       + B0map(iy, ix);

        % Índice espectral más cercano en frange
        [~, if_grasa] = min(abs(frange - f_grasa_local));
        [~, if_agua]  = min(abs(frange - f_agua_local));

        % Supresión local del pulso SPSP en esa columna X y esa frecuencia
        sup_grasa_local = abs(mopt(if_grasa, ix));
        sup_agua_local  = abs(mopt(if_agua,  ix));

        if roi_grasa(iy, ix)
            factor_spsp(iy, ix) = 1 - sup_grasa_local;
        elseif roi_interior(iy, ix)
            factor_spsp(iy, ix) = 1 - sup_agua_local;
        else
            factor_spsp(iy, ix) = 0; % fondo = negro
        end
    end
end

% Evitar valores fuera de rango por seguridad numérica
factor_spsp = max(0, min(1, factor_spsp));

imagen_spsp = factor_spsp .* phantom;

%% BLOQUE 4 - Figura resultado SPSP con B0
figure('Name','Resultado SPSP con B0', 'Position', [50 50 1200 420]);

subplot(1,3,1);
imagesc(phantom); colormap(gca,'gray'); axis image; colorbar;
title('Sin supresión');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,2);
imagesc(imagen_spsp); colormap(gca,'gray'); axis image; colorbar;
title(sprintf('Con SPSP + B0\nSup. grasa: %.1f%%', ...
    100*(1 - mean(imagen_spsp(roi_grasa)) / mean(phantom(roi_grasa)))));
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,3);
imagesc(phantom - imagen_spsp); colormap(gca,'hot'); axis image; colorbar;
title('Señal suprimida por SPSP');
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
hold off;

%% BLOQUE 5 - Mapa del factor de supresión
figure('Name','Mapa factor supresión SPSP con B0');
imagesc(factor_spsp); colormap(gca,'parula'); axis image; colorbar;
title('Factor de supresión SPSP por píxel (1=sin cambio, 0=suprimido)');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

%% BLOQUE 6 - B0map superpuesto al phantom
figure('Name','B0map sobre phantom - SPSP');
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
    sup_zona_gas = 1 - mean(factor_spsp(mascara_grasa_gas));
end

if ~any(mascara_grasa_normal(:))
    warning('La máscara de grasa en zona normal está vacía.');
    sup_zona_normal = NaN;
else
    sup_zona_normal = 1 - mean(factor_spsp(mascara_grasa_normal));
end

fprintf('\n=== Efecto del B0 en la supresión de grasa (SPSP) ===\n');
fprintf('Supresión grasa zona afectada por B0: %.1f%%\n', sup_zona_gas*100);
fprintf('Supresión grasa zona normal:          %.1f%%\n', sup_zona_normal*100);
fprintf('-> Diferencia: %.1f pp\n', (sup_zona_normal - sup_zona_gas)*100);

%% BLOQUE 8 - Cuantificación global
s_orig_grasa = mean(phantom(roi_grasa));
s_spsp_grasa = mean(imagen_spsp(roi_grasa));
s_orig_agua  = mean(phantom(roi_interior));
s_spsp_agua  = mean(imagen_spsp(roi_interior));

fprintf('\n=== Resultado global SPSP con B0 (1.5T) ===\n');
fprintf('%-20s %-12s %-15s\n', 'ROI', 'Señal orig', 'Señal SPSP');
fprintf('%-20s %-12.4f %-15.4f  ->  sup=%.1f%%\n', ...
    'Grasa (capa 3)', s_orig_grasa, s_spsp_grasa, ...
    100*(1-s_spsp_grasa/s_orig_grasa));
fprintf('%-20s %-12.4f %-15.4f  ->  cambio=%.1f%%\n', ...
    'Agua (capa 2)', s_orig_agua, s_spsp_agua, ...
    100*(1-s_spsp_agua/s_orig_agua));
fprintf('\nNota: con B0 inhomogéneo, cada píxel presenta una frecuencia local distinta.\n');
fprintf('En SPSP, el efecto depende tanto del desplazamiento de frecuencia como\n');
fprintf('de la posición espacial considerada en el diseño del pulso.\n');