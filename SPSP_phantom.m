% spsp_phantom.m
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

%% BLOQUE 1 - Parámetros
fatfreq  = -(Sys.fat * 1e-6 * Sys.gamma * Sys.field);
freqband = 15;
fatval   = find(abs(frange - fatfreq) < freqband);
watval   = find(abs(frange - 0)       < freqband);

fprintf('\nFrecuencia de grasa: %.1f Hz\n', fatfreq);

%% BLOQUE 2 - Respuesta espectral de SPSP
% mopt viene de test_2D_spec.m, tamaño (Nf x Nx) = (41 x 128)
% A diferencia de CHESS, cada columna ix tiene un perfil espectral distinto
% porque SPSP sí tiene selectividad espacial

% Supresión en banda de grasa para cada posición X
sup_spsp_grasa_por_x = mean(abs(mopt(fatval, :)), 1);   % (1 x 128)
sup_spsp_agua_por_x  = mean(abs(mopt(watval, :)), 1);   % (1 x 128)

% Visualizar cómo varía la supresión a lo largo de X
% SPSP debería tener supresión alta en las X de grasa y baja en las de agua

figure('Name','Supresión SPSP a lo largo de X');
plot(x, sup_spsp_grasa_por_x, 'r-', 'LineWidth', 2, 'DisplayName', 'Banda grasa');
hold on;
plot(x, sup_spsp_agua_por_x,  'b-', 'LineWidth', 2, 'DisplayName', 'Banda agua');

xline(x(1),   'k--', 'HandleVisibility','off');
xline(x(10),  'k--', 'HandleVisibility','off');
xline(x(35),  'k--', 'HandleVisibility','off');
xline(x(75),  'k--', 'HandleVisibility','off');
xline(x(99),  'k--', 'HandleVisibility','off');
xline(x(118), 'k--', 'HandleVisibility','off');

hold off;
xlabel('X [cm]');
ylabel('Supresión media');
title('Supresión SPSP según posición X');
legend;
grid on;

%dfat = 1:10 → aprox -12.5 a -10.7 cm
%dfat = 35:75 → aprox -5.9 a 2.0 cm
%dfat = 99:118 → aprox 6.8 a 10.5 cm

%% BLOQUE 3 - Aplicar SPSP al phantom
% Para cada columna ix, aplicamos la supresión específica de esa posición X
% Esto es la diferencia clave con CHESS: el factor no es uniforme en X

factor_spsp = ones(128,128);

for ix = 1:128
    col_grasa    = roi_grasa(:, ix);
    col_interior = roi_interior(:, ix);
    col_fondo    = roi_fondo(:, ix);

    factor_spsp(col_grasa,    ix) = 1 - sup_spsp_grasa_por_x(ix);
    factor_spsp(col_interior, ix) = 1 - sup_spsp_agua_por_x(ix);
    factor_spsp(col_fondo,    ix) = 0;
end

imagen_spsp = factor_spsp .* phantom;

% %% BLOQUE 4 - Figura resultado SPSP
% figure('Name','Resultado SPSP en phantom', 'Position', [50 50 1200 420]);
% 
% subplot(1,3,1);
% imagesc(phantom); colormap(gca,'gray'); axis image; colorbar;
% title('Sin supresión');
% xlabel('X [px]'); ylabel('Y [px]');
% hold on;
% contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
% hold off;
% 
% subplot(1,3,2);
% imagesc(imagen_spsp); colormap(gca,'gray'); axis image; colorbar;
% title(sprintf('Con SPSP\nSup. grasa: %.1f%%', ...
%     100*(1 - mean(imagen_spsp(roi_grasa)) / mean(phantom(roi_grasa)))));
% xlabel('X [px]');
% hold on;
% contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
% hold off;
% 
% subplot(1,3,3);
% imagesc(phantom - imagen_spsp); colormap(gca,'hot'); axis image; colorbar;
% title('Señal suprimida por SPSP');
% xlabel('X [px]');
% hold on;
% contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
% hold off;

%% BLOQUE 4 - Figura resultado SPSP
figure('Name','Resultado SPSP en phantom', 'Position', [50 50 1200 420]);
subplot(1,3,1);
imagesc(phantom); colormap(gca,'gray'); axis image; colorbar;
title('Sin supresión');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,2);
imagesc(imagen_spsp); colormap(gca,'gray'); axis image; colorbar;
title(sprintf('Con SPSP\nSup. grasa: %.1f%%', ...
    100*(1 - mean(imagen_spsp(roi_grasa)) / mean(phantom(roi_grasa)))));
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;

subplot(1,3,3);
diff_img = phantom - imagen_spsp;
imagesc(diff_img); colormap(gca,'turbo'); axis image; colorbar;

% Ajuste de rango: mínimo en 0, máximo en el percentil 99 para no saturar
valid_vals = diff_img(diff_img > 0);
if ~isempty(valid_vals)
    clim([0, prctile(valid_vals, 99)]);
end

title('Señal suprimida por SPSP');
xlabel('X [px]');
hold on;
contour(double(roi_grasa), 1, 'w-', 'LineWidth', 1.5);
hold off;


%% BLOQUE 5 - Mapa del factor de supresión
% Muestra visualmente cómo varía la supresión en el espacio 2D
% Esto NO existe en CHESS porque allí el factor es uniforme
figure('Name','Mapa factor supresión SPSP');
imagesc(factor_spsp); colormap(gca,'parula'); axis image; colorbar;
title('Factor de supresión SPSP por píxel (1=sin cambio, 0=suprimido)');
xlabel('X [px]'); ylabel('Y [px]');
hold on;
contour(double(roi_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;
% el factor varía en X porque SPSP tiene selectividad espacial

%% BLOQUE 6 - Cuantificación
s_orig_grasa  = mean(phantom(roi_grasa));
s_spsp_grasa  = mean(imagen_spsp(roi_grasa));
s_orig_agua   = mean(phantom(roi_interior));
s_spsp_agua   = mean(imagen_spsp(roi_interior));

fprintf('\n=== Resultado SPSP ===\n');
fprintf('%-20s %-12s %-15s\n', 'ROI', 'Señal orig', 'Señal SPSP');
fprintf('%-20s %-12.4f %-15.4f  ->  sup=%.1f%%\n', ...
    'Grasa (capa 3)', s_orig_grasa, s_spsp_grasa, ...
    100*(1-s_spsp_grasa/s_orig_grasa));
fprintf('%-20s %-12.4f %-15.4f  ->  cambio=%.1f%%\n', ...
    'Agua (capa 2)', s_orig_agua, s_spsp_agua, ...
    100*(1-s_spsp_agua/s_orig_agua));
fprintf('\n-> Ideal: supresión alta en grasa + cambio ~0%% en agua\n');