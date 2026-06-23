clear;
clc;
close all;

%% ============================================================
% RESUMEN Y VISUALIZACIÓN DE BIOMARCADORES DE LA ONDA T
%
% Este script carga los resultados del análisis morfológico de la onda T
% y genera tablas resumen, figuras de onda T media, mapas de correlación
% y diagramas de caja de los biomarcadores calculados.
%
% Entrada:
%   - ondasT_features_completas.mat
%
% Salidas:
%   - features_ondaT.xlsx
%   - tablas_correlaciones_ondaT.xlsx
%   - correlaciones_ondaT.mat
%   - figuras de onda T media por grupo
%   - mapas de correlación
%   - boxplots de biomarcadores y correlaciones
%
% Grupos analizados:
%   - SANO
%   - FA_PAROXISTICA_RS
%
% Requisitos:
%   - Resultados generados previamente por el script de análisis de onda T.
%% ============================================================

%% CONFIGURACIÓN GENERAL

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\analisis_ondaTHIBRIDOnuevo';
archivo_mat = fullfile(carpeta_out, 'ondasT_features_completas.mat');

carpeta_boxplots = fullfile(carpeta_out, 'boxplots_correlaciones_T');
if ~exist(carpeta_boxplots, 'dir')
    mkdir(carpeta_boxplots);
end

carpeta_heatmaps = fullfile(carpeta_out, 'heatmaps_correlacion_T');
if ~exist(carpeta_heatmaps, 'dir')
    mkdir(carpeta_heatmaps);
end

carpeta_pmedia = fullfile(carpeta_out, 'figuras_T_media');
if ~exist(carpeta_pmedia, 'dir')
    mkdir(carpeta_pmedia);
end

if ~exist(archivo_mat, 'file')
    error('No se encuentra el archivo %s. Ejecuta primero el script principal de onda T.', archivo_mat);
end

load(archivo_mat, ...
    'tabla_features', ...
    'resultados_sanos', ...
    'resultados_parox', ...
    'correlaciones_intra_sanos', ...
    'correlaciones_intra_parox');

%% ETIQUETAS ABREVIADAS

grupo_sano = 'SANO';
grupo_parox = 'FA_PA_RS';
grupo_inter = 'SANO vs FA_PA_RS';

%% FIGURA DE ONDA T MEDIA POR GRUPO

if ~isempty(resultados_sanos) || ~isempty(resultados_parox)

    if ~isempty(resultados_sanos)
        n_muestras = numel(resultados_sanos(1).T_media);
    else
        n_muestras = numel(resultados_parox(1).T_media);
    end

    Fs = 500;
    t_ms = ((1:n_muestras) - ceil(n_muestras/2)) / Fs * 1000;

    Tmed_sanos = [];
    for i = 1:numel(resultados_sanos)
        if isfield(resultados_sanos(i), 'T_media') && ~isempty(resultados_sanos(i).T_media)
            Tmed_sanos = [Tmed_sanos; resultados_sanos(i).T_media(:)']; %#ok<AGROW>
        end
    end

    Tmed_parox = [];
    for i = 1:numel(resultados_parox)
        if isfield(resultados_parox(i), 'T_media') && ~isempty(resultados_parox(i).T_media)
            Tmed_parox = [Tmed_parox; resultados_parox(i).T_media(:)']; %#ok<AGROW>
        end
    end

    %% Escala común para ambos grupos

    valores_y = [Tmed_sanos(:); Tmed_parox(:)];
    valores_y = valores_y(isfinite(valores_y));

    if ~isempty(valores_y)
        ymin = min(valores_y);
        ymax = max(valores_y);
        margen = 0.05 * (ymax - ymin);

        if margen == 0
            margen = 0.01;
        end

        ylim_comun = [ymin - margen, ymax + margen];
    else
        ylim_comun = [];
    end

    %% Figura con misma escala en ambos subplots

    f = figure('Visible','off');

    subplot(1,2,1)
    hold on

    if ~isempty(Tmed_sanos)
        plot(t_ms, Tmed_sanos', ...
            'Color', [0.75 0.75 0.75], ...
            'LineWidth', 0.6);

        plot(t_ms, mean(Tmed_sanos, 1, 'omitnan'), ...
            'k', ...
            'LineWidth', 2);
    end

    hold off
    title('T media - SANO', 'Interpreter','none');
    ylabel('Amplitud (mV)', 'Interpreter','none', 'FontSize', 14);
    xlabel('Tiempo relativo (ms)', 'Interpreter','none');
    grid on
    box on

    if ~isempty(ylim_comun)
        ylim(ylim_comun);
    end

    subplot(1,2,2)
    hold on

    if ~isempty(Tmed_parox)
        plot(t_ms, Tmed_parox', ...
            'Color', [0.75 0.75 0.75], ...
            'LineWidth', 0.6);

        plot(t_ms, mean(Tmed_parox, 1, 'omitnan'), ...
            'k', ...
            'LineWidth', 2);
    end

    hold off
    title('T media - FA_PA_RS', 'Interpreter','none');
    ylabel('Amplitud (mV)', 'Interpreter','none', 'FontSize', 14);
    xlabel('Tiempo relativo (ms)', 'Interpreter','none');
    grid on
    box on

    if ~isempty(ylim_comun)
        ylim(ylim_comun);
    end

    sgtitle('Ondas T medias por paciente en cada grupo', 'Interpreter','none');

    saveas(f, fullfile(carpeta_pmedia, 'T_media_por_grupo.png'));
    close(f)
end

%% 1) CORRELACIONES INTRAGRUPO

correlacion_media_sanos = [];
correlacion_media_parox = [];
correlacion_std_sanos = [];
correlacion_std_parox = [];

% Sanos
for i = 1:numel(resultados_sanos)-1
    for j = i+1:numel(resultados_sanos)

        if isfield(resultados_sanos(i), 'T_media') && isfield(resultados_sanos(j), 'T_media')
            correlacion_media_sanos(end + 1, 1) = corr(resultados_sanos(i).T_media', resultados_sanos(j).T_media');
        end

        if isfield(resultados_sanos(i), 'T_std') && isfield(resultados_sanos(j), 'T_std')
            correlacion_std_sanos(end + 1, 1) = corr(resultados_sanos(i).T_std', resultados_sanos(j).T_std');
        end
    end
end

% Paroxísticos
for i = 1:numel(resultados_parox)-1
    for j = i+1:numel(resultados_parox)

        if isfield(resultados_parox(i), 'T_media') && isfield(resultados_parox(j), 'T_media')
            correlacion_media_parox(end + 1, 1) = corr(resultados_parox(i).T_media', resultados_parox(j).T_media');
        end

        if isfield(resultados_parox(i), 'T_std') && isfield(resultados_parox(j), 'T_std')
            correlacion_std_parox(end + 1, 1) = corr(resultados_parox(i).T_std', resultados_parox(j).T_std');
        end
    end
end

%% 2) CORRELACIONES INTERGRUPOS

correlacion_cruzada_media = [];
correlacion_cruzada_std = [];

for i = 1:numel(resultados_sanos)
    for j = 1:numel(resultados_parox)

        if isfield(resultados_sanos(i), 'T_media') && isfield(resultados_parox(j), 'T_media')
            correlacion_cruzada_media(end + 1, 1) = corr(resultados_sanos(i).T_media', resultados_parox(j).T_media');
        end

        if isfield(resultados_sanos(i), 'T_std') && isfield(resultados_parox(j), 'T_std')
            correlacion_cruzada_std(end + 1, 1) = corr(resultados_sanos(i).T_std', resultados_parox(j).T_std');
        end
    end
end

%% MAPAS DE CORRELACIÓN

% --- SANO vs SANO (T_media)
M_sanos_media = nan(numel(resultados_sanos));

for i = 1:numel(resultados_sanos)
    for j = 1:numel(resultados_sanos)
        M_sanos_media(i,j) = corr(resultados_sanos(i).T_media', resultados_sanos(j).T_media');
    end
end

f = figure('Visible','off');
imagesc(M_sanos_media);
axis square
colorbar
colormap(jet)
caxis([-1 1])
title('Correlación entre pacientes SANO (T media)', 'Interpreter','none')
xlabel('Paciente', 'Interpreter','none')
ylabel('Paciente', 'Interpreter','none', 'FontSize', 14)
saveas(f, fullfile(carpeta_heatmaps, 'heatmap_sanos_media.png'));
close(f)

% --- FA_PA_RS vs FA_PA_RS (T_media)
M_parox_media = nan(numel(resultados_parox));

for i = 1:numel(resultados_parox)
    for j = 1:numel(resultados_parox)
        M_parox_media(i,j) = corr(resultados_parox(i).T_media', resultados_parox(j).T_media');
    end
end

f = figure('Visible','off');
imagesc(M_parox_media);
axis square
colorbar
colormap(jet)
caxis([-1 1])
title('Correlación entre FA_PA_RS (T media)', 'Interpreter','none')
xlabel('Paciente', 'Interpreter','none')
ylabel('Paciente', 'Interpreter','none', 'FontSize', 14)
saveas(f, fullfile(carpeta_heatmaps, 'heatmap_parox_media.png'));
close(f)

% --- SANO vs FA_PA_RS (T_media)
M_cruzada_media = nan(numel(resultados_sanos), numel(resultados_parox));

for i = 1:numel(resultados_sanos)
    for j = 1:numel(resultados_parox)
        M_cruzada_media(i,j) = corr(resultados_sanos(i).T_media', resultados_parox(j).T_media');
    end
end

f = figure('Visible','off');
imagesc(M_cruzada_media);
axis square
colorbar
colormap(jet)
caxis([-1 1])
title('Correlación SANO vs FA_PA_RS (T media)', 'Interpreter','none')
xlabel('FA_PA_RS', 'Interpreter','none')
ylabel('SANO', 'Interpreter','none', 'FontSize', 14)
saveas(f, fullfile(carpeta_heatmaps, 'heatmap_cruzada_media.png'));
close(f)

% --- SANO vs SANO (T_std)
M_sanos_std = nan(numel(resultados_sanos));

for i = 1:numel(resultados_sanos)
    for j = 1:numel(resultados_sanos)
        M_sanos_std(i,j) = corr(resultados_sanos(i).T_std', resultados_sanos(j).T_std');
    end
end

% --- FA_PA_RS vs FA_PA_RS (T_std)
M_parox_std = nan(numel(resultados_parox));

for i = 1:numel(resultados_parox)
    for j = 1:numel(resultados_parox)
        M_parox_std(i,j) = corr(resultados_parox(i).T_std', resultados_parox(j).T_std');
    end
end

% --- SANO vs FA_PA_RS (T_std)
M_cruzada_std = nan(numel(resultados_sanos), numel(resultados_parox));

for i = 1:numel(resultados_sanos)
    for j = 1:numel(resultados_parox)
        M_cruzada_std(i,j) = corr(resultados_sanos(i).T_std', resultados_parox(j).T_std');
    end
end


%% CORRELACIÓN MEDIA DE CADA PACIENTE CON EL RESTO DE SU GRUPO

% ============================================================
% T_media: un valor por paciente
% ============================================================

corr_grupo_media_sanos = nan(numel(resultados_sanos), 1);

for i = 1:numel(resultados_sanos)
    fila = M_sanos_media(i, :);
    fila(i) = NaN;  % excluir correlación consigo mismo
    corr_grupo_media_sanos(i) = mean(fila, 'omitnan');
end

corr_grupo_media_parox = nan(numel(resultados_parox), 1);

for i = 1:numel(resultados_parox)
    fila = M_parox_media(i, :);
    fila(i) = NaN;  % excluir correlación consigo mismo
    corr_grupo_media_parox(i) = mean(fila, 'omitnan');
end

% ============================================================
% T_std: un valor por paciente
% ============================================================

corr_grupo_std_sanos = nan(numel(resultados_sanos), 1);

for i = 1:numel(resultados_sanos)
    fila = M_sanos_std(i, :);
    fila(i) = NaN;  % excluir correlación consigo mismo
    corr_grupo_std_sanos(i) = mean(fila, 'omitnan');
end

corr_grupo_std_parox = nan(numel(resultados_parox), 1);

for i = 1:numel(resultados_parox)
    fila = M_parox_std(i, :);
    fila(i) = NaN;  % excluir correlación consigo mismo
    corr_grupo_std_parox(i) = mean(fila, 'omitnan');
end

% Tabla con un valor de similitud por paciente para T_media y T_std.
pacientes_sanos = string({resultados_sanos.registro})';
pacientes_parox = string({resultados_parox.registro})';

tabla_corr_grupo = table( ...
    [pacientes_sanos; pacientes_parox], ...
    [repmat("SANO", numel(resultados_sanos), 1); ...
     repmat("FA_PAROXISTICA_RS", numel(resultados_parox), 1)], ...
    [corr_grupo_media_sanos; corr_grupo_media_parox], ...
    [corr_grupo_std_sanos; corr_grupo_std_parox], ...
    'VariableNames', { ...
        'Paciente', ...
        'Grupo', ...
        'CorrGrupoMedia', ...
        'CorrGrupoStd'} ...
);

% Añadir ambos valores a tabla_features
tabla_features.CorrGrupoMedia = nan(height(tabla_features), 1);
tabla_features.CorrGrupoStd = nan(height(tabla_features), 1);

for i = 1:height(tabla_corr_grupo)

    idx = string(tabla_features.Paciente) == tabla_corr_grupo.Paciente(i) & ...
          string(tabla_features.Grupo) == tabla_corr_grupo.Grupo(i);

    tabla_features.CorrGrupoMedia(idx) = tabla_corr_grupo.CorrGrupoMedia(i);
    tabla_features.CorrGrupoStd(idx) = tabla_corr_grupo.CorrGrupoStd(i);
end

%% TABLA RESUMEN POR GRUPO

tabla_resumen_corr_grupo = table( ...
    {'SANO'; 'FA_PA_RS'}, ...
    [sum(isfinite(corr_grupo_media_sanos)); ...
     sum(isfinite(corr_grupo_media_parox))], ...
    {sprintf('%.3f ± %.3f', ...
        mean(corr_grupo_media_sanos, 'omitnan'), ...
        std(corr_grupo_media_sanos, 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', ...
        mean(corr_grupo_media_parox, 'omitnan'), ...
        std(corr_grupo_media_parox, 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', ...
        mean(corr_grupo_std_sanos, 'omitnan'), ...
        std(corr_grupo_std_sanos, 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', ...
        mean(corr_grupo_std_parox, 'omitnan'), ...
        std(corr_grupo_std_parox, 0, 'omitnan'))}, ...
    'VariableNames', { ...
        'Grupo', ...
        'Npacientes', ...
        'CorrGrupoMedia', ...
        'CorrGrupoStd'} ...
);

%% LIMPIEZA DE VALORES NO VÁLIDOS

correlaciones_intra_sanos = correlaciones_intra_sanos(~isnan(correlaciones_intra_sanos));
correlaciones_intra_parox = correlaciones_intra_parox(~isnan(correlaciones_intra_parox));

correlacion_media_sanos = correlacion_media_sanos(~isnan(correlacion_media_sanos));
correlacion_media_parox = correlacion_media_parox(~isnan(correlacion_media_parox));

correlacion_std_sanos = correlacion_std_sanos(~isnan(correlacion_std_sanos));
correlacion_std_parox = correlacion_std_parox(~isnan(correlacion_std_parox));

correlacion_cruzada_media = correlacion_cruzada_media(~isnan(correlacion_cruzada_media));
correlacion_cruzada_std = correlacion_cruzada_std(~isnan(correlacion_cruzada_std));

%% 3) TABLA INTRAPACIENTE

idx_sanos = strcmp(tabla_features.Grupo, 'SANO');
idx_parox = strcmp(tabla_features.Grupo, 'FA_PAROXISTICA_RS');

tabla_intrapaciente = table( ...
    {'SANO'; 'FA_PA_RS'}, ...
    [sum(idx_sanos); sum(idx_parox)], ...
    {sprintf('%.3f ± %.3f', mean(tabla_features.CorrIntraMedia(idx_sanos), 'omitnan'), std(tabla_features.CorrIntraMedia(idx_sanos), 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(tabla_features.CorrIntraMedia(idx_parox), 'omitnan'), std(tabla_features.CorrIntraMedia(idx_parox), 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(tabla_features.CorrIntraStd(idx_sanos), 'omitnan'), std(tabla_features.CorrIntraStd(idx_sanos), 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(tabla_features.CorrIntraStd(idx_parox), 'omitnan'), std(tabla_features.CorrIntraStd(idx_parox), 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(tabla_features.AmpMedia(idx_sanos), 'omitnan'), std(tabla_features.AmpMedia(idx_sanos), 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(tabla_features.AmpMedia(idx_parox), 'omitnan'), std(tabla_features.AmpMedia(idx_parox), 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(tabla_features.AmpStd(idx_sanos), 'omitnan'), std(tabla_features.AmpStd(idx_sanos), 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(tabla_features.AmpStd(idx_parox), 'omitnan'), std(tabla_features.AmpStd(idx_parox), 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(tabla_features.StdMedia(idx_sanos), 'omitnan'), std(tabla_features.StdMedia(idx_sanos), 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(tabla_features.StdMedia(idx_parox), 'omitnan'), std(tabla_features.StdMedia(idx_parox), 0, 'omitnan'))}, ...
    'VariableNames', { ...
    'Grupo', ...
    'Npacientes', ...
    'CorrIntraMedia', ...
    'CorrIntraStd', ...
    'AmpMedia', ...
    'AmpStd', ...
    'StdMedia'} ...
);

%% 4) TABLA INTRAGRUPO

tabla_intragrupo = table( ...
    {'SANO'; 'FA_PA_RS'}, ...
    [numel(correlacion_media_sanos); numel(correlacion_media_parox)], ...
    {sprintf('%.3f ± %.3f', mean(correlacion_media_sanos, 'omitnan'), std(correlacion_media_sanos, 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(correlacion_media_parox, 'omitnan'), std(correlacion_media_parox, 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(correlacion_std_sanos, 'omitnan'), std(correlacion_std_sanos, 0, 'omitnan')); ...
     sprintf('%.3f ± %.3f', mean(correlacion_std_parox, 'omitnan'), std(correlacion_std_parox, 0, 'omitnan'))}, ...
    'VariableNames', {'Grupo','Ncomparaciones','CorrMediaGrupo','CorrStdGrupo'} ...
);

%% 5) TABLA INTERGRUPOS

tabla_intergrupos = table( ...
    {'SANO vs FA_PA_RS'}, ...
    numel(correlacion_cruzada_media), ...
    {sprintf('%.3f ± %.3f', mean(correlacion_cruzada_media, 'omitnan'), std(correlacion_cruzada_media, 0, 'omitnan'))}, ...
    {sprintf('%.3f ± %.3f', mean(correlacion_cruzada_std, 'omitnan'), std(correlacion_cruzada_std, 0, 'omitnan'))}, ...
    'VariableNames', {'Grupo','Ncomparaciones','CorrMediaGrupo','CorrStdGrupo'} ...
);

%% 6) MOSTRAR TABLAS

disp(' ');
disp('TABLA INTRAPACIENTE');
disp(tabla_intrapaciente);

disp(' ');
disp('TABLA INTRAGRUPO');
disp(tabla_intragrupo);

disp(' ');
disp('TABLA INTERGRUPOS');
disp(tabla_intergrupos);


disp(' ');
disp('CORRELACIONES DE CADA PACIENTE CON EL RESTO DE SU GRUPO');
disp(tabla_corr_grupo);

disp(' ');
disp('RESUMEN DE CorrGrupoMedia Y CorrGrupoStd POR GRUPO');
disp(tabla_resumen_corr_grupo);

%% 7) BOXPLOTS

% Intrapaciente - CorrIntraMedia
f = figure('Visible','off');

datos_intra = [tabla_features.CorrIntraMedia(idx_sanos); tabla_features.CorrIntraMedia(idx_parox)];
grupos_intra = [repmat({'SANO'}, numel(tabla_features.CorrIntraMedia(idx_sanos)), 1); ...
                repmat({'FA_PA_RS'}, numel(tabla_features.CorrIntraMedia(idx_parox)), 1)];

boxplot(datos_intra, grupos_intra);
title('Media de correlación de ondas T por paciente', 'Interpreter','none');
ylabel('Correlación', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_intrapaciente_media.png'));
close(f)

% Intrapaciente - CorrIntraStd
f = figure('Visible','off');

datos_intra_std = [tabla_features.CorrIntraStd(idx_sanos); tabla_features.CorrIntraStd(idx_parox)];
grupos_intra_std = [repmat({'SANO'}, numel(tabla_features.CorrIntraStd(idx_sanos)), 1); ...
                    repmat({'FA_PA_RS'}, numel(tabla_features.CorrIntraStd(idx_parox)), 1)];

boxplot(datos_intra_std, grupos_intra_std);
title('Desviación de correlación de ondas T por paciente', 'Interpreter','none');
ylabel('Correlación', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_intrapaciente_std.png'));
close(f)

% AmpMedia
f = figure('Visible','off');

datos_amp_media = [tabla_features.AmpMedia(idx_sanos); tabla_features.AmpMedia(idx_parox)];
grupos_amp_media = [repmat({'SANO'}, sum(idx_sanos), 1); ...
                    repmat({'FA_PA_RS'}, sum(idx_parox), 1)];

boxplot(datos_amp_media, grupos_amp_media);
title('Amplitud media de la onda T', 'Interpreter','none');
ylabel('AmpMedia (mV)', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_AmpMedia.png'));
close(f)

% AmpStd
f = figure('Visible','off');

datos_amp_std = [tabla_features.AmpStd(idx_sanos); tabla_features.AmpStd(idx_parox)];
grupos_amp_std = [repmat({'SANO'}, sum(idx_sanos), 1); ...
                  repmat({'FA_PA_RS'}, sum(idx_parox), 1)];

boxplot(datos_amp_std, grupos_amp_std);
title('Variabilidad de amplitud de la onda T', 'Interpreter','none');
ylabel('AmpStd (mV)', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_AmpStd.png'));
close(f)

% StdMedia
f = figure('Visible','off');

datos_std_media = [tabla_features.StdMedia(idx_sanos); tabla_features.StdMedia(idx_parox)];
grupos_std_media = [repmat({'SANO'}, sum(idx_sanos), 1); ...
                    repmat({'FA_PA_RS'}, sum(idx_parox), 1)];

boxplot(datos_std_media, grupos_std_media);
title('Variabilidad media de la morfología de la onda T', 'Interpreter','none');
ylabel('StdMedia (mV)', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_StdMedia.png'));
close(f)

% CorrGrupoMedia: T_media, un valor por paciente
f = figure('Visible','off');

datos_corr_grupo_media = [corr_grupo_media_sanos; corr_grupo_media_parox];
grupos_corr_grupo_media = [ ...
    repmat({'SANO'}, numel(corr_grupo_media_sanos), 1); ...
    repmat({'FA_PA_RS'}, numel(corr_grupo_media_parox), 1)];

boxplot(datos_corr_grupo_media, grupos_corr_grupo_media);
title('Similitud de la onda T media con el resto del grupo', 'Interpreter','none');
ylabel('Correlación media', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_CorrGrupoMedia.png'));
close(f)

% CorrGrupoStd: T_std, un valor por paciente
f = figure('Visible','off');

datos_corr_grupo_std = [corr_grupo_std_sanos; corr_grupo_std_parox];
grupos_corr_grupo_std = [ ...
    repmat({'SANO'}, numel(corr_grupo_std_sanos), 1); ...
    repmat({'FA_PA_RS'}, numel(corr_grupo_std_parox), 1)];

boxplot(datos_corr_grupo_std, grupos_corr_grupo_std);
title('Similitud de T_std con el resto del grupo', 'Interpreter','none');
ylabel('Correlación media', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_CorrGrupoStd.png'));
close(f)

% Correlación T_media: intragrupo e intergrupos
f = figure('Visible','off');

datos_media = [ ...
    correlacion_media_sanos; ...
    correlacion_media_parox; ...
    correlacion_cruzada_media];

grupos_media = [ ...
    repmat({'SANO'}, numel(correlacion_media_sanos), 1); ...
    repmat({'FA_PA_RS'}, numel(correlacion_media_parox), 1); ...
    repmat({'SANO vs FA_PA_RS'}, numel(correlacion_cruzada_media), 1)];

boxplot(datos_media, grupos_media);
title('Correlación de T_media entre pacientes', 'Interpreter','none');
ylabel('Correlación', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_correlacion_T_media_intra_inter.png'));
close(f)

% Correlación T_std: intragrupo e intergrupos
f = figure('Visible','off');

datos_std = [ ...
    correlacion_std_sanos; ...
    correlacion_std_parox; ...
    correlacion_cruzada_std];

grupos_std = [ ...
    repmat({'SANO'}, numel(correlacion_std_sanos), 1); ...
    repmat({'FA_PA_RS'}, numel(correlacion_std_parox), 1); ...
    repmat({'SANO vs FA_PA_RS'}, numel(correlacion_cruzada_std), 1)];

boxplot(datos_std, grupos_std);
title('Correlación de T_std entre pacientes', 'Interpreter','none');
ylabel('Correlación', 'Interpreter','none', 'FontSize', 14);
xlabel('Grupo', 'Interpreter','none');
grid on

saveas(f, fullfile(carpeta_boxplots, 'boxplot_correlacion_T_std_intra_inter.png'));
close(f)


%% 8) GUARDAR TABLAS

archivo_excel = fullfile(carpeta_out, 'tablas_correlaciones_ondaT.xlsx');

% En el Excel de correlaciones se guardan:
%   - Intrapaciente
%   - Intragrupo
%   - Intergrupos
%   - Resumen_CorrGrupoMedia, con CorrGrupoMedia y CorrGrupoStd
writetable(tabla_intrapaciente, archivo_excel, 'Sheet', 'Intrapaciente');
writetable(tabla_intragrupo, archivo_excel, 'Sheet', 'Intragrupo');
writetable(tabla_intergrupos, archivo_excel, 'Sheet', 'Intergrupos');
writetable(tabla_resumen_corr_grupo, archivo_excel, 'Sheet', 'Resumen_CorrGrupoMedia');

% Guardar en features_ondaT.xlsx una fila por paciente,
% incluyendo CorrGrupoMedia y CorrGrupoStd
archivo_features = fullfile(carpeta_out, 'features_ondaT.xlsx');
writetable(tabla_features, archivo_features);

fprintf('Excel de correlaciones guardado en:\n%s\n', archivo_excel);
fprintf('Features por paciente guardadas en:\n%s\n', archivo_features);

%% 9) GUARDAR RESULTADOS DE CORRELACIONES

save(fullfile(carpeta_out, 'correlaciones_ondaT.mat'), ...
    'tabla_features', ...
    'tabla_intrapaciente', 'tabla_corr_grupo', ...
    'tabla_resumen_corr_grupo', ...
    'tabla_intragrupo', 'tabla_intergrupos', ...
    'correlaciones_intra_sanos', 'correlaciones_intra_parox', ...
    'correlacion_media_sanos', 'correlacion_media_parox', ...
    'correlacion_std_sanos', 'correlacion_std_parox', ...
    'correlacion_cruzada_media', 'correlacion_cruzada_std', ...
    'corr_grupo_media_sanos', 'corr_grupo_media_parox', ...
    'corr_grupo_std_sanos', 'corr_grupo_std_parox', ...
    'M_sanos_media', 'M_parox_media', 'M_cruzada_media', ...
    'M_sanos_std', 'M_parox_std', 'M_cruzada_std');

fprintf('Resultados de correlaciones guardados en .mat\n');

fprintf('\nFIN SCRIPT DE TABLAS, HEATMAPS Y BOXPLOTS DE ONDA T\n');