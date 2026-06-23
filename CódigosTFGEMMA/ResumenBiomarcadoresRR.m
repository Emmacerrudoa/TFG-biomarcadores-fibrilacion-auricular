clear
clc
close all

%% ============================================================
% VISUALIZACIÓN Y RESUMEN DE BIOMARCADORES RR
%
% Este script carga los biomarcadores RR calculados previamente y genera
% tablas resumen, boxplots por grupo, histogramas medios de intervalos RR
% y figuras de Poincaré por registro.
%
% Entradas:
%   - biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia.xlsx
%   - resultados_RR_4grupos.mat
%
% Salidas:
%   - biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia_y_corr.xlsx
%   - resumen_biomarcadores_RR_por_grupo.xlsx
%   - boxplots de biomarcadores RR
%   - histograma medio RR por grupo
%   - figuras de Poincaré por registro
%
% Requisitos:
%   - Resultados generados previamente por el script de análisis RR.
%   - MATLAB con funciones de representación gráfica.
%% ============================================================

%% CONFIGURACIÓN GENERAL

% Modificar esta ruta según la ubicación local de los resultados RR.
carpeta_out = 'C:\Users\Emma\Documents\MATLAB\analisisRR_hibrido';

archivo_excel = fullfile(carpeta_out, 'biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia.xlsx');
archivo_mat   = fullfile(carpeta_out, 'resultados_RR_4grupos.mat');

carpeta_boxplots = fullfile(carpeta_out, 'boxplots_biomarcadores');
if ~exist(carpeta_boxplots, 'dir')
    mkdir(carpeta_boxplots);
end

carpeta_poincare = fullfile(carpeta_out, 'figuras_poincare_desde_mat');
if ~exist(carpeta_poincare, 'dir')
    mkdir(carpeta_poincare);
end

%% 1. CARGAR ARCHIVOS DE ENTRADA

if ~exist(archivo_excel, 'file')
    error('No se encuentra el Excel: %s', archivo_excel);
end

if ~exist(archivo_mat, 'file')
    error('No se encuentra el MAT: %s', archivo_mat);
end

T_bio = readtable(archivo_excel);
S = load(archivo_mat);

if ~isfield(S, 'resultados_por_grupo') || ~isfield(S, 'edges_rr')
    error('El archivo .mat no contiene todas las variables necesarias.');
end

resultados_por_grupo = S.resultados_por_grupo;
edges_rr = S.edges_rr;

if isempty(T_bio)
    error('El Excel está vacío.');
end

%% 2. DEFINIR ORDEN DE LOS GRUPOS SEGÚN EL ARCHIVO .MAT

orden_grupos_mat = cell(1, numel(resultados_por_grupo));

for g = 1:numel(resultados_por_grupo)
    orden_grupos_mat{g} = resultados_por_grupo(g).grupo;
end

orden_grupos_abrev = cell(size(orden_grupos_mat));

for g = 1:numel(orden_grupos_mat)
    orden_grupos_abrev{g} = nombre_grupo_abreviado(orden_grupos_mat{g});
end

%% 3. RECALCULAR CORRELACIONES DE HISTOGRAMAS POR GRUPO

correlaciones_grupos = struct([]);

for g = 1:numel(resultados_por_grupo)

    nombre_grupo = resultados_por_grupo(g).grupo;
    Rg = resultados_por_grupo(g).resultados;

    if isempty(Rg)
        correlaciones_grupos(g).grupo = nombre_grupo;
        correlaciones_grupos(g).H = [];
        correlaciones_grupos(g).Rhist = [];
        correlaciones_grupos(g).nombres = {};
        correlaciones_grupos(g).media_hist = [];
        correlaciones_grupos(g).std_hist = [];
        continue
    end

    n = numel(Rg);
    L = length(Rg(1).hist_rr);

    H = nan(n, L);
    nombres = cell(n,1);

    for i = 1:n
        if length(Rg(i).hist_rr) == L
            H(i,:) = Rg(i).hist_rr;
        else
            warning('Histograma de %s tiene longitud inesperada (%d vs %d), fila omitida.', ...
                Rg(i).registro, length(Rg(i).hist_rr), L);
        end
        nombres{i} = Rg(i).registro;
    end

    media_hist = mean(H, 1, 'omitnan');
    std_hist = std(H, 0, 1, 'omitnan');

    if n > 1
        Rhist = corrcoef(H');
    else
        Rhist = NaN;
    end

    correlaciones_grupos(g).grupo = nombre_grupo;
    correlaciones_grupos(g).H = H;
    correlaciones_grupos(g).Rhist = Rhist;
    correlaciones_grupos(g).nombres = nombres;
    correlaciones_grupos(g).media_hist = media_hist;
    correlaciones_grupos(g).std_hist = std_hist;
end

%% 4. RECALCULAR CORRELACIÓN MEDIA DEL HISTOGRAMA RR POR PACIENTE

CorrelacionHistRRGrupo = nan(height(T_bio),1);

for g = 1:numel(resultados_por_grupo)

    nombre_grupo = resultados_por_grupo(g).grupo;
    Rg = resultados_por_grupo(g).resultados;
    Rhist_g = correlaciones_grupos(g).Rhist;

    if isempty(Rg)
        continue
    end

    if numel(Rg) == 1 || (numel(Rhist_g) == 1 && isnan(Rhist_g))
        idx_tabla = strcmp(T_bio.Grupo, nombre_grupo) & strcmp(T_bio.Paciente, Rg(1).registro);
        CorrelacionHistRRGrupo(idx_tabla) = NaN;
        continue
    end

    for i = 1:numel(Rg)
        fila_corr = Rhist_g(i,:);
        fila_corr(i) = NaN;
        valor_medio = mean(fila_corr, 'omitnan');

        idx_tabla = strcmp(T_bio.Grupo, nombre_grupo) & strcmp(T_bio.Paciente, Rg(i).registro);
        CorrelacionHistRRGrupo(idx_tabla) = valor_medio;
    end
end

T_bio.CorrelacionHistRRGrupo = CorrelacionHistRRGrupo;

%% 5. TABLA RESUMEN FINAL POR GRUPO

Grupo = {};
N_pacientes = [];
N_intervalos_RR = [];

RR_mean_txt = {};
SDNN_txt = {};
RMSSD_txt = {};
SDSD_txt = {};
pNN50_txt = {};
pNN20_txt = {};
CV_RR_txt = {};
SD1_txt = {};
SD2_txt = {};
SD1_SD2_ratio_txt = {};
CorrHist_txt = {};
SampEn_txt = {};
LF_txt = {};
HF_txt = {};
LF_HF_txt = {};
LFnu_txt = {};
HFnu_txt = {};

for g = 1:numel(orden_grupos_mat)

    grupo_actual = orden_grupos_mat{g};
    idx = strcmp(T_bio.Grupo, grupo_actual);
    Tg = T_bio(idx, :);

    if isempty(Tg)
        continue
    end

    Grupo{end+1,1} = nombre_grupo_abreviado(grupo_actual);
    N_pacientes(end+1,1) = height(Tg);
    N_intervalos_RR(end+1,1) = sum(Tg.NumeroRR, 'omitnan');

    RR_mean_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.RR_mean, 'omitnan'), std(Tg.RR_mean, 0, 'omitnan'));

    SDNN_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SDNN, 'omitnan'), std(Tg.SDNN, 0, 'omitnan'));

    RMSSD_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.RMSSD, 'omitnan'), std(Tg.RMSSD, 0, 'omitnan'));

    SDSD_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SDSD, 'omitnan'), std(Tg.SDSD, 0, 'omitnan'));

    pNN50_txt{end+1,1} = sprintf('%.2f ± %.2f', ...
        mean(Tg.pNN50, 'omitnan'), std(Tg.pNN50, 0, 'omitnan'));

    pNN20_txt{end+1,1} = sprintf('%.2f ± %.2f', ...
        mean(Tg.pNN20, 'omitnan'), std(Tg.pNN20, 0, 'omitnan'));

    CV_RR_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.CV_RR, 'omitnan'), std(Tg.CV_RR, 0, 'omitnan'));

    SD1_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SD1, 'omitnan'), std(Tg.SD1, 0, 'omitnan'));

    SD2_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SD2, 'omitnan'), std(Tg.SD2, 0, 'omitnan'));

    SD1_SD2_ratio_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SD1_SD2_ratio, 'omitnan'), std(Tg.SD1_SD2_ratio, 0, 'omitnan'));

    CorrHist_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.CorrelacionHistRRGrupo, 'omitnan'), std(Tg.CorrelacionHistRRGrupo, 0, 'omitnan'));

    SampEn_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.SampEn, 'omitnan'), std(Tg.SampEn, 0, 'omitnan'));

    LF_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.LF, 'omitnan'), std(Tg.LF, 0, 'omitnan'));

    HF_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.HF, 'omitnan'), std(Tg.HF, 0, 'omitnan'));

    LF_HF_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.LF_HF, 'omitnan'), std(Tg.LF_HF, 0, 'omitnan'));

    LFnu_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.LFnu, 'omitnan'), std(Tg.LFnu, 0, 'omitnan'));

    HFnu_txt{end+1,1} = sprintf('%.3f ± %.3f', ...
        mean(Tg.HFnu, 'omitnan'), std(Tg.HFnu, 0, 'omitnan'));
end

T_resumen = table( ...
    Grupo, N_pacientes, N_intervalos_RR, ...
    RR_mean_txt, SDNN_txt, RMSSD_txt, SDSD_txt, pNN50_txt, pNN20_txt, ...
    CV_RR_txt, SD1_txt, SD2_txt, SD1_SD2_ratio_txt, ...
    SampEn_txt, LF_txt, HF_txt, LF_HF_txt, LFnu_txt, HFnu_txt, CorrHist_txt, ...
    'VariableNames', {'Grupo','N_pacientes','N_intervalos_RR','RR_mean','SDNN','RMSSD','SDSD', ...
    'pNN50','pNN20','CV_RR','SD1','SD2','SD1_SD2_ratio', ...
    'SampEn','LF','HF','LF_HF','LFnu','HFnu','CorrelacionHistRRGrupo'} ...
);

writetable(T_bio, fullfile(carpeta_out, ...
    'biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia_y_corr.xlsx'));

disp(T_resumen)

writetable(T_resumen, fullfile(carpeta_out, ...
    'resumen_biomarcadores_RR_por_grupo.xlsx'));

%% 6. BOXPLOTS DE BIOMARCADORES POR GRUPO CON UNIDADES

variables_bio = { ...
    'RR_mean', 'SDNN', 'RMSSD', 'SDSD', 'pNN50', 'pNN20', ...
    'CV_RR', 'SD1', 'SD2', 'SD1_SD2_ratio', ...
    'SampEn', 'LF', 'HF', 'LF_HF', 'LFnu', ...
    'HFnu', 'CorrelacionHistRRGrupo'};

for v = 1:numel(variables_bio)

    nombre_var = variables_bio{v};

    datos_var = T_bio.(nombre_var);
    idx_validos = ~isnan(datos_var);

    grupos_orig = T_bio.Grupo(idx_validos);
    datos_plot = datos_var(idx_validos);

    grupos_abrev = cell(size(grupos_orig));

    for ii = 1:numel(grupos_orig)
        grupos_abrev{ii} = nombre_grupo_abreviado(grupos_orig{ii});
    end

    grupos_cat = categorical(grupos_abrev, orden_grupos_abrev, 'Ordinal', true);

    [nombre_label, unidad] = etiqueta_variable_unidad_RR(nombre_var);

    if isempty(unidad)
        etiqueta_y = nombre_label;
    else
        etiqueta_y = [nombre_label ' (' unidad ')'];
    end

    f = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

    boxplot(datos_plot, grupos_cat, ...
        'Symbol', 'r+', ...
        'Widths', 0.55);

    grid on
    box on

    ax = gca;

    set(ax, ...
        'FontSize', 17, ...
        'LineWidth', 1.1);

    t = title(ax, ['Distribución de ' nombre_label ' por grupo'], ...
        'Interpreter','none');
    t.FontSize = 20;
    t.FontWeight = 'bold';

    yl = ylabel(ax, etiqueta_y, ...
        'Interpreter','none');
    yl.FontSize = 24;
    yl.FontWeight = 'normal';

    xl = xlabel(ax, 'Grupo', ...
        'Interpreter','none');
    xl.FontSize = 17;
    xl.FontWeight = 'normal';

    saveas(f, fullfile(carpeta_boxplots, ['boxplot_' nombre_var '.png']));
    close(f)
end

%% 7. HISTOGRAMA MEDIO RR POR GRUPO

f = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);
hold on

centros_bins = edges_rr(1:end-1) + diff(edges_rr)/2;

nombres_leyenda = {};

for g = 1:numel(correlaciones_grupos)

    if isempty(correlaciones_grupos(g).media_hist)
        continue
    end

    plot(centros_bins, correlaciones_grupos(g).media_hist, ...
        'LineWidth', 2);

    nombres_leyenda{end+1} = nombre_grupo_abreviado(correlaciones_grupos(g).grupo); %#ok<SAGROW>
end

if ~isempty(nombres_leyenda)
    legend(nombres_leyenda, ...
        'Location', 'northeast', ...
        'Interpreter', 'none', ...
        'FontSize', 11);
end

grid on
box on

ax = gca;

set(ax, ...
    'FontSize', 16, ...
    'LineWidth', 1.1);

xl = xlabel(ax, 'Intervalo RR (s)', ...
    'Interpreter','none');
xl.FontSize = 22;
xl.FontWeight = 'normal';

yl = ylabel(ax, 'Probabilidad', ...
    'Interpreter','none');
yl.FontSize = 24;
yl.FontWeight = 'normal';

t = title(ax, 'Histograma medio de intervalos RR por grupo', ...
    'Interpreter','none');
t.FontSize = 20;
t.FontWeight = 'bold';

saveas(f, fullfile(carpeta_out, 'histogramas_medios_RR_por_grupo.png'));
close(f)

%% 8. FIGURAS DE POINCARÉ DESDE EL ARCHIVO .MAT

for g = 1:numel(resultados_por_grupo)

    nombre_grupo = resultados_por_grupo(g).grupo;
    Rg = resultados_por_grupo(g).resultados;

    for i = 1:numel(Rg)

        if ~isfield(Rg(i), 'RR') || isempty(Rg(i).RR)
            continue
        end

        RR = Rg(i).RR(:);

        if numel(RR) < 50
            continue
        end

        RRn  = [];
        RRn1 = [];

        if isfield(Rg(i), 'segmentos_RR') && ~isempty(Rg(i).segmentos_RR)

            for si = 1:numel(Rg(i).segmentos_RR)

                seg = Rg(i).segmentos_RR{si}(:);

                if numel(seg) >= 2
                    RRn  = [RRn;  seg(1:end-1)]; %#ok<AGROW>
                    RRn1 = [RRn1; seg(2:end)];   %#ok<AGROW>
                end
            end

        else
            warning('Poincare omitido para %s (%s): segmentos_RR no disponibles.', ...
                Rg(i).registro, nombre_grupo);
            continue
        end

        if isempty(RRn) || isempty(RRn1)
            warning('Poincare omitido para %s (%s): no hay pares RR válidos.', ...
                Rg(i).registro, nombre_grupo);
            continue
        end

        xp = (RRn + RRn1) / sqrt(2);
        yp = (RRn1 - RRn) / sqrt(2);

        cx = mean(xp, 'omitnan');
        cy = mean(yp, 'omitnan');

        SD1 = Rg(i).SD1;
        SD2 = Rg(i).SD2;

        nombre_grupo_fig = nombre_grupo_largo(nombre_grupo);

        f = figure('Visible','off', ...
                   'Color','w', ...
                   'Position', [100 100 1500 700]);

        %% POINCARÉ CLÁSICO

        ax1 = subplot(1,2,1);

        plot(ax1, RRn, RRn1, 'k.', 'MarkerSize', 10)
        hold(ax1, 'on')

        plot(ax1, mean(RR, 'omitnan'), mean(RR, 'omitnan'), ...
            'ro', ...
            'MarkerFaceColor', 'r', ...
            'MarkerSize', 9)

        xlabel(ax1, 'RR_n (s)', ...
            'FontSize', 34, ...
            'Interpreter', 'tex', ...
            'FontWeight', 'normal')

        ylabel(ax1, 'RR_{n+1} (s)', ...
            'FontSize', 34, ...
            'Interpreter', 'tex', ...
            'FontWeight', 'normal')

        title(ax1, 'Poincaré clásico', ...
            'FontSize', 28, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none')

        axis(ax1, 'equal')
        grid(ax1, 'on')
        box(ax1, 'on')

        set(ax1, ...
            'FontSize', 24, ...      % números de los ejes
            'LineWidth', 1.2)

        %% POINCARÉ ROTADO

        ax2 = subplot(1,2,2);

        plot(ax2, xp, yp, 'b.', 'MarkerSize', 10)
        hold(ax2, 'on')

        theta = linspace(0, 2*pi, 200);

        plot(ax2, cx + SD2*cos(theta), cy + SD1*sin(theta), ...
            'r', ...
            'LineWidth', 2.2)

        xlabel(ax2, 'Eje largo, SD2 (s)', ...
            'FontSize', 34, ...
            'Interpreter', 'none', ...
            'FontWeight', 'normal')

        ylabel(ax2, 'Eje corto, SD1 (s)', ...
            'FontSize', 34, ...
            'Interpreter', 'none', ...
            'FontWeight', 'normal')

        title(ax2, sprintf('Poincaré rotado: SD1 = %.3f s, SD2 = %.3f s', SD1, SD2), ...
            'FontSize', 28, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none')

        axis(ax2, 'equal')
        grid(ax2, 'on')
        box(ax2, 'on')

        set(ax2, ...
            'FontSize', 24, ...      % números de los ejes
            'LineWidth', 1.2)

        %% TÍTULO GENERAL

        sgtitle(sprintf('%s - %s', nombre_grupo_fig, Rg(i).registro), ...
            'Interpreter', 'none', ...
            'FontSize', 30, ...
            'FontWeight', 'bold')

        %% GUARDAR FIGURA

        nombre_fig = ['Poincare_' nombre_grupo_abreviado(nombre_grupo) '_' Rg(i).registro '.png'];
        ruta_fig = fullfile(carpeta_poincare, nombre_fig);

        exportgraphics(f, ruta_fig, 'Resolution', 300);

        close(f)
    end
end

%% ============================================================
% FUNCIÓN LOCAL: ETIQUETAS Y UNIDADES
%% ============================================================

function [nombre_label, unidad] = etiqueta_variable_unidad_RR(nombre_var)

    switch nombre_var

        case 'RR_mean'
            nombre_label = 'RR mean';
            unidad = 's';

        case 'SDNN'
            nombre_label = 'SDNN';
            unidad = 's';

        case 'RMSSD'
            nombre_label = 'RMSSD';
            unidad = 's';

        case 'SDSD'
            nombre_label = 'SDSD';
            unidad = 's';

        case 'pNN50'
            nombre_label = 'pNN50';
            unidad = '%';

        case 'pNN20'
            nombre_label = 'pNN20';
            unidad = '%';

        case 'CV_RR'
            nombre_label = 'CV RR';
            unidad = '';

        case 'SD1'
            nombre_label = 'SD1';
            unidad = 's';

        case 'SD2'
            nombre_label = 'SD2';
            unidad = 's';

        case 'SD1_SD2_ratio'
            nombre_label = 'SD1/SD2';
            unidad = '';

        case 'SampEn'
            nombre_label = 'SampEn';
            unidad = '';

        case 'LF'
            nombre_label = 'LF';
            unidad = 's^2';

        case 'HF'
            nombre_label = 'HF';
            unidad = 's^2';

        case 'LF_HF'
            nombre_label = 'LF/HF';
            unidad = '';

        case 'LFnu'
            nombre_label = 'LFnu';
            unidad = '%';

        case 'HFnu'
            nombre_label = 'HFnu';
            unidad = '%';

        case 'CorrelacionHistRRGrupo'
            nombre_label = 'Correlación histograma RR';
            unidad = '';

        otherwise
            nombre_label = nombre_var;
            unidad = '';
    end
end

%% ============================================================
% FUNCIÓN LOCAL: ABREVIAR NOMBRE DEL GRUPO
%% ============================================================

function nombre = nombre_grupo_abreviado(grupo)

    switch grupo
        case 'SANO'
            nombre = 'SANO';

        case {'FA_PAROXISTICA_SR', 'FA_PAROXISTICA_RS'}
            nombre = 'FA_PA_RS';

        case {'FA_PAROXISTICA_AF', 'FA_PAROXISTICA_FA'}
            nombre = 'FA_PA_FA';

        case 'FA_PERSISTENTE'
            nombre = 'FA_PE';

        otherwise
            nombre = grupo;
    end
end

%% ============================================================
% FUNCIÓN LOCAL: NOMBRE LARGO DE GRUPO
%% ============================================================

function nombre = nombre_grupo_largo(grupo)

    switch grupo

        case 'SANO'
            nombre = 'Sano';

        case {'FA_PAROXISTICA_SR', 'FA_PAROXISTICA_RS'}
            nombre = 'FA paroxística (RS)';

        case {'FA_PAROXISTICA_AF', 'FA_PAROXISTICA_FA'}
            nombre = 'FA paroxística (FA)';

        case 'FA_PERSISTENTE'
            nombre = 'FA persistente';

        otherwise
            nombre = grupo;
    end
end