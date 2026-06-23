clear
clc
close all

%% ============================================================
% UMAP:
% SANO vs FA PAROXÍSTICA EN RITMO SINUSAL
%
% Biomarcadores seleccionados:
%   1. SD1_SD2_ratio
%   2. LFnu
%   3. P_CorrIntraMedia
%   4. RMSSD
%   5. DF_residual_Hz
%   6. T_AmpMedia
%   7. T_CorrIntraMedia
%   8. LF
%   9. SD2
%  10. CV_RR
%
% Criterio:
%   - Se partió de los biomarcadores significativos tras la corrección FDR.
%   - Se descartaron biomarcadores redundantes con |rho| >= 0.90.
%% ============================================================

rng(42)

%% RUTA DE LA TABLA MAESTRA

archivo_tabla = ...
    'C:\Users\Emma\Documents\MATLAB\RESULTADOS_ESTADISTICOS_TFG\tabla_maestra_biomarcadores_TFGnfft.xlsx';

carpeta_salida = fileparts(archivo_tabla);

%% LEER TABLA

T = readtable(archivo_tabla);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

T.Grupo = strtrim(string(T.Grupo));

%% SELECCIONAR LOS DOS GRUPOS

idx_sano = T.Grupo == "SANO";

idx_FA_RS = ...
    T.Grupo == "FA_PAROXISTICA_RS" | ...
    T.Grupo == "FA_PAROXISTICA_SR";

idx_grupos = idx_sano | idx_FA_RS;

T_umap = T(idx_grupos, :);

%% BIOMARCADORES SELECCIONADOS

biomarcadores = { ...
    'SD1_SD2_ratio', ...
    'LFnu', ...
    'P_CorrIntraMedia', ...
    'RMSSD', ...
    'DF_residual_Hz', ...
    'T_AmpMedia', ...
    'T_CorrIntraMedia', ...
    'LF', ...
    'SD2', ...
    'CV_RR'};

%% COMPROBAR QUE LAS COLUMNAS EXISTEN

for i = 1:numel(biomarcadores)

    if ~ismember(biomarcadores{i}, ...
            T_umap.Properties.VariableNames)

        error('No se encuentra la columna: %s', ...
            biomarcadores{i});
    end
end

%% CREAR MATRIZ DE DATOS

X = T_umap{:, biomarcadores};

grupo = T_umap.Grupo;

%% ELIMINAR FILAS CON NaN O Inf

filas_validas = all(isfinite(X), 2);

X = X(filas_validas, :);
grupo = grupo(filas_validas);

%% NOMBRES MOSTRADOS EN LA FIGURA

grupo_figura = strings(size(grupo));

grupo_figura(grupo == "SANO") = ...
    "Sano";

grupo_figura( ...
    grupo == "FA_PAROXISTICA_RS" | ...
    grupo == "FA_PAROXISTICA_SR") = ...
    "FA paroxística en RS";

%% MOSTRAR EL NÚMERO DE PACIENTES

fprintf('\nPacientes incluidos en el UMAP:\n');

fprintf('Sanos: %d\n', ...
    sum(grupo_figura == "Sano"));

fprintf('FA paroxística en RS: %d\n', ...
    sum(grupo_figura == "FA paroxística en RS"));

fprintf('Total: %d\n\n', size(X, 1));

%% COMPROBAR VARIABILIDAD DE LAS COLUMNAS

desviaciones = std(X, 0, 1);

if any(desviaciones == 0)

    columnas_constantes = biomarcadores(desviaciones == 0);

    error( ...
        'Hay biomarcadores sin variabilidad: %s', ...
        strjoin(columnas_constantes, ', '));
end

%% NORMALIZAR

Xz = zscore(X);

%% APLICAR UMAP

[reduction, umap_obj] = run_umap(Xz, ...
    'n_components', 2, ...
    'n_neighbors', 10, ...
    'min_dist', 0.4, ...
    'metric', 'euclidean', ...
    'randomize', false);

%% ÍNDICES PARA LA FIGURA

idx1 = grupo_figura == "Sano";
idx2 = grupo_figura == "FA paroxística en RS";

%% REPRESENTACIÓN

f = figure( ...
    'Color', 'w', ...
    'Position', [100 100 900 700]);

hold on

scatter( ...
    reduction(idx1, 1), ...
    reduction(idx1, 2), ...
    85, ...
    [0.00 0.45 0.74], ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.7);

scatter( ...
    reduction(idx2, 1), ...
    reduction(idx2, 2), ...
    85, ...
    [0.85 0.33 0.10], ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.7);

xlabel('UMAP 1', ...
    'FontSize', 17);

ylabel('UMAP 2', ...
    'FontSize', 17);

title( ...
    {'UMAP: sujetos sanos frente a FA paroxística en RS'}, ...
    'FontSize', 17, ...
    'FontWeight', 'normal');

legend( ...
    {'Sano', 'FA paroxística en RS'}, ...
    'Location', 'best', ...
    'FontSize', 13);

ax = gca;

set(ax, ...
    'FontSize', 15, ...
    'LineWidth', 1.1, ...
    'TickDir', 'out');

grid on
box on
axis equal

hold off

%% GUARDAR FIGURA

archivo_salida = fullfile( ...
    carpeta_salida, ...
    'UMAP_2D_SANO_vs_FA_PAROXISTICA_RS.png');

exportgraphics(f, archivo_salida, ...
    'Resolution', 300);

fprintf('Figura guardada en:\n%s\n', archivo_salida);