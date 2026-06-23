clear
clc
close all

%% ============================================================
% UMAP:
% FA PAROXÍSTICA DURANTE EL EPISODIO vs FA PERSISTENTE
%
% Biomarcadores seleccionados:
%   1. RR_mean
%   2. SDNN
%   3. RMSSD
%   4. pNN20
%   5. CV_RR
%   6. SD1_SD2_ratio
%   7. SampEn
%   8. LF_HF
%   9. DF_completo_Hz
%  10. DF_residual_Hz
%
% Criterio:
%   - Debido al reducido número de biomarcadores significativos tras FDR,
%     se partió de todos los biomarcadores incluidos en el análisis estadístico.
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

idx_FA_paroxistica = ...
    T.Grupo == "FA_PAROXISTICA_FA" | ...
    T.Grupo == "FA_PAROXISTICA_AF";

idx_FA_persistente = ...
    T.Grupo == "FA_PERSISTENTE";

idx_grupos = idx_FA_paroxistica | idx_FA_persistente;

T_umap = T(idx_grupos, :);

%% BIOMARCADORES SELECCIONADOS

biomarcadores = { ...
    'RR_mean', ...
    'SDNN', ...
    'RMSSD', ...
    'pNN20', ...
    'CV_RR', ...
    'SD1_SD2_ratio', ...
    'SampEn', ...
    'LF_HF', ...
    'DF_completo_Hz', ...
    'DF_residual_Hz'};

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

grupo_figura( ...
    grupo == "FA_PAROXISTICA_FA" | ...
    grupo == "FA_PAROXISTICA_AF") = ...
    "FA paroxística en episodio";

grupo_figura(grupo == "FA_PERSISTENTE") = ...
    "FA persistente";

%% MOSTRAR EL NÚMERO DE PACIENTES

fprintf('\nPacientes incluidos en el UMAP:\n');

fprintf('FA paroxística en episodio: %d\n', ...
    sum(grupo_figura == "FA paroxística en episodio"));

fprintf('FA persistente: %d\n', ...
    sum(grupo_figura == "FA persistente"));

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

idx1 = grupo_figura == "FA paroxística en episodio";
idx2 = grupo_figura == "FA persistente";

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
    {'UMAP: FA paroxística en FA frente a FA persistente'}, ...
    'FontSize', 17, ...
    'FontWeight', 'normal');

legend( ...
    {'FA paroxística en episodio', 'FA persistente'}, ...
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
    'UMAP_2D_FA_PAROXISTICA_FA_vs_FA_PERSISTENTE.png');

exportgraphics(f, archivo_salida, ...
    'Resolution', 300);

fprintf('Figura guardada en:\n%s\n', archivo_salida);