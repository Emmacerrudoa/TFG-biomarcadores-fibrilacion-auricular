clear; 
clc; 
close all;

%% ============================================================
% HEATMAP DE CORRELACIÓN ENTRE BIOMARCADORES EN FA
% Comparación: FA PAROXÍSTICA EN EPISODIO vs FA PERSISTENTE
%
% Debido al reducido número de biomarcadores significativos tras FDR,
% se incluyen todos los biomarcadores utilizados en el análisis estadístico
% de esta comparación.
%
% El objetivo es estudiar la correlación entre biomarcadores y detectar
% posibles variables redundantes antes de aplicar UMAP.
%
% Solo se guarda el heatmap en formato PNG.
%% ============================================================

%% RUTAS

% Modificar esta ruta según la ubicación local de la tabla maestra.
archivo_tabla = ...
    'C:\Users\Emma\Documents\MATLAB\RESULTADOS_ESTADISTICOS_TFG\tabla_maestra_biomarcadores_TFG.xlsx';

carpeta_salida = fileparts(archivo_tabla);

%% LEER TABLA

T = readtable(archivo_tabla);

T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

T.Grupo = strtrim(string(T.Grupo));

%% DEFINIR COMPARACIÓN

nombre_comp = "FA_PAROXISTICA_FA_vs_FA_PERSISTENTE";
grupos = ["FA_PAROXISTICA_FA", "FA_PERSISTENTE"];

fprintf('\n============================================\n');
fprintf('Comparación: %s\n', nombre_comp);
fprintf('============================================\n');

%% 1) Filtrar filas de los grupos correspondientes

idx_grupo = ismember(T.Grupo, grupos);
Tsub = T(idx_grupo, :);

fprintf('\nFilas totales usadas: %d\n', height(Tsub));

for g = 1:numel(grupos)

    fprintf('  %-25s %d sujetos\n', ...
        char(grupos(g)), ...
        sum(Tsub.Grupo == grupos(g)));
end

%% 2) Seleccionar biomarcadores disponibles durante FA

biomarcadores = [
    "RR_mean"
    "SDNN"
    "RMSSD"
    "SDSD"
    "pNN50"
    "pNN20"
    "CV_RR"
    "SD1"
    "SD2"
    "SD1_SD2_ratio"
    "SampEn"
    "LF"
    "HF"
    "LF_HF"
    "LFnu"
    "HFnu"
    "DF_completo_Hz"
    "DF_residual_Hz"
];

%% 3) Comprobar que existan en la tabla

biomarcadores = biomarcadores(ismember(biomarcadores, string(Tsub.Properties.VariableNames)));

if numel(biomarcadores) < 2
    error('Hay menos de 2 biomarcadores de FA presentes en la tabla maestra.');
end

%% 4) Crear matriz numérica

X = Tsub{:, cellstr(biomarcadores)};

fprintf('\nValores válidos por biomarcador:\n');

for i = 1:numel(biomarcadores)

    n_validos = sum(isfinite(X(:, i)));

    fprintf('  %-30s %d de %d\n', ...
        char(biomarcadores(i)), ...
        n_validos, ...
        size(X, 1));
end

idx_ok = all(isfinite(X), 2);

fprintf('\nPacientes con todos los biomarcadores completos: %d de %d\n', ...
    sum(idx_ok), ...
    size(X, 1));

X = X(idx_ok, :);

if size(X, 1) < 3
    error('Hay menos de 3 observaciones válidas.');
end

fprintf('Observaciones usadas: %d\n', size(X, 1));
fprintf('Biomarcadores usados: %d\n', numel(biomarcadores));
disp(biomarcadores)

%% 5) Calcular correlación de Spearman

C = corr(X, 'Type', 'Spearman', 'Rows', 'pairwise');

%% 6) Crear etiquetas más limpias

etiquetas = strrep(cellstr(biomarcadores), '_', ' ');

%% 7) Ajustar tamaño de figura según número de biomarcadores

n = numel(biomarcadores);

if n > 20
    ancho = 1900;
    alto = 1500;
    fuente = 8;
elseif n > 12
    ancho = 1500;
    alto = 1100;
    fuente = 10;
else
    ancho = 950;
    alto = 800;
    fuente = 12;
end

%% 8) Dibujar heatmap

fig = figure('Color', 'w', 'Position', [100 100 ancho alto]);

h = heatmap(C);

h.XDisplayLabels = etiquetas;
h.YDisplayLabels = etiquetas;

h.Title = ['Correlación entre biomarcadores para UMAP - ' ...
    strrep(char(nombre_comp), '_', ' ')];

h.XLabel = 'Biomarcadores';
h.YLabel = 'Biomarcadores';
h.Colormap = parula;
h.ColorLimits = [-1 1];
h.CellLabelFormat = '%.2f';
h.FontSize = fuente;

%% 9) Guardar solo el PNG

nombre_base = 'Heatmap_FA_PAROXISTICA_FA_vs_FA_PERSISTENTE';

archivo_png = fullfile(carpeta_salida, [nombre_base '.png']);

exportgraphics(fig, archivo_png, 'Resolution', 300);

close(fig)

%% 10) Mensaje final

fprintf('\nHeatmap guardado:\n');
fprintf('- PNG: %s\n', archivo_png);

fprintf('\nFIN. Heatmap PNG generado.\n');