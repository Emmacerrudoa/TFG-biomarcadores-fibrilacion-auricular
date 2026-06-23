clear; 
clc; 
close all;

%% ============================================================
% HEATMAP DE CORRELACIÓN ENTRE BIOMARCADORES EN RS
% Comparación: SANO vs FA PAROXÍSTICA EN RITMO SINUSAL
%
% Selecciona automáticamente los biomarcadores marcados como Usar_UMAP = Sí
% en la tabla de resultados de Mann-Whitney.
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

archivo_mannwhitney = fullfile(carpeta_salida, ...
    'resultados_MannWhitney_biomarcadores_correcion.xlsx');

%% LEER TABLAS

T = readtable(archivo_tabla);
Rmw = readtable(archivo_mannwhitney);

T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
Rmw.Properties.VariableNames = matlab.lang.makeValidName(Rmw.Properties.VariableNames);

T.Grupo = strtrim(string(T.Grupo));

Rmw.Comparacion = strtrim(string(Rmw.Comparacion));
Rmw.Biomarcador = strtrim(string(Rmw.Biomarcador));
Rmw.Usar_UMAP = strtrim(string(Rmw.Usar_UMAP));

%% DEFINIR COMPARACIÓN

nombre_comp = "SANO_vs_FA_PAROXISTICA_RS";
grupos = ["SANO", "FA_PAROXISTICA_RS"];

fprintf('\n============================================\n');
fprintf('Comparación: %s\n', nombre_comp);
fprintf('============================================\n');

%% 1) Seleccionar biomarcadores marcados para UMAP

idx_umap = Rmw.Comparacion == nombre_comp & ...
          (upper(Rmw.Usar_UMAP) == "SÍ" | upper(Rmw.Usar_UMAP) == "SI");

biomarcadores = unique(Rmw.Biomarcador(idx_umap), 'stable');

if numel(biomarcadores) < 2
    error('No hay suficientes biomarcadores para hacer heatmap.');
end

%% 2) Filtrar filas de los grupos correspondientes

idx_grupo = ismember(T.Grupo, grupos);
Tsub = T(idx_grupo, :);

%% 3) Quedarse solo con biomarcadores que existan en la tabla maestra

biomarcadores = biomarcadores(ismember(biomarcadores, string(Tsub.Properties.VariableNames)));

if numel(biomarcadores) < 2
    error('Hay menos de 2 biomarcadores presentes en la tabla maestra.');
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

if n > 12
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

nombre_base = 'Heatmap_SANO_vs_FA_PAROXISTICA_RS';

archivo_png = fullfile(carpeta_salida, [nombre_base '.png']);

exportgraphics(fig, archivo_png, 'Resolution', 300);

close(fig)

%% 10) Mensaje final

fprintf('\nHeatmap guardado:\n');
fprintf('- PNG: %s\n', archivo_png);

fprintf('\nFIN. Heatmap PNG generado.\n');