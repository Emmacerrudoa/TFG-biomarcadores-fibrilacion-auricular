clear
clc
close all

%% ============================================================
% UMAP 2D:
% SUJETOS SANOS
%
% Comparación:
%   1. Minuto 1:   0 a  60 s
%   2. Minuto 3: 120 a 180 s
%
% Biomarcadores seleccionados:
%   1. RR_mean
%   2. RMSSD
%   3. pNN20
%   4. CV_RR
%   5. SD2
%   6. SD1_SD2
%   7. DF_residual_Hz
%   8. P_CorrIntraMedia
%   9. P_AmpMedia
%  10. P_StdMedia
%  11. T_CorrIntraMedia
%  12. T_AmpMedia
%  13. T_StdMedia
%% ============================================================

rng(42)

%% RUTA DE LA TABLA

archivo_tabla = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_SANO\biomarcadores_sanos_60s.xlsx';

carpeta_salida = fileparts(archivo_tabla);

%% LEER TABLA

if ~isfile(archivo_tabla)
    error('No se encuentra el archivo:\n%s', archivo_tabla);
end

T = readtable(archivo_tabla);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

T.Paciente = strtrim(string(T.Paciente));
T.Momento = strtrim(string(T.Momento));

%% COMPROBAR COLUMNAS OBLIGATORIAS

columnas_obligatorias = { ...
    'Paciente', ...
    'Momento', ...
    'Orden_momento'};

for i = 1:numel(columnas_obligatorias)

    if ~ismember(columnas_obligatorias{i}, ...
            T.Properties.VariableNames)

        error('No se encuentra la columna obligatoria: %s', ...
            columnas_obligatorias{i});
    end
end

%% BIOMARCADORES SELECCIONADOS

biomarcadores = { ...
    'RR_mean', ...
    'RMSSD', ...
    'pNN20', ...
    'CV_RR', ...
    'SD2', ...
    'SD1_SD2', ...
    'DF_residual_Hz', ...
    'P_CorrIntraMedia', ...
    'P_AmpMedia', ...
    'P_StdMedia', ...
    'T_CorrIntraMedia', ...
    'T_AmpMedia', ...
    'T_StdMedia'};

%% COMPROBAR QUE LOS BIOMARCADORES EXISTEN

for i = 1:numel(biomarcadores)

    if ~ismember(biomarcadores{i}, ...
            T.Properties.VariableNames)

        error('No se encuentra la columna: %s', ...
            biomarcadores{i});
    end
end

%% CONSERVAR SOLO PACIENTES CON LOS DOS MOMENTOS

pacientes = unique(T.Paciente, 'stable');

pacientes_validos = strings(0,1);

for i = 1:numel(pacientes)

    paciente_actual = pacientes(i);

    idx_paciente = T.Paciente == paciente_actual;

    n_minuto1 = sum(idx_paciente & T.Orden_momento == 1);
    n_minuto3 = sum(idx_paciente & T.Orden_momento == 2);

    if n_minuto1 == 1 && n_minuto3 == 1

        pacientes_validos(end+1,1) = ...
            paciente_actual; %#ok<SAGROW>
    end
end

if isempty(pacientes_validos)

    error('No se encontraron pacientes con los dos momentos completos.');
end

T_umap = T(ismember(T.Paciente, pacientes_validos), :);

%% CREAR MATRIZ DE DATOS

X = T_umap{:, biomarcadores};

momento = T_umap.Orden_momento;
paciente = T_umap.Paciente;

%% ELIMINAR PAREJAS CON NaN O Inf

filas_validas = all(isfinite(X), 2);

pacientes_con_datos_validos = unique( ...
    paciente(filas_validas), ...
    'stable');

pacientes_finales = strings(0,1);

for i = 1:numel(pacientes_con_datos_validos)

    paciente_actual = pacientes_con_datos_validos(i);

    idx_paciente = ...
        paciente == paciente_actual & filas_validas;

    n_minuto1 = sum(idx_paciente & momento == 1);
    n_minuto3 = sum(idx_paciente & momento == 2);

    if n_minuto1 == 1 && n_minuto3 == 1

        pacientes_finales(end+1,1) = ...
            paciente_actual; %#ok<SAGROW>
    end
end

idx_final = ...
    ismember(paciente, pacientes_finales) & ...
    filas_validas;

X = X(idx_final, :);
momento = momento(idx_final);
paciente = paciente(idx_final);

if isempty(X)

    error('No quedan pacientes con las dos ventanas válidas.');
end

%% ÍNDICES DE LOS DOS MOMENTOS

idx1 = momento == 1;
idx2 = momento == 2;

%% MOSTRAR NÚMERO DE OBSERVACIONES

fprintf('\nPacientes incluidos en el UMAP:\n');

fprintf('Minuto 1: %d\n', sum(idx1));
fprintf('Minuto 3: %d\n', sum(idx2));

fprintf('Pacientes completos: %d\n', ...
    numel(pacientes_finales));

fprintf('Total de puntos: %d\n\n', ...
    size(X,1));

%% COMPROBAR VARIABILIDAD

desviaciones = std(X, 0, 1);

if any(desviaciones == 0)

    columnas_constantes = ...
        biomarcadores(desviaciones == 0);

    error('Hay biomarcadores sin variabilidad: %s', ...
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

%% REPRESENTACIÓN

f = figure( ...
    'Color', 'w', ...
    'Position', [100 100 900 700]);

hold on

scatter( ...
    reduction(idx1,1), ...
    reduction(idx1,2), ...
    85, ...
    [0.00 0.45 0.74], ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.7);

scatter( ...
    reduction(idx2,1), ...
    reduction(idx2,2), ...
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
    'UMAP temporal en sujetos sanos', ...
    'FontSize', 17, ...
    'FontWeight', 'normal');

legend( ...
    { ...
    'Minuto 1 (0 a 60 s)', ...
    'Minuto 3 (120 a 180 s)'}, ...
    'Location', 'best', ...
    'FontSize', 13);

set(gca, ...
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
    'UMAP_2D_SANO_MINUTO_1_vs_MINUTO_3.png');

exportgraphics( ...
    f, ...
    archivo_salida, ...
    'Resolution', 300);

fprintf('Figura guardada en:\n%s\n', ...
    archivo_salida);
