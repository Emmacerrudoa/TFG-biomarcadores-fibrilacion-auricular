clear
clc
close all

%% ============================================================
% HEATMAP DE CORRELACIÓN PARA UMAP
% TRANSICIONES FA -> RS
% VENTANAS DE 10 s, 3 MOMENTOS
%
% Momentos:
%   MOMENTO 1: -180 a -170 s
%   MOMENTO 2:  -90 a  -80 s
%   MOMENTO 3:  -10 a    0 s
%
% Se utilizan todos los biomarcadores disponibles en ventanas de FA:
%
% Biomarcadores RR:
%   RR_mean
%   SDNN
%   RMSSD
%   SDSD
%   pNN20
%   pNN50
%   CV_RR
%   SD1
%   SD2
%   SD1_SD2
%
% Frecuencia dominante:
%   DF_completo_Hz
%   DF_residual_Hz
%
% No se seleccionan biomarcadores por significación estadística,
% porque en FA -> RS solo pocos biomarcadores fueron significativos
% tras FDR y el heatmap perdería utilidad.
%
% Solo se conservan transiciones completas con:
%   - una fila para el momento 1
%   - una fila para el momento 2
%   - una fila para el momento 3
%
% Si una de las tres ventanas tiene NaN o Inf, se elimina también
% el resto de ventanas de la misma transición.
%
% La correlación se calcula con Spearman usando las tres ventanas.
% Solo se guarda el heatmap en formato PNG.
%% ============================================================

%% ============================
% RUTAS
%% ============================

carpeta_base = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_FA_A_RS_3ventanas';

archivo_biomarcadores = fullfile( ...
    carpeta_base, ...
    'biomarcadores_transiciones_FA_a_RS_10s_3ventanas.xlsx');

archivo_png = fullfile( ...
    carpeta_base, ...
    'Heatmap_todos_biomarcadores_FA_a_RS_10s_3ventanas.png');

%% ============================
% COMPROBAR ARCHIVO
%% ============================

if ~isfile(archivo_biomarcadores)

    error('No se encuentra el archivo de biomarcadores:\n%s', ...
        archivo_biomarcadores);
end

%% ============================
% LEER TABLA
%% ============================

T = readtable(archivo_biomarcadores);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

%% ============================
% COMPROBAR COLUMNAS OBLIGATORIAS
%% ============================

columnas_obligatorias = { ...
    'ID_transicion', ...
    'Momento', ...
    'Orden_momento'};

for i = 1:numel(columnas_obligatorias)

    if ~ismember(columnas_obligatorias{i}, ...
            T.Properties.VariableNames)

        error('Falta la columna %s en la tabla.', ...
            columnas_obligatorias{i});
    end
end

T.ID_transicion = strtrim(string(T.ID_transicion));
T.Momento = strtrim(string(T.Momento));

%% ============================
% BIOMARCADORES
%% ============================

biomarcadores = { ...
    'RR_mean', ...
    'SDNN', ...
    'RMSSD', ...
    'SDSD', ...
    'pNN20', ...
    'pNN50', ...
    'CV_RR', ...
    'SD1', ...
    'SD2', ...
    'SD1_SD2', ...
    'DF_completo_Hz', ...
    'DF_residual_Hz' ...
};

%% ============================
% COMPROBAR BIOMARCADORES
%% ============================

for i = 1:numel(biomarcadores)

    biom = biomarcadores{i};

    if ~ismember(biom, T.Properties.VariableNames)

        error('No se encuentra el biomarcador %s en la tabla.', biom);
    end

    if ~isnumeric(T.(biom))

        error('La columna %s no es numérica.', biom);
    end
end

fprintf('\nBiomarcadores incluidos: %d\n', ...
    numel(biomarcadores));

disp(string(biomarcadores(:)))

%% ============================
% CONSERVAR SOLO TRANSICIONES
% CON LOS TRES MOMENTOS COMPLETOS
%% ============================

transiciones = unique(T.ID_transicion, 'stable');
transiciones_validas = strings(0,1);

for i = 1:numel(transiciones)

    id = transiciones(i);

    idx_id = T.ID_transicion == id;

    n_momento1 = sum(idx_id & T.Orden_momento == 1);
    n_momento2 = sum(idx_id & T.Orden_momento == 2);
    n_momento3 = sum(idx_id & T.Orden_momento == 3);

    if n_momento1 == 1 && ...
            n_momento2 == 1 && ...
            n_momento3 == 1

        transiciones_validas(end+1,1) = id; %#ok<SAGROW>
    end
end

if isempty(transiciones_validas)

    error('No se encontraron transiciones con los tres momentos completos.');
end

Tsub = T(ismember(T.ID_transicion, transiciones_validas), :);

fprintf('\nTransiciones completas antes de eliminar NaN/Inf: %d\n', ...
    numel(transiciones_validas));

%% ============================
% CREAR MATRIZ NUMÉRICA
%% ============================

X = Tsub{:, biomarcadores};

%% ============================
% ELIMINAR FILAS CON NaN O Inf
%% ============================

idx_filas_validas = all(isfinite(X), 2);

Tsub = Tsub(idx_filas_validas, :);

%% ============================
% VOLVER A EXIGIR TRANSICIONES COMPLETAS
%
% Si una de las tres ventanas tiene NaN o Inf, se elimina también
% el resto de ventanas de la misma transición.
%% ============================

transiciones_despues_nan = unique(Tsub.ID_transicion, 'stable');
transiciones_finales = strings(0,1);

for i = 1:numel(transiciones_despues_nan)

    id = transiciones_despues_nan(i);

    idx_id = Tsub.ID_transicion == id;

    n_momento1 = sum(idx_id & Tsub.Orden_momento == 1);
    n_momento2 = sum(idx_id & Tsub.Orden_momento == 2);
    n_momento3 = sum(idx_id & Tsub.Orden_momento == 3);

    if n_momento1 == 1 && ...
            n_momento2 == 1 && ...
            n_momento3 == 1

        transiciones_finales(end+1,1) = id; %#ok<SAGROW>
    end
end

if isempty(transiciones_finales)

    error(['No quedan transiciones completas después de eliminar ' ...
        'NaN o Inf.']);
end

idx_final = ismember(Tsub.ID_transicion, transiciones_finales);

Tsub = Tsub(idx_final, :);

Tsub = sortrows( ...
    Tsub, ...
    {'ID_transicion', 'Orden_momento'});

X = Tsub{:, biomarcadores};

if size(X,1) < 3

    error('Hay menos de tres observaciones válidas para calcular correlaciones.');
end

fprintf('Transiciones completas finales: %d\n', ...
    numel(transiciones_finales));

fprintf('Observaciones utilizadas: %d\n', ...
    size(X,1));

fprintf('FA alejada: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 1));

fprintf('FA intermedia: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 2));

fprintf('FA previa a RS: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 3));

%% ============================
% COMPROBAR VARIABILIDAD
%% ============================

desviaciones = std(X, 0, 1, 'omitnan');

if any(desviaciones == 0)

    columnas_constantes = ...
        biomarcadores(desviaciones == 0);

    error('Hay biomarcadores sin variabilidad: %s', ...
        strjoin(columnas_constantes, ', '));
end

%% ============================
% CORRELACIÓN DE SPEARMAN
%% ============================

C = corr( ...
    X, ...
    'Type', 'Spearman', ...
    'Rows', 'pairwise');

%% ============================
% ETIQUETAS DE LA FIGURA
%% ============================

etiquetas = { ...
    'RR mean', ...
    'SDNN', ...
    'RMSSD', ...
    'SDSD', ...
    'pNN20', ...
    'pNN50', ...
    'CV RR', ...
    'SD1', ...
    'SD2', ...
    'SD1/SD2', ...
    'DF completo Hz', ...
    'DF residual Hz' ...
};

n = numel(biomarcadores);

if n > 18

    ancho = 1800;
    alto = 1450;
    fuente = 8;

elseif n > 12

    ancho = 1500;
    alto = 1150;
    fuente = 9;

else

    ancho = 1150;
    alto = 900;
    fuente = 11;
end

%% ============================
% DIBUJAR HEATMAP
%% ============================

fig = figure( ...
    'Color', 'w', ...
    'Position', [100 100 ancho alto]);

h = heatmap( ...
    etiquetas, ...
    etiquetas, ...
    C);

h.Title = ...
    'Correlación entre biomarcadores en transiciones FA a RS';

h.XLabel = 'Biomarcadores';
h.YLabel = 'Biomarcadores';

h.Colormap = parula;
h.ColorLimits = [-1 1];
h.CellLabelFormat = '%.2f';
h.FontSize = fuente;

%% ============================
% GUARDAR PNG
%% ============================

drawnow

frame = getframe(fig);

imwrite(frame.cdata, archivo_png);

close(fig)

fprintf('\nHeatmap guardado en:\n%s\n', archivo_png);
fprintf('\nFIN.\n');