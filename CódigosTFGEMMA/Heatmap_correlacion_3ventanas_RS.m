clear
clc
close all

%% ============================================================
% HEATMAP DE CORRELACIÓN PARA UMAP
% TRANSICIONES RS -> FA
% VENTANAS DE 10 s, 3 MOMENTOS
%
% Momentos:
%   MOMENTO 1: -180 a -170 s
%   MOMENTO 2:  -90 a  -80 s
%   MOMENTO 3:  -10 a    0 s
%
% Se seleccionan únicamente los biomarcadores significativos
% después de la corrección FDR en comparaciones que incluyen
% el momento 3:
%
%   - Momento 3 - Momento 1
%   - Momento 3 - Momento 2
%
% No se usan los biomarcadores significativos solo en:
%
%   - Momento 2 - Momento 1
%
% porque esa comparación no incluye la ventana inmediatamente
% previa al inicio de la FA.
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
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_RS_A_FA_3ventanas';

archivo_biomarcadores = fullfile( ...
    carpeta_base, ...
    'biomarcadores_transiciones_RS_a_FA_10s_3ventanas.xlsx');

archivo_wilcoxon = fullfile( ...
    carpeta_base, ...
    'resultados_Wilcoxon_FA_paroxistica_en_RS_10s_3ventanas_FDR.xlsx');

archivo_png = fullfile( ...
    carpeta_base, ...
    'Heatmap_biomarcadores_significativos_FDR_RS_a_FA_10s_3ventanas.png');

%% ============================
% COMPROBAR ARCHIVOS
%% ============================

if ~isfile(archivo_biomarcadores)

    error('No se encuentra el archivo de biomarcadores:\n%s', ...
        archivo_biomarcadores);
end

if ~isfile(archivo_wilcoxon)

    error('No se encuentra el archivo de Wilcoxon:\n%s', ...
        archivo_wilcoxon);
end

%% ============================
% LEER TABLAS
%% ============================

T = readtable(archivo_biomarcadores);
R = readtable(archivo_wilcoxon);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

R.Properties.VariableNames = ...
    matlab.lang.makeValidName(R.Properties.VariableNames);

%% ============================
% COMPROBAR COLUMNAS OBLIGATORIAS
%% ============================

columnas_T = { ...
    'ID_transicion', ...
    'Momento', ...
    'Orden_momento'};

for i = 1:numel(columnas_T)

    if ~ismember(columnas_T{i}, ...
            T.Properties.VariableNames)

        error('Falta la columna %s en la tabla de biomarcadores.', ...
            columnas_T{i});
    end
end

columnas_R = { ...
    'Biomarcador', ...
    'Comparacion', ...
    'p_ajustado_FDR', ...
    'Usar_UMAP_FDR'};

for i = 1:numel(columnas_R)

    if ~ismember(columnas_R{i}, ...
            R.Properties.VariableNames)

        error('Falta la columna %s en la tabla estadística.', ...
            columnas_R{i});
    end
end

%% ============================
% NORMALIZAR FORMATOS
%% ============================

T.ID_transicion = ...
    strtrim(string(T.ID_transicion));

T.Momento = ...
    strtrim(string(T.Momento));

R.Biomarcador = ...
    strtrim(string(R.Biomarcador));

R.Comparacion = ...
    strtrim(string(R.Comparacion));

R.Usar_UMAP_FDR = ...
    upper(strtrim(string(R.Usar_UMAP_FDR)));

%% ============================
% SELECCIONAR BIOMARCADORES
% SIGNIFICATIVOS TRAS FDR
%% ============================

idx_fdr = ...
    (R.Usar_UMAP_FDR == "SÍ" | ...
     R.Usar_UMAP_FDR == "SI") & ...
    (R.Comparacion == "Momento 3 - Momento 1" | ...
     R.Comparacion == "Momento 3 - Momento 2");

R_fdr = R(idx_fdr, :);

if isempty(R_fdr)

    error(['No hay biomarcadores significativos tras FDR ' ...
        'en comparaciones que incluyan el momento 3.']);
end

R_fdr = sortrows( ...
    R_fdr, ...
    'p_ajustado_FDR', ...
    'ascend');

biomarcadores = ...
    unique(string(R_fdr.Biomarcador), 'stable');

%% ============================
% COMPROBAR BIOMARCADORES
%% ============================

existen = ismember( ...
    biomarcadores, ...
    string(T.Properties.VariableNames));

if any(~existen)

    fprintf('\nBiomarcadores significativos que no aparecen en la tabla:\n');
    disp(biomarcadores(~existen))
end

biomarcadores = biomarcadores(existen);

if isempty(biomarcadores)

    error(['Ningún biomarcador significativo tras FDR ' ...
        'aparece en la tabla de biomarcadores.']);
end

for i = 1:numel(biomarcadores)

    biom = char(biomarcadores(i));

    if ~isnumeric(T.(biom))

        error('La columna %s no es numérica.', biom);
    end
end

fprintf('\nBiomarcadores incluidos: %d\n', ...
    numel(biomarcadores));

disp(biomarcadores)

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

X = Tsub{:, cellstr(biomarcadores)};

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

X = Tsub{:, cellstr(biomarcadores)};

if size(X,1) < 3

    error('Hay menos de tres observaciones válidas para calcular correlaciones.');
end

fprintf('Transiciones completas finales: %d\n', ...
    numel(transiciones_finales));

fprintf('Observaciones utilizadas: %d\n', ...
    size(X,1));

fprintf('RS alejado: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 1));

fprintf('RS intermedio: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 2));

fprintf('RS previo: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 3));

%% ============================
% COMPROBAR VARIABILIDAD
%% ============================

desviaciones = std(X, 0, 1, 'omitnan');

if any(desviaciones == 0)

    columnas_constantes = ...
        biomarcadores(desviaciones == 0);

    error('Hay biomarcadores sin variabilidad: %s', ...
        strjoin(cellstr(columnas_constantes), ', '));
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

etiquetas = crear_etiquetas_biomarcadores(biomarcadores);

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
    'Correlación entre biomarcadores significativos en transiciones RS a FA';

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

%% ============================================================
% FUNCIÓN PARA CREAR ETIQUETAS LEGIBLES
%% ============================================================

function etiquetas = crear_etiquetas_biomarcadores(biomarcadores)

etiquetas = strings(size(biomarcadores));

for i = 1:numel(biomarcadores)

    biom = string(biomarcadores(i));

    switch biom

        case "RR_mean"
            etiquetas(i) = "RR mean";

        case "SDNN"
            etiquetas(i) = "SDNN";

        case "RMSSD"
            etiquetas(i) = "RMSSD";

        case "SDSD"
            etiquetas(i) = "SDSD";

        case "pNN20"
            etiquetas(i) = "pNN20";

        case "pNN50"
            etiquetas(i) = "pNN50";

        case "CV_RR"
            etiquetas(i) = "CV RR";

        case "SD1"
            etiquetas(i) = "SD1";

        case "SD2"
            etiquetas(i) = "SD2";

        case "SD1_SD2"
            etiquetas(i) = "SD1/SD2";

        case "DF_completo_Hz"
            etiquetas(i) = "DF completo Hz";

        case "DF_residual_Hz"
            etiquetas(i) = "DF residual Hz";

        case "P_CorrIntraMedia"
            etiquetas(i) = "P CorrIntraMedia";

        case "P_CorrIntraStd"
            etiquetas(i) = "P CorrIntraStd";

        case "P_AmpMedia"
            etiquetas(i) = "P AmpMedia";

        case "P_AmpStd"
            etiquetas(i) = "P AmpStd";

        case "P_StdMedia"
            etiquetas(i) = "P StdMedia";

        case "T_CorrIntraMedia"
            etiquetas(i) = "T CorrIntraMedia";

        case "T_CorrIntraStd"
            etiquetas(i) = "T CorrIntraStd";

        case "T_AmpMedia"
            etiquetas(i) = "T AmpMedia";

        case "T_AmpStd"
            etiquetas(i) = "T AmpStd";

        case "T_StdMedia"
            etiquetas(i) = "T StdMedia";

        otherwise

            etiquetas(i) = strrep(biom, "_", " ");
    end
end

etiquetas = cellstr(etiquetas);

end