clear
clc
close all

%% ============================================================
% HEATMAP DE CORRELACIÓN PARA UMAP
% FA PERSISTENTE
%
% Se utilizan todos los biomarcadores disponibles en ventanas de 60 s:
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
%   DF_completo_Hz
%   DF_residual_Hz
%
% No se seleccionan por significación estadística, porque FA persistente
% se utiliza como control de estabilidad temporal.
%
% Solo se conservan pacientes con:
%   - una fila para el minuto 1
%   - una fila para el minuto 3
%
% La correlación se calcula con Spearman usando las dos ventanas.
% Solo se guarda el heatmap en formato PNG.
%% ============================================================

%% ============================
% RUTAS
%% ============================

carpeta_base = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_FA_PERSISTENTE';

archivo_biomarcadores = fullfile( ...
    carpeta_base, ...
    'biomarcadores_FA_persistente_60s.xlsx');

archivo_png = fullfile( ...
    carpeta_base, ...
    'Heatmap_correlacion_biomarcadores_FA_persistente.png');

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

columnas_obligatorias = { ...
    'Paciente', ...
    'Momento', ...
    'Orden_momento'};

for i = 1:numel(columnas_obligatorias)

    if ~ismember(columnas_obligatorias{i}, ...
            T.Properties.VariableNames)

        error('Falta la columna %s en la tabla.', ...
            columnas_obligatorias{i});
    end
end

T.Paciente = strtrim(string(T.Paciente));
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

        error('La columna %s no es numerica.', biom);
    end
end

fprintf('\nBiomarcadores incluidos: %d\n', ...
    numel(biomarcadores));

disp(string(biomarcadores(:)))

%% ============================
% CONSERVAR SOLO PACIENTES
% CON LOS DOS MOMENTOS COMPLETOS
%% ============================

pacientes = unique(T.Paciente, 'stable');
pacientes_validos = strings(0,1);

for i = 1:numel(pacientes)

    paciente = pacientes(i);

    idx_paciente = T.Paciente == paciente;

    n_minuto1 = sum(idx_paciente & T.Orden_momento == 1);
    n_minuto3 = sum(idx_paciente & T.Orden_momento == 2);

    if n_minuto1 == 1 && n_minuto3 == 1

        pacientes_validos(end+1,1) = paciente; %#ok<SAGROW>
    end
end

if isempty(pacientes_validos)

    error('No se encontraron pacientes con los dos momentos completos.');
end

Tsub = T(ismember(T.Paciente, pacientes_validos), :);

fprintf('\nPacientes completos antes de eliminar NaN/Inf: %d\n', ...
    numel(pacientes_validos));

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
% VOLVER A EXIGIR PAREJAS COMPLETAS
%
% Si una de las dos ventanas tiene NaN o Inf, se elimina también
% la otra ventana del mismo paciente.
%% ============================

pacientes_despues_nan = unique(Tsub.Paciente, 'stable');
pacientes_finales = strings(0,1);

for i = 1:numel(pacientes_despues_nan)

    paciente = pacientes_despues_nan(i);

    idx_paciente = Tsub.Paciente == paciente;

    n_minuto1 = sum(idx_paciente & Tsub.Orden_momento == 1);
    n_minuto3 = sum(idx_paciente & Tsub.Orden_momento == 2);

    if n_minuto1 == 1 && n_minuto3 == 1

        pacientes_finales(end+1,1) = paciente; %#ok<SAGROW>
    end
end

idx_final = ismember(Tsub.Paciente, pacientes_finales);

Tsub = Tsub(idx_final, :);

X = Tsub{:, biomarcadores};

if size(X,1) < 3

    error('Hay menos de tres observaciones validas para calcular correlaciones.');
end

fprintf('Pacientes completos finales: %d\n', ...
    numel(pacientes_finales));

fprintf('Observaciones utilizadas: %d\n', ...
    size(X,1));

fprintf('Minuto 1: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 1));

fprintf('Minuto 3: %d observaciones\n', ...
    sum(Tsub.Orden_momento == 2));

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

if n > 12

    ancho = 1500;
    alto = 1100;
    fuente = 10;

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
    'Correlación entre biomarcadores en FA persistente';

h.XLabel = 'Biomarcadores';
h.YLabel = 'Biomarcadores';

h.Colormap = parula;
h.ColorLimits = [-1 1];
h.CellLabelFormat = '%.2f';
h.FontSize = fuente;

%% ============================
% GUARDAR PNG
%% ============================

exportgraphics( ...
    fig, ...
    archivo_png, ...
    'Resolution', 300);

close(fig)

fprintf('\nHeatmap guardado en:\n%s\n', archivo_png);
fprintf('\nFIN.\n');