clear
clc
close all

%% ============================================================
% CREAR TABLA MAESTRA DE BIOMARCADORES POR PACIENTE
%
% Incluye:
%   - Biomarcadores RR
%   - Frecuencia dominante DF por paciente
%   - Biomarcadores de onda P
%   - Biomarcadores de onda T
%
% Las tablas de entrada ya contienen una fila por paciente.
%
% No incluye:
%   - Columnas de control
%   - Potencia del pico de DF
%   - Desviación intrapaciente de DF
%   - Área y energía de onda P/T
%   - Correlaciones interpaciente
%   - Resultados estadísticos
%
% Salidas:
%   - tabla_maestra_biomarcadores_TFG.xlsx
%   - resumen_grupos_tabla_maestra.xlsx
%% ============================================================

%% ============================
% RUTAS
%% ============================

% Modificar esta ruta según la ubicación local de las carpetas de resultados.
ruta_base = 'C:\Users\Emma\Documents\MATLAB';

carpeta_RR = fullfile(ruta_base, 'analisisRR_hibrido');

carpeta_DF = fullfile(ruta_base, ...
    'DF_120s_completo_y_residual');

carpeta_P = fullfile(ruta_base, ...
    'analisis_ondaPHIBRIDOnuevo');

carpeta_T = fullfile(ruta_base, ...
    'analisis_ondaTHIBRIDOnuevo');

archivo_RR = fullfile(carpeta_RR, ...
    'biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia_y_corr.xlsx');

archivo_DF = fullfile(carpeta_DF, ...
    'resultados_DF_por_paciente.xlsx');

archivo_P = fullfile(carpeta_P, ...
    'features_ondaP.xlsx');

archivo_T = fullfile(carpeta_T, ...
    'features_ondaT.xlsx');

carpeta_salida = fullfile(ruta_base, ...
    'RESULTADOS_ESTADISTICOS_TFG');

if ~exist(carpeta_salida, 'dir')
    mkdir(carpeta_salida);
end

archivo_tabla_maestra = fullfile(carpeta_salida, ...
    'tabla_maestra_biomarcadores_TFG.xlsx');

archivo_resumen_grupos = fullfile(carpeta_salida, ...
    'resumen_grupos_tabla_maestra.xlsx');

%% ============================
% LEER TABLAS
%% ============================

fprintf('\nLeyendo tablas...\n');

T_RR = leer_y_preparar_tabla(archivo_RR, "RR");
T_DF = leer_y_preparar_tabla(archivo_DF, "DF");
T_P  = leer_y_preparar_tabla(archivo_P,  "P");
T_T  = leer_y_preparar_tabla(archivo_T,  "T");

%% ============================
% PREPARAR RR
%% ============================

fprintf('\nPreparando RR...\n');

cols_control_RR = { ...
    'NumeroRR', ...
    'NumeroVentanasValidas', ...
    'NumeroSegmentos', ...
    'Numero_dRR', ...
    'CorrelacionHistRRGrupo' ...
};

T_RR = quitar_columnas_si_existen(T_RR, cols_control_RR);

%% ============================
% PREPARAR DF POR PACIENTE
%% ============================

fprintf('\nPreparando DF por Paciente...\n');

% En la tabla de DF, la columna Registro contiene el identificador
% del paciente y se renombra como Paciente para poder unir las tablas.
cols_DF = { ...
    'Registro', ...
    'Grupo', ...
    'DF_completo_media_Hz', ...
    'DF_residual_media_Hz' ...
};

T_DF = T_DF(:, cols_DF);

T_DF.Properties.VariableNames{ ...
    strcmp(T_DF.Properties.VariableNames, 'Registro')} = ...
    'Paciente';

T_DF.Paciente = strtrim(string(T_DF.Paciente));
T_DF.Grupo = strtrim(string(T_DF.Grupo));

% Nombres finales para la tabla maestra.

T_DF.Properties.VariableNames{ ...
    strcmp(T_DF.Properties.VariableNames, ...
    'DF_completo_media_Hz')} = ...
    'DF_completo_Hz';

T_DF.Properties.VariableNames{ ...
    strcmp(T_DF.Properties.VariableNames, ...
    'DF_residual_media_Hz')} = ...
    'DF_residual_Hz';

%% ============================
% PREPARAR ONDA P Y ONDA T
%% ============================

fprintf('\nPreparando onda P y onda T...\n');

cols_control_P = { ...
    'NumOndasP', ...
    'N_ventanas_totales', ...
    'N_ventanas_usadas', ...
    'N_descartadas_sin_R', ...
    'N_descartadas_pocas_R', ...
    'N_descartadas_senal' ...
};

cols_control_T = { ...
    'NumOndasT', ...
    'N_ventanas_totales', ...
    'N_ventanas_usadas', ...
    'N_descartadas_sin_R', ...
    'N_descartadas_pocas_R', ...
    'N_descartadas_senal' ...
};

T_P = quitar_columnas_si_existen(T_P, cols_control_P);
T_T = quitar_columnas_si_existen(T_T, cols_control_T);

%% ============================
% AÑADIR PREFIJOS
%% ============================

% RR se deja sin prefijo.
% Las columnas DF ya empiezan por DF_.
% Las columnas de onda P y onda T se identifican con los prefijos P_ y T_.

T_P = anadir_prefijo_biomarcadores(T_P, 'P_');
T_T = anadir_prefijo_biomarcadores(T_T, 'T_');

%% ============================
% CONTROL DE DUPLICADOS
%% ============================

claves_RR = strcat(T_RR.Paciente, "|", T_RR.Grupo);
claves_DF = strcat(T_DF.Paciente, "|", T_DF.Grupo);
claves_P  = strcat(T_P.Paciente,  "|", T_P.Grupo);
claves_T  = strcat(T_T.Paciente,  "|", T_T.Grupo);

if numel(unique(claves_RR)) ~= height(T_RR)
    error('Hay pacientes duplicados en la tabla RR.');
end

if numel(unique(claves_DF)) ~= height(T_DF)
    error('Hay pacientes duplicados en la tabla DF.');
end

if numel(unique(claves_P)) ~= height(T_P)
    error('Hay pacientes duplicados en la tabla de onda P.');
end

if numel(unique(claves_T)) ~= height(T_T)
    error('Hay pacientes duplicados en la tabla de onda T.');
end

fprintf('\nControl de duplicados superado.\n');

%% ============================
% UNIR TABLAS
%% ============================

fprintf('\nUniendo tablas...\n');

T_maestra = T_RR;

T_maestra = outerjoin(T_maestra, T_DF, ...
    'Keys', {'Paciente', 'Grupo'}, ...
    'MergeKeys', true);

T_maestra = outerjoin(T_maestra, T_P, ...
    'Keys', {'Paciente', 'Grupo'}, ...
    'MergeKeys', true);

T_maestra = outerjoin(T_maestra, T_T, ...
    'Keys', {'Paciente', 'Grupo'}, ...
    'MergeKeys', true);

T_maestra.Paciente = string(T_maestra.Paciente);
T_maestra.Grupo = string(T_maestra.Grupo);

T_maestra = sortrows(T_maestra, ...
    {'Grupo', 'Paciente'});

%% Colocar Paciente y Grupo al principio

vars = T_maestra.Properties.VariableNames;

vars_inicio = { ...
    'Paciente', ...
    'Grupo'};

vars_resto = setdiff(vars, ...
    vars_inicio, ...
    'stable');

T_maestra = T_maestra(:, ...
    [vars_inicio, vars_resto]);

%% ============================
% COMPROBACIÓN FINAL DE COLUMNAS
%% ============================

fprintf('\nColumnas finales de la tabla maestra:\n');

disp(T_maestra.Properties.VariableNames');

%% ============================
% GUARDAR TABLA MAESTRA
%% ============================

writetable(T_maestra, archivo_tabla_maestra);

fprintf('\nTabla maestra guardada en:\n%s\n', ...
    archivo_tabla_maestra);

%% ============================
% RESUMEN DE GRUPOS
%% ============================

fprintf('\nResumen de pacientes por grupo:\n');

grupos = unique(T_maestra.Grupo);

ResumenGrupos = table();

for i = 1:numel(grupos)

    grupo = grupos(i);

    idx = T_maestra.Grupo == grupo;

    n = sum(idx);

    fprintf('  %s: %d pacientes/filas\n', ...
        grupo, n);

    fila = table();

    fila.Grupo = grupo;
    fila.N = n;

    ResumenGrupos = [ResumenGrupos; fila]; %#ok<AGROW>
end

writetable(ResumenGrupos, archivo_resumen_grupos);

fprintf('\nResumen de grupos guardado en:\n%s\n', ...
    archivo_resumen_grupos);

fprintf('\nFIN: tabla maestra limpia creada.\n');

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function T = leer_y_preparar_tabla(archivo, nombre_bloque)

    if ~isfile(archivo)

        error( ...
            'No se encuentra el archivo del bloque %s:\n%s', ...
            nombre_bloque, ...
            archivo);
    end

    fprintf('  Leyendo bloque %s:\n    %s\n', ...
        nombre_bloque, ...
        archivo);

    T = readtable(archivo);

    T.Properties.VariableNames = ...
        matlab.lang.makeValidName( ...
        T.Properties.VariableNames);

    %% Buscar columna de paciente

    posibles_paciente = { ...
        'Paciente', ...
        'Registro' ...
    };

    col_paciente = encontrar_columna( ...
        T, ...
        posibles_paciente);

    if isempty(col_paciente)

        fprintf('\nColumnas encontradas en %s:\n', archivo);

        disp(T.Properties.VariableNames');

        error([ ...
            'El archivo %s no tiene columna Paciente ', ...
            'ni columna Registro.'], ...
            archivo);
    end

    %% Renombrar Registro como Paciente si hace falta

    if ~strcmp(col_paciente, 'Paciente')

        T.Properties.VariableNames{ ...
            strcmp(T.Properties.VariableNames, ...
            col_paciente)} = ...
            'Paciente';
    end

    %% Comprobar Grupo

    if ~ismember('Grupo', ...
            T.Properties.VariableNames)

        fprintf('\nColumnas encontradas en %s:\n', archivo);

        disp(T.Properties.VariableNames');

        error('El archivo %s no tiene columna Grupo.', ...
            archivo);
    end

    %% Preparar claves

    T.Paciente = strtrim(string(T.Paciente));
    T.Grupo = strtrim(string(T.Grupo));

    %% Eliminar filas sin paciente o grupo

    idx_validas = ...
        strlength(T.Paciente) > 0 & ...
        strlength(T.Grupo) > 0;

    T = T(idx_validas, :);

end

function col = encontrar_columna(T, posibles)

    col = '';

    vars = T.Properties.VariableNames;

    for i = 1:numel(posibles)

        idx = strcmpi(vars, posibles{i});

        if any(idx)

            col = vars{find(idx, 1)};
            return
        end
    end

end

function T = quitar_columnas_si_existen(T, columnas)

    for i = 1:numel(columnas)

        col = columnas{i};

        if ismember(col, ...
                T.Properties.VariableNames)

            T.(col) = [];
        end
    end

end

function T = anadir_prefijo_biomarcadores(T, prefijo)

    vars = T.Properties.VariableNames;

    for i = 1:numel(vars)

        v = vars{i};

        if strcmp(v, 'Paciente') || ...
                strcmp(v, 'Grupo')

            continue
        end

        if ~startsWith(v, prefijo)

            T.Properties.VariableNames{i} = ...
                matlab.lang.makeValidName( ...
                [prefijo, v]);
        end
    end

end