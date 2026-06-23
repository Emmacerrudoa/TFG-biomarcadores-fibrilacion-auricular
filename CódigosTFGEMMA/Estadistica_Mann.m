clear
clc
close all

%% ============================================================
% ANÁLISIS MANN-WHITNEY DE BIOMARCADORES POR GRUPO
%
% Entrada:
%   tabla_maestra_biomarcadores_TFG.xlsx
%
% Comparaciones:
%   1) SANO vs FA_PAROXISTICA_RS
%   2) FA_PAROXISTICA_FA vs FA_PERSISTENTE
%
% Salida:
%   resultados_MannWhitney_biomarcadores_correccion.xlsx
%
% Nota:
%   - En MATLAB, la prueba de Mann-Whitney se realiza con ranksum.
%   - Se ignoran los valores NaN.
%   - Las columnas no numéricas no se analizan.
%   - Las columnas de identificación, control o procesado no se analizan.
%   - En FA activa vs FA persistente no se analizan biomarcadores P/T,
%     porque P y T solo se calcularon en ritmo sinusal.
%   - Los valores de p se ajustan mediante el método Benjamini-Hochberg.
%   - La corrección FDR se aplica de forma independiente en cada comparación.
%% ============================================================

%% ============================
% RUTAS
%% ============================

% Modificar esta ruta según la ubicación local de la tabla maestra.
carpeta_base = 'C:\Users\Emma\Documents\MATLAB\RESULTADOS_ESTADISTICOS_TFG';

archivo_entrada = fullfile(carpeta_base, 'tabla_maestra_biomarcadores_TFG.xlsx');
archivo_salida  = fullfile(carpeta_base, 'resultados_MannWhitney_biomarcadores_correccion.xlsx');

alpha = 0.05;

%% ============================
% LEER TABLA MAESTRA
%% ============================

if ~isfile(archivo_entrada)
    error('No se encuentra la tabla maestra:\n%s', archivo_entrada);
end

T = readtable(archivo_entrada);

T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

if ~ismember('Paciente', T.Properties.VariableNames)
    error('La tabla maestra no tiene columna Paciente.');
end

if ~ismember('Grupo', T.Properties.VariableNames)
    error('La tabla maestra no tiene columna Grupo.');
end

T.Paciente = strtrim(string(T.Paciente));
T.Grupo = strtrim(string(T.Grupo));

fprintf('\nTabla maestra cargada correctamente.\n');
fprintf('Filas: %d\n', height(T));
fprintf('Columnas: %d\n', width(T));

%% ============================
% DEFINIR COMPARACIONES
%% ============================

comparaciones = {
    "SANO",              "FA_PAROXISTICA_RS", "SANO_vs_FA_PAROXISTICA_RS"
    "FA_PAROXISTICA_FA", "FA_PERSISTENTE",    "FA_PAROXISTICA_FA_vs_FA_PERSISTENTE"
};

%% ============================
% IDENTIFICAR BIOMARCADORES
%% ============================

vars = T.Properties.VariableNames;

% Columnas que no son biomarcadores.
% Son columnas de identificación, control de calidad, recuentos o procesado.
% Según la tabla maestra final, estas son las columnas de control reales.
cols_excluir = { ...
    'Paciente', ...
    'Grupo', ...
    'P_CorrGrupoMedia', ...
    'P_CorrGrupoStd', ...
    'T_CorrGrupoMedia', ...
    'T_CorrGrupoStd', ...
    'P_CorrInterGrupoMedia', ...
    'T_CorrInterGrupoMedia', ...
    'P_CorrInterGrupoStd', ...
    'T_CorrInterGrupoStd',   
};

biomarcadores = {};

for i = 1:numel(vars)

    v = vars{i};

    if ismember(v, cols_excluir)
        continue
    end

    if isnumeric(T.(v))
        biomarcadores{end+1} = v; %#ok<SAGROW>
    end
end

fprintf('\nBiomarcadores numéricos detectados: %d\n', numel(biomarcadores));
disp(string(biomarcadores(:)))

%% ============================
% PRUEBA MANN-WHITNEY
%% ============================

Resultados = table;

for c = 1:size(comparaciones, 1)

    grupo1 = comparaciones{c, 1};
    grupo2 = comparaciones{c, 2};
    nombre_comp = comparaciones{c, 3};

    fprintf('\n============================================\n');
    fprintf('Comparacion: %s\n', nombre_comp);
    fprintf('============================================\n');

    idx1 = T.Grupo == grupo1;
    idx2 = T.Grupo == grupo2;

    fprintf('  %s: %d filas\n', grupo1, sum(idx1));
    fprintf('  %s: %d filas\n', grupo2, sum(idx2));

    for b = 1:numel(biomarcadores)

        biom = biomarcadores{b};

        % En FA activa vs FA persistente no se analizan biomarcadores P/T,
        % porque P y T solo se calcularon en ritmo sinusal.
        if nombre_comp == "FA_PAROXISTICA_FA_vs_FA_PERSISTENTE" && ...
           (startsWith(biom, 'P_') || startsWith(biom, 'T_'))
            continue
        end

        x1 = T.(biom)(idx1);
        x2 = T.(biom)(idx2);

        x1 = x1(isfinite(x1));
        x2 = x2(isfinite(x2));

        n1 = numel(x1);
        n2 = numel(x2);

        media1 = NaN;
        media2 = NaN;
        std1 = NaN;
        std2 = NaN;
        mediana1 = NaN;
        mediana2 = NaN;
        iqr1 = NaN;
        iqr2 = NaN;
        p = NaN;
        significativo = "No datos suficientes";
        direccion = "";

        if n1 >= 3 && n2 >= 3

            media1 = mean(x1, 'omitnan');
            media2 = mean(x2, 'omitnan');

            std1 = std(x1, 0, 'omitnan');
            std2 = std(x2, 0, 'omitnan');

            mediana1 = median(x1, 'omitnan');
            mediana2 = median(x2, 'omitnan');

            iqr1 = iqr(x1);
            iqr2 = iqr(x2);

            p = ranksum(x1, x2);

            if p < alpha
                significativo = "Sí";
            else
                significativo = "No";
            end

            if mediana2 > mediana1
                direccion = grupo2 + " > " + grupo1;
            elseif mediana2 < mediana1
                direccion = grupo2 + " < " + grupo1;
            else
                direccion = grupo2 + " = " + grupo1;
            end
        end

        nueva_fila = table;

        nueva_fila.Comparacion = string(nombre_comp);
        nueva_fila.Biomarcador = string(biom);

        nueva_fila.Grupo_1 = string(grupo1);
        nueva_fila.Grupo_2 = string(grupo2);

        nueva_fila.N_Grupo_1 = n1;
        nueva_fila.N_Grupo_2 = n2;

        nueva_fila.Media_Grupo_1 = media1;
        nueva_fila.DE_Grupo_1 = std1;
        nueva_fila.Mediana_Grupo_1 = mediana1;
        nueva_fila.IQR_Grupo_1 = iqr1;

        nueva_fila.Media_Grupo_2 = media2;
        nueva_fila.DE_Grupo_2 = std2;
        nueva_fila.Mediana_Grupo_2 = mediana2;
        nueva_fila.IQR_Grupo_2 = iqr2;

        nueva_fila.p_value = p;
        nueva_fila.Significativo = string(significativo);
        nueva_fila.Direccion_Mediana = string(direccion);

        Resultados = [Resultados; nueva_fila]; %#ok<AGROW>

    end
end

%% ============================
% CORRECCIÓN FDR
% BENJAMINI-HOCHBERG
%% ============================

Resultados.p_FDR = NaN(height(Resultados), 1);
Resultados.Significativo_FDR = strings(height(Resultados), 1);

for c = 1:size(comparaciones, 1)

    nombre_comp = comparaciones{c, 3};

    idx_comp = Resultados.Comparacion == nombre_comp;

    p_originales = Resultados.p_value(idx_comp);

    idx_validos = isfinite(p_originales);
    p_validos = p_originales(idx_validos);

    p_ajustados_comp = NaN(size(p_originales));

    if ~isempty(p_validos)

        m = numel(p_validos);

        [p_ordenados, orden] = sort(p_validos, 'ascend');

        rangos = (1:m)';

        p_ordenados = p_ordenados(:);
        p_bh_ordenados = p_ordenados .* m ./ rangos;

        for i = m-1:-1:1
            p_bh_ordenados(i) = min( ...
                p_bh_ordenados(i), ...
                p_bh_ordenados(i+1));
        end

        p_bh_ordenados = min(p_bh_ordenados, 1);

        p_bh_validos = NaN(m, 1);
        p_bh_validos(orden) = p_bh_ordenados;

        p_ajustados_comp(idx_validos) = p_bh_validos;
    end

    Resultados.p_FDR(idx_comp) = p_ajustados_comp;

end

for i = 1:height(Resultados)

    if ~isfinite(Resultados.p_FDR(i))
        Resultados.Significativo_FDR(i) = "No datos suficientes";
    elseif Resultados.p_FDR(i) < alpha
        Resultados.Significativo_FDR(i) = "Sí";
    else
        Resultados.Significativo_FDR(i) = "No";
    end
end

%% ============================
% ORDENAR RESULTADOS
%% ============================

Resultados = sortrows(Resultados, {'Comparacion', 'p_FDR', 'p_value'});

%% ============================
% AÑADIR COLUMNA PARA UMAP
%% ============================

Resultados.Usar_UMAP = strings(height(Resultados),1);

for i = 1:height(Resultados)

    if Resultados.Significativo_FDR(i) == "Sí"
        Resultados.Usar_UMAP(i) = "Sí";
    else
        Resultados.Usar_UMAP(i) = "No";
    end
end

%% ============================
% REORDENAR COLUMNAS
%% ============================

vars = Resultados.Properties.VariableNames;

vars_inicio = { ...
    'Comparacion', ...
    'Biomarcador', ...
    'Grupo_1', ...
    'Grupo_2', ...
    'N_Grupo_1', ...
    'N_Grupo_2', ...
    'Media_Grupo_1', ...
    'DE_Grupo_1', ...
    'Mediana_Grupo_1', ...
    'IQR_Grupo_1', ...
    'Media_Grupo_2', ...
    'DE_Grupo_2', ...
    'Mediana_Grupo_2', ...
    'IQR_Grupo_2', ...
    'p_value', ...
    'p_FDR', ...
    'Significativo', ...
    'Significativo_FDR', ...
    'Direccion_Mediana', ...
    'Usar_UMAP' ...
};

vars_inicio = vars_inicio(ismember(vars_inicio, vars));
vars_resto = setdiff(vars, vars_inicio, 'stable');

Resultados = Resultados(:, [vars_inicio, vars_resto]);

%% ============================
% GUARDAR EXCEL
%% ============================

writetable(Resultados, archivo_salida);

fprintf('\nResultados Mann-Whitney guardados en:\n%s\n', archivo_salida);

%% ============================
% RESUMEN FINAL
%% ============================

fprintf('\nResumen de biomarcadores significativos:\n');

for c = 1:size(comparaciones, 1)

    nombre_comp = comparaciones{c, 3};

    idx_comp = Resultados.Comparacion == nombre_comp;
    idx_sig = idx_comp & Resultados.Significativo_FDR == "Sí";

    fprintf('\n%s\n', nombre_comp);
    fprintf('  Significativos después de FDR: %d de %d\n', ...
        sum(idx_sig), sum(idx_comp));

    if any(idx_sig)
        disp(Resultados(idx_sig, ...
            {'Biomarcador', 'p_value', 'p_FDR', 'Direccion_Mediana'}));
    end
end

fprintf('\nFIN Mann-Whitney.\n');