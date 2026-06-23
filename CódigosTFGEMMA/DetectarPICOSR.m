function [locs_R_final, motivo] = DetectarPICOSR(x, locs_R_anot, Fs)

% COMBINAR_R_HIBRIDO_CONSERVADOR
%
% Método híbrido conservador para obtener picos R.
%
% Estrategia:
%   1) Usa las anotaciones como base principal.
%   2) Ajusta las anotaciones al pico real del QRS.
%   3) Limpia posibles falsos R procedentes de las anotaciones.
%   4) Ejecuta detector automático auxiliar.
%   5) Añade detecciones automáticas solo si rellenan huecos largos
%      y generan RR fisiológicos.
%
% Entrada:
%   x          -> señal ECG filtrada
%   locs_R_anot -> picos R anotados en coordenadas locales
%   Fs         -> frecuencia de muestreo
%
% Salida:
%   locs_R_final -> picos R finales
%   motivo       -> mensaje informativo si algo falla

motivo = '';
locs_R_final = [];

x = x(:);
N = length(x);

if isempty(x) || any(~isfinite(x)) || isempty(Fs) || Fs <= 0
    motivo = 'Señal vacía o Fs no válido';
    return
end

%% 1) Limpiar anotaciones

locs_R_anot = limpiar_locs_local(locs_R_anot, N);

%% 2) Ajustar anotaciones y limpiar falsos R anotados

if ~isempty(locs_R_anot)

    locs_R_base = ajustar_R_al_pico(x, locs_R_anot, Fs);
    locs_R_base = limpiar_locs_local(locs_R_base, N);
    locs_R_base = fusionar_locs_cercanos_local(locs_R_base, Fs, 0.12);

    % Las anotaciones .qrs se usan como base, pero no se aceptan directamente.
    % Se eliminan posibles falsos R que generan patrones RR sospechosos
    % y presentan amplitud QRS baja respecto a sus vecinos.
    locs_R_base = limpiar_R_anotados_falsos(x, locs_R_base, Fs);
    locs_R_base = limpiar_locs_local(locs_R_base, N);
    locs_R_base = fusionar_locs_cercanos_local(locs_R_base, Fs, 0.12);

else

    locs_R_base = [];

end

%% 3) Detectar R automáticamente como apoyo

[locs_R_det, ~, motivo_det] = detectar_R(x, Fs);

locs_R_det = limpiar_locs_local(locs_R_det, N);

if ~isempty(locs_R_det)

    locs_R_det = ajustar_R_al_pico(x, locs_R_det, Fs);
    locs_R_det = limpiar_locs_local(locs_R_det, N);
    locs_R_det = fusionar_locs_cercanos_local(locs_R_det, Fs, 0.12);

end

%% 4) Si no hay anotaciones, usar detector

if isempty(locs_R_base)

    if isempty(locs_R_det)
        motivo = ['Sin R anotados y detector vacío: ' motivo_det];
        return
    end

    locs_R_final = locs_R_det;
    locs_R_final = limpiar_locs_local(locs_R_final, N);
    locs_R_final = fusionar_locs_cercanos_local(locs_R_final, Fs, 0.12);

    if numel(locs_R_final) < 3
        motivo = 'Menos de 3 picos R usando solo detector';
    end

    return
end

%% 5) Añadir solo detecciones que rellenen huecos largos

locs_R_final = locs_R_base(:);

tol_cerca = round(0.18 * Fs);

RR_min = 0.30;
RR_max = 1.50;

% Hueco mínimo entre dos R anotados para considerar la inclusión
% de una detección automática intermedia.
RR_gap_min = 0.90;

for i = 1:numel(locs_R_det)

    r_new = locs_R_det(i);

    if isempty(locs_R_final)
        locs_R_final = r_new;
        continue
    end

    dist_min = min(abs(locs_R_final - r_new));

    if dist_min <= tol_cerca
        continue
    end

    locs_ord = sort(locs_R_final(:));

    prev_idx = find(locs_ord < r_new, 1, 'last');
    next_idx = find(locs_ord > r_new, 1, 'first');

    aceptar = false;

    if ~isempty(prev_idx) && ~isempty(next_idx)

        r_prev = locs_ord(prev_idx);
        r_next = locs_ord(next_idx);

        rr_gap = (r_next - r_prev) / Fs;
        rr1 = (r_new - r_prev) / Fs;
        rr2 = (r_next - r_new) / Fs;

        if rr_gap > RR_gap_min && ...
           rr1 >= RR_min && rr1 <= RR_max && ...
           rr2 >= RR_min && rr2 <= RR_max

            aceptar = true;
        end

    elseif ~isempty(prev_idx)

        r_prev = locs_ord(prev_idx);
        rr1 = (r_new - r_prev) / Fs;

        if rr1 >= RR_min && rr1 <= RR_max
            aceptar = true;
        end

    elseif ~isempty(next_idx)

        r_next = locs_ord(next_idx);
        rr2 = (r_next - r_new) / Fs;

        if rr2 >= RR_min && rr2 <= RR_max
            aceptar = true;
        end
    end

    if aceptar
        locs_R_final(end+1,1) = r_new; %#ok<AGROW>
        locs_R_final = sort(locs_R_final(:));
    end
end

%% 6) Limpieza final y ajuste final

locs_R_final = limpiar_locs_local(locs_R_final, N);
locs_R_final = fusionar_locs_cercanos_local(locs_R_final, Fs, 0.12);

locs_R_final = ajustar_R_al_pico(x, locs_R_final, Fs);
locs_R_final = limpiar_locs_local(locs_R_final, N);
locs_R_final = fusionar_locs_cercanos_local(locs_R_final, Fs, 0.12);

% Limpieza final por si tras el ajuste queda algún falso R
locs_R_final = limpiar_R_anotados_falsos(x, locs_R_final, Fs);
locs_R_final = limpiar_locs_local(locs_R_final, N);
locs_R_final = fusionar_locs_cercanos_local(locs_R_final, Fs, 0.12);

if numel(locs_R_final) < 3
    motivo = 'Menos de 3 picos R tras combinación híbrida';
    return
end

%% 7) Control final informativo de RR

RR_final = diff(locs_R_final) / Fs;
RR_final = RR_final(isfinite(RR_final));

if ~isempty(RR_final)

    porc_fuera = 100 * sum(RR_final < RR_min | RR_final > RR_max) / numel(RR_final);

    if porc_fuera > 20
        motivo = sprintf('Aviso: %.1f%% de RR fuera de 0.30-1.50 s tras combinación híbrida', porc_fuera);
    end
end

end

function locs_out = limpiar_R_anotados_falsos(x, locs_R, Fs)

% LIMPIAR_R_ANOTADOS_FALSOS
%
% Elimina falsos R procedentes de anotaciones o lista final.
%
% Criterios:
%   1) Si hay dos R demasiado cercanos, conserva el de mayor amplitud.
%   2) Si un R central genera patrón corto-largo/largo-corto y es débil,
%      se elimina.
%   3) Si un R tiene amplitud QRS muy baja respecto a la mediana global,
%      se elimina.

x = x(:);
N = length(x);

locs_R = limpiar_locs_local(locs_R, N);
locs_R = sort(locs_R(:));

if numel(locs_R) < 5
    locs_out = locs_R;
    return
end

RR_min = 0.30;

% Umbral relativo para descartar picos con amplitud QRS baja
% respecto a la amplitud global de la ventana.
factor_amp_global = 0.40;

cambios = true;

while cambios

    cambios = false;

    if numel(locs_R) < 5
        break
    end

    %% Calcular amplitudes locales de QRS

    amp_R = nan(numel(locs_R),1);

    for i = 1:numel(locs_R)
        amp_R(i) = amplitud_QRS_local(x, locs_R(i), Fs);
    end

    amp_valid = amp_R(isfinite(amp_R));

    if numel(amp_valid) < 5
        break
    end

    amp_ref_global = median(amp_valid, 'omitnan');

    if ~isfinite(amp_ref_global) || amp_ref_global <= 0
        break
    end

    %% 1) Eliminar pares demasiado cercanos

    RR = diff(locs_R) / Fs;

    idx_corto = find(RR < RR_min, 1, 'first');

    if ~isempty(idx_corto)

        i1 = idx_corto;
        i2 = idx_corto + 1;

        amp1 = amp_R(i1);
        amp2 = amp_R(i2);

        if ~isfinite(amp1)
            amp1 = 0;
        end

        if ~isfinite(amp2)
            amp2 = 0;
        end

        if amp1 >= amp2
            locs_R(i2) = [];
        else
            locs_R(i1) = [];
        end

        cambios = true;
        continue
    end

    %% 2) Eliminar R central débil con patrón corto-largo o largo-corto

    eliminar_idx = [];

    for i = 2:numel(locs_R)-1

        rr_prev = (locs_R(i) - locs_R(i-1)) / Fs;
        rr_next = (locs_R(i+1) - locs_R(i)) / Fs;

        amp_prev = amp_R(i-1);
        amp_i    = amp_R(i);
        amp_next = amp_R(i+1);

        amp_vecinos = median([amp_prev amp_next], 'omitnan');

        if ~isfinite(amp_i) || ~isfinite(amp_vecinos) || amp_vecinos <= 0
            continue
        end

        patron_corto_largo = rr_prev < 0.55 && rr_next > 0.90;
        patron_largo_corto = rr_prev > 0.90 && rr_next < 0.55;

        R_central_debil = amp_i < 0.60 * amp_vecinos;

        if (patron_corto_largo || patron_largo_corto) && R_central_debil
            eliminar_idx = i;
            break
        end
    end

    if ~isempty(eliminar_idx)
        locs_R(eliminar_idx) = [];
        cambios = true;
        continue
    end

    %% 3) Eliminar R con amplitud global demasiado baja

    idx_bajos = find(amp_R < factor_amp_global * amp_ref_global);

    if ~isempty(idx_bajos)

        % Eliminar el más pequeño primero
        [~, jmin] = min(amp_R(idx_bajos));
        eliminar_idx = idx_bajos(jmin);

        locs_R(eliminar_idx) = [];
        cambios = true;
        continue
    end

end

locs_out = limpiar_locs_local(locs_R, N);

end

function amp = amplitud_QRS_local(x, R, Fs)

x = x(:);
N = length(x);

amp = NaN;

if isempty(x) || ~isfinite(R) || R < 1 || R > N
    return
end

vent_qrs  = round(0.04 * Fs);  % 40 ms alrededor del R
vent_base = round(0.12 * Fs);  % 120 ms para estimar línea local

ini_qrs = max(1, R - vent_qrs);
fin_qrs = min(N, R + vent_qrs);

ini_base = max(1, R - vent_base);
fin_base = min(N, R + vent_base);

seg_qrs = x(ini_qrs:fin_qrs);
seg_base = x(ini_base:fin_base);

if isempty(seg_qrs) || isempty(seg_base) || ...
   any(~isfinite(seg_qrs)) || any(~isfinite(seg_base))
    return
end

baseline = median(seg_base, 'omitnan');

amp = max(abs(seg_qrs - baseline), [], 'omitnan');

end

function [locs_R, pks_R, motivo] = detectar_R(x, Fs)

x = x(:);
motivo = '';

locs_R = [];
pks_R = [];

if isempty(x) || all(~isfinite(x))
    motivo = 'Señal vacía o no finita';
    return
end

x0 = x - median(x, 'omitnan');

dx = diff(x0);
dx(end+1) = dx(end);

sq = dx.^2;

win = round(0.15 * Fs);
win = max(win,1);

x_int = movmean(sq, win);

if isempty(x_int) || all(~isfinite(x_int))
    motivo = 'Señal integrada vacía o no finita';
    return
end

max_int = max(x_int);

if ~isfinite(max_int) || max_int <= 0
    motivo = 'Máximo de señal integrada no válido';
    return
end

umbral = 0.35 * max_int;

try
    [~, locs] = findpeaks(x_int, ...
        'MinPeakHeight', umbral, ...
        'MinPeakDistance', round(0.30 * Fs));
catch
    locs = [];
end

if isempty(locs)
    motivo = 'findpeaks no detecta picos';
    return
end

locs_R = zeros(size(locs));
pks_R = zeros(size(locs));

ventana_busqueda = round(0.06 * Fs);

for k = 1:length(locs)

    ini = max(1, locs(k) - ventana_busqueda);
    fin = min(length(x0), locs(k) + ventana_busqueda);

    segmento = x0(ini:fin);

    if isempty(segmento) || any(~isfinite(segmento))
        continue
    end

    segmento_c = segmento - median(segmento, 'omitnan');

    [~, idx_peak] = max(abs(segmento_c));

    locs_R(k) = ini + idx_peak - 1;
    pks_R(k) = x0(locs_R(k));
end

idx_valid = locs_R > 0 & isfinite(locs_R);

locs_R = locs_R(idx_valid);
pks_R = pks_R(idx_valid);

[locs_R, ia] = unique(round(locs_R));
pks_R = pks_R(ia);

%% Limpieza por amplitud absoluta

if ~isempty(pks_R)

    amp_abs = abs(pks_R);
    amp_ref = prctile(amp_abs, 75);

    if isfinite(amp_ref) && amp_ref > 0

        idx_ok = amp_abs >= 0.25 * amp_ref;

        locs_R = locs_R(idx_ok);
        pks_R = pks_R(idx_ok);
    end
end

if numel(locs_R) < 3
    motivo = 'Menos de 3 picos R tras limpieza';
end

end

function locs = limpiar_locs_local(locs, N)

if isempty(locs)
    locs = [];
    return
end

locs = round(locs(:));
locs = locs(isfinite(locs));
locs = unique(locs);
locs = locs(locs >= 1 & locs <= N);

end

%% ============================================================
% FUNCIÓN LOCAL: AJUSTAR PICOS R AL PICO REAL DEL QRS
% Versión adaptativa según polaridad dominante del QRS
% ============================================================

function locs_R_aj = ajustar_R_al_pico(x, locs_R, Fs)

x = x(:);
locs_R = round(locs_R(:));

locs_R = unique(locs_R);
locs_R = locs_R(locs_R >= 1 & locs_R <= length(x));

if isempty(locs_R)
    locs_R_aj = [];
    return
end

%% Ventana de búsqueda alrededor del R original

busqueda_ms = 80;
busqueda_muestras = round(busqueda_ms/1000 * Fs);

%% ============================================================
% 1) ESTIMAR POLARIDAD DOMINANTE DEL QRS
% ============================================================

amps_qrs = nan(numel(locs_R), 1);

for i = 1:numel(locs_R)

    r0 = locs_R(i);

    ini = max(1, r0 - busqueda_muestras);
    fin = min(length(x), r0 + busqueda_muestras);

    segmento = x(ini:fin);

    if isempty(segmento) || all(~isfinite(segmento))
        continue
    end

    % Buscar el punto de mayor amplitud absoluta solo para estimar polaridad
    [~, idx_abs] = max(abs(segmento));

    amps_qrs(i) = segmento(idx_abs);

end

polaridad = median(amps_qrs, 'omitnan');

if ~isfinite(polaridad)
    locs_R_aj = [];
    return
end

%% ============================================================
% 2) AJUSTAR CADA R SEGÚN LA POLARIDAD DOMINANTE
% ============================================================

locs_R_aj = nan(size(locs_R));

for i = 1:numel(locs_R)

    r0 = locs_R(i);

    ini = max(1, r0 - busqueda_muestras);
    fin = min(length(x), r0 + busqueda_muestras);

    segmento = x(ini:fin);

    if isempty(segmento) || all(~isfinite(segmento))
        continue
    end

    if polaridad >= 0
        % QRS predominantemente positivo: ajustar al máximo positivo
        [~, idx_local] = max(segmento);
    else
        % QRS predominantemente negativo: ajustar al mínimo negativo
        [~, idx_local] = min(segmento);
    end

    locs_R_aj(i) = ini + idx_local - 1;

end

%% ============================================================
% 3) LIMPIEZA FINAL
% ============================================================

locs_R_aj = locs_R_aj(isfinite(locs_R_aj));
locs_R_aj = unique(round(locs_R_aj));
locs_R_aj = locs_R_aj(locs_R_aj >= 1 & locs_R_aj <= length(x));

end

function locs_out = fusionar_locs_cercanos_local(locs, Fs, tol_s)

locs = round(locs(:));
locs = locs(isfinite(locs));
locs = unique(locs);

locs_out = [];

if isempty(locs)
    return
end

tol = round(tol_s * Fs);

i = 1;

while i <= numel(locs)

    grupo = locs(i);
    j = i + 1;

    while j <= numel(locs) && abs(locs(j) - grupo(end)) <= tol
        grupo(end+1,1) = locs(j); %#ok<AGROW>
        j = j + 1;
    end

    centro = median(grupo);
    [~, idx_best] = min(abs(grupo - centro));

    locs_out(end+1,1) = grupo(idx_best); %#ok<AGROW>

    i = j;
end

locs_out = round(locs_out(:));
locs_out = unique(locs_out);

end
