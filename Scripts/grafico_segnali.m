% Grafico segnali prima e dopo i filtri
clear; clc; close all;

% Percorso file
basePath = '../Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
if ~exist(basePath, 'dir')
    basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
end

soggetto_dir = 'subject_01';
cartellaSoggetto = fullfile(basePath, soggetto_dir, filesep);

file_eda = fullfile(cartellaSoggetto, 'EDA.csv');
file_bvp = fullfile(cartellaSoggetto, 'BVP.csv');
file_temp = fullfile(cartellaSoggetto, 'TEMP.csv');
file_acc = fullfile(cartellaSoggetto, 'ACC.csv');

if ~isfile(file_eda) || ~isfile(file_bvp) || ~isfile(file_temp) || ~isfile(file_acc)
    error('File non trovati. Controlla il path del dataset.');
end

% Frequenze di campionamento originali
fs_eda_orig = 4;
fs_temp_orig = 4;
fs_acc_orig = 32;
fs_bvp = 64;

% Caricamento
eda_orig = readmatrix(file_eda);
bvp_orig = readmatrix(file_bvp);
temp_orig = readmatrix(file_temp);
acc_orig = readmatrix(file_acc);

% Pulizia dati iniziali
if length(eda_orig) > 2 && eda_orig(2) == fs_eda_orig
    eda_orig = eda_orig(3:end); 
end
if length(bvp_orig) > 2 && bvp_orig(2) == fs_bvp
    bvp_orig = bvp_orig(3:end); 
end
if length(temp_orig) > 2 && temp_orig(2) == fs_temp_orig
    temp_orig = temp_orig(3:end); 
end
if size(acc_orig, 1) > 2 && acc_orig(2, 1) == fs_acc_orig
    acc_orig = acc_orig(3:end, :); 
end

% Rimozione eventuali NaN
eda_orig(isnan(eda_orig)) = 0;
bvp_orig(isnan(bvp_orig)) = 0;
temp_orig(isnan(temp_orig)) = 0;
acc_orig(isnan(acc_orig)) = 0;

fs_target = fs_bvp;

%% EDA
eda_ricampionato = resample(eda_orig, fs_target, fs_eda_orig);

% Filtro Low-pass (1 Hz)
[b_eda, a_eda] = butter(2, 1/(fs_target/2), 'low');
eda_filtrato = filtfilt(b_eda, a_eda, eda_ricampionato);

% Decomposizione Tonica/Fasica (0.05 Hz)
[b_tonico, a_tonico] = butter(2, 0.05/(fs_target/2), 'low');
eda_tonico = filtfilt(b_tonico, a_tonico, eda_filtrato);
eda_fasico = eda_filtrato - eda_tonico;

%% BVP
% Filtro Passa-Banda (0.5-5 Hz)
[b_bvp, a_bvp] = butter(4, [0.5 5]/(fs_bvp/2), 'bandpass');
bvp_filtrato = filtfilt(b_bvp, a_bvp, bvp_orig);

%% TEMP
temp_ricampionato = resample(temp_orig, fs_target, fs_temp_orig);
temp_med = medfilt1(temp_ricampionato, 15);

% Filtro Low-pass (0.1 Hz)
[b_temp, a_temp] = butter(2, 0.1/(fs_target/2), 'low');
temp_filtrato = filtfilt(b_temp, a_temp, temp_med);

%% Processing ACC
acc_mag = sqrt(acc_orig(:,1).^2 + acc_orig(:,2).^2 + acc_orig(:,3).^2);
acc_ricampionato = resample(acc_mag, fs_target, fs_acc_orig);
acc_filtrato = medfilt1(acc_ricampionato, 15);

%% Visualizzazione
% Definizione dei vettori temporali per il plotting
t_eda_orig = (0:length(eda_orig)-1) / fs_eda_orig;
t_eda_ricampionato = (0:length(eda_ricampionato)-1) / fs_target;

t_bvp_orig = (0:length(bvp_orig)-1) / fs_bvp;

t_temp_orig = (0:length(temp_orig)-1) / fs_temp_orig;
t_temp_ricampionato = (0:length(temp_ricampionato)-1) / fs_target;

t_acc_orig = (0:length(acc_mag)-1) / fs_acc_orig;
t_acc_ricampionato = (0:length(acc_ricampionato)-1) / fs_target;

% Plot EDA
figure('Name', 'Segnale EDA prima e dopo il Filtraggio', 'NumberTitle', 'off');
subplot(3, 1, 1);
plot(t_eda_orig, eda_orig);
title('EDA Originale');
xlabel('Tempo (s)'); ylabel('Ampiezza');

subplot(3, 1, 2);
plot(t_eda_ricampionato, eda_ricampionato);
title('EDA Ricampionato (64Hz)');
xlabel('Tempo (s)'); ylabel('Ampiezza');

subplot(3, 1, 3);
plot(t_eda_ricampionato, eda_filtrato, 'b'); hold on;
plot(t_eda_ricampionato, eda_tonico, 'r', 'LineWidth', 1.5);
plot(t_eda_ricampionato, eda_fasico, 'g');
legend('EDA Filtrato', 'Tonico', 'Fasico');
title('EDA Filtrato e Componenti Tonico/Fasico');
xlabel('Tempo (s)'); ylabel('Ampiezza');

% Plot BVP
figure('Name', 'Segnale BVP prima e dopo il Filtraggio', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_bvp_orig, bvp_orig);
title('BVP Originale');
xlabel('Tempo (s)'); ylabel('Ampiezza');

subplot(2, 1, 2);
plot(t_bvp_orig, bvp_filtrato);
title('BVP Filtrato Passa-Banda');
xlabel('Tempo (s)'); ylabel('Ampiezza');

% Plot TEMP
figure('Name', 'Segnale TEMP prima e dopo il Filtraggio', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_temp_orig, temp_orig);
title('TEMP Originale');
xlabel('Tempo (s)'); ylabel('Ampiezza');

subplot(2, 1, 2);
plot(t_temp_ricampionato, temp_filtrato);
title('TEMP Ricampionato e Filtrato');
xlabel('Tempo (s)'); ylabel('Ampiezza');

% Plot ACC
figure('Name', 'Segnale ACC prima e dopo il filtraggio', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_acc_orig, acc_mag);
title('ACC Originale (Magnitude)');
xlabel('Tempo (s)'); ylabel('Ampiezza');

subplot(2, 1, 2);
plot(t_acc_ricampionato, acc_filtrato);
title('ACC Ricampionato e con Filtro Mediano (Magnitude)');
xlabel('Tempo (s)'); ylabel('Ampiezza');
