% Pulisco la console
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
fs_target = fs_bvp;

% Caricamento
eda_orig = readmatrix(file_eda);
bvp_orig = readmatrix(file_bvp);
temp_orig = readmatrix(file_temp);
acc_orig = readmatrix(file_acc);

% Aggiustamenti soliti
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

% Rimozione NaN
eda_orig(isnan(eda_orig)) = 0;
bvp_orig(isnan(bvp_orig)) = 0;
temp_orig(isnan(temp_orig)) = 0;
acc_orig(isnan(acc_orig)) = 0;

%% EDA
eda_ricampionato = resample(eda_orig, fs_target, fs_eda_orig);
[b_eda, a_eda] = butter(2, 1/(fs_target/2), 'low');
eda_pulito = filtfilt(b_eda, a_eda, eda_ricampionato);

[b_tonico, a_tonico] = butter(2, 0.05/(fs_target/2), 'low');
eda_tonico = filtfilt(b_tonico, a_tonico, eda_pulito);
eda_fasico = eda_pulito - eda_tonico;

%% BVP
[b_bvp, a_bvp] = butter(4, [0.5 5]/(fs_target/2), 'bandpass');
bvp_pulito = filtfilt(b_bvp, a_bvp, bvp_orig);

% Derivate
bvp_d1 = diff([bvp_pulito(1); bvp_pulito]) * fs_target;
bvp_d2 = diff([bvp_d1(1); bvp_d1]) * fs_target;

%% TEMP
temp_ricampionata = resample(temp_orig, fs_target, fs_temp_orig);
temp_med = medfilt1(temp_ricampionata, 15);
[b_temp, a_temp] = butter(2, 0.1/(fs_target/2), 'low');
temp_pulita = filtfilt(b_temp, a_temp, temp_med);

%% ACC
acc_magnitudo = sqrt(acc_orig(:,1).^2 + acc_orig(:,2).^2 + acc_orig(:,3).^2);
acc_ricampionato = resample(acc_magnitudo, fs_target, fs_acc_orig);
acc_pulito = medfilt1(acc_ricampionato, 15);

%% Vettori Temporali per la visualizzazione
t_eda_orig = (0:length(eda_orig)-1) / fs_eda_orig;
t_eda_ricampionato = (0:length(eda_ricampionato)-1) / fs_target;

t_bvp_orig = (0:length(bvp_orig)-1) / fs_bvp;

t_temp_orig = (0:length(temp_orig)-1) / fs_temp_orig;
t_temp_ricampionato = (0:length(temp_ricampionata)-1) / fs_target;

t_acc_orig = (0:length(acc_orig)-1) / fs_acc_orig;
t_acc_ricampionato = (0:length(acc_ricampionato)-1) / fs_target;

%% Visualizzazione 
% Limiti per mostrare meglio i dettagli dei filtri
x_lim_bvp = [300, 320];

figure('Name', 'Effetti del Preprocessing sui Segnali (Prima e Dopo)', 'NumberTitle', 'off', 'Position', [100 100 1200 800]);

% EDA Plot
subplot(4, 2, 1);
plot(t_eda_orig, eda_orig, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_eda_ricampionato, eda_pulito, 'b', 'LineWidth', 1.5);
title('EDA: Grezzo vs Elaborato (Passa-Basso)');
legend('Grezzo (4Hz)', 'Pulito (64Hz)');
xlabel('Tempo (s)'); ylabel('\muS');
xlim([0 t_eda_orig(end)]);

subplot(4, 2, 2);
plot(t_eda_ricampionato, eda_tonico, 'r', 'LineWidth', 1.5); hold on;
plot(t_eda_ricampionato, eda_fasico, 'g', 'LineWidth', 1.2);
title('EDA: Componente Tonica e Fasica');
legend('Tonica', 'Fasica');
xlabel('Tempo (s)'); ylabel('\muS');
xlim([0 t_eda_orig(end)]);

% BVP Plot
subplot(4, 2, 3);
plot(t_bvp_orig, bvp_orig, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_bvp_orig, bvp_pulito, 'b', 'LineWidth', 1.2);
title('BVP: Grezzo vs Elaborato (Passa-Banda)');
legend('Grezzo', 'Filtrato');
xlabel('Tempo (s)'); ylabel('Ampiezza');
xlim(x_lim_bvp);

subplot(4, 2, 4);
plot(t_bvp_orig, bvp_pulito, 'b'); hold on;
plot(t_bvp_orig, bvp_d1/max(abs(bvp_d1))*max(bvp_pulito), 'r--');
title('BVP: Segnale Pulito e sua Derivata (Normalizzata)');
legend('BVP Pulito', 'Derivata (1a)');
xlabel('Tempo (s)'); ylabel('Ampiezza');
xlim(x_lim_bvp);

% TEMP Plot
subplot(4, 2, 5);
plot(t_temp_orig, temp_orig, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_temp_ricampionato, temp_pulita, 'r', 'LineWidth', 1.5);
title('TEMP: Grezzo vs Elaborato (Mediano + Passa-Basso)');
legend('Grezzo (4Hz)', 'Pulito (64Hz)');
xlabel('Tempo (s)'); ylabel('°C');
xlim([0 t_temp_orig(end)]);

subplot(4, 2, 6);
plot(t_temp_orig, temp_orig, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_temp_ricampionato, temp_pulita, 'r', 'LineWidth', 1.5);
title('TEMP: Dettaglio dell''effetto filtro (Zoom su variazioni)');
xlabel('Tempo (s)'); ylabel('°C');
xlim([400, 1000]);

% ACC Plot
subplot(4, 2, 7);
plot(t_acc_orig, acc_magnitudo, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_acc_ricampionato, acc_pulito, 'k', 'LineWidth', 1.2);
title('ACC: Magnitudo Grezza vs Filtrata (Filtro Mediano)');
legend('Grezzo (32Hz)', 'Pulito (64Hz)');
xlabel('Tempo (s)'); ylabel('g');
xlim([0 t_acc_orig(end)]);

subplot(4, 2, 8);
plot(t_acc_orig, acc_magnitudo, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_acc_ricampionato, acc_pulito, 'k', 'LineWidth', 1.2);
title('ACC: Dettaglio Artefatti rimossi (Zoom)');
xlabel('Tempo (s)'); ylabel('g');
xlim([300, 350]);

sgtitle('Fase di Preprocessing del Segnale: Da Grezzo alle Feature ML');
