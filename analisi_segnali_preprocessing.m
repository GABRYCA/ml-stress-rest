% =========================================================================
% Rilevazione dello Stress (Empatica E4) - Baseline del Paper
% Esame di Interfacce Uomo-Macchina
% Grafici - GC
% =========================================================================

clear; clc; close all;

%% Configurazione Iniziale

basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
soggetti = dir([basePath 'subject_*']);

% Frequenze di campionamento secondo Empatica E4
fs_eda_orig = 4;
fs_bvp = 64;

% Parametri per il windowing (segmentazione)
% Finestre da 60 secondi
dimensioneFinestra = 60;

featureTotali = [];
labelTotali = [];

%% Loop di Elaborazione Dati (Preprocessing e Feature Extraction)
numSoggettiDaElaborare = min(length(soggetti), 29); 
fprintf('Inizio dell''estrazione delle 27 feature per %d soggetti...\n', numSoggettiDaElaborare);

for s = 1:numSoggettiDaElaborare
    cartellaSoggetto = fullfile(basePath, soggetti(s).name, filesep);

    file_eda = fullfile(cartellaSoggetto, 'EDA.csv');
    file_bvp = fullfile(cartellaSoggetto, 'BVP.csv');

    if ~isfile(file_eda) || ~isfile(file_bvp)
        continue; 
    end

    fprintf('Elaborazione del %s...\n', soggetti(s).name);

    eda = readmatrix(file_eda);
    bvp = readmatrix(file_bvp);

    eda(isnan(eda)) = 0; bvp(isnan(bvp)) = 0;    

    % Preprocessamento BVP
    % Filtro Chebyshev II di quarto ordine, attenuazione ferma-banda di 20 dB, passa-banda 0.5-5 Hz
    [b_bvp, a_bvp] = cheby2(4, 20, [0.5 5]/(fs_bvp/2), 'bandpass');
    bvp_pulito = filtfilt(b_bvp, a_bvp, bvp);

    % Prima e seconda derivata BVP
    bvp_d1 = diff([bvp_pulito(1); bvp_pulito]) * fs_bvp;
    bvp_d2 = diff([bvp_d1(1); bvp_d1]) * fs_bvp;

    % Preprocessamento EDa
    % Sovracampionamento da 4Hz a 64Hz
    fs_eda = fs_bvp; 
    eda_ricampionato = resample(eda, fs_eda, fs_eda_orig);

    % Filtro gaussiano (40 punti, sigma 400ms = 0.4s)
    sigma_samples = 0.4 * fs_eda; 
    alpha = (40 - 1) / (2 * sigma_samples);
    win_gauss = gausswin(40, alpha);
    win_gauss = win_gauss / sum(win_gauss); 
    eda_pulito = conv(eda_ricampionato, win_gauss, 'same');

    % Estrazione segnale differenziale per le SCR (finestra Bartlett di 20 punti)
    % Ossia estrazione fasica dall'EDA
    eda_diff = diff([eda_pulito(1); eda_pulito]);
    win_bartlett = bartlett(20);
    win_bartlett = win_bartlett / sum(win_bartlett);
    fasico = conv(eda_diff, win_bartlett, 'same');
    
    %% Plot dei Grafici
    % Creazione dell'asse del tempo
    tempo_bvp = (0:length(bvp)-1) / fs_bvp;
    tempo_eda_orig = (0:length(eda)-1) / fs_eda_orig;
    tempo_eda_ric = (0:length(eda_ricampionato)-1) / fs_eda;
    
    % Per una migliore visualizzazione, plottiamo un segmento di 60 secondi
    % (ad es. tra il minuto 5 e il minuto 6 per evitare la fase di adattamento iniziale)
    inizio_plot = 60 * 5; 
    fine_plot = inizio_plot + 60;
    
    idx_bvp = tempo_bvp >= inizio_plot & tempo_bvp <= fine_plot;
    idx_eda_orig = tempo_eda_orig >= inizio_plot & tempo_eda_orig <= fine_plot;
    idx_eda_ric = tempo_eda_ric >= inizio_plot & tempo_eda_ric <= fine_plot;
    
    figure('Name', ['Segnali BVP - ' soggetti(s).name], 'Position', [100 100 1000 600]);
    
    % BVP Originale
    subplot(3,1,1);
    plot(tempo_bvp(idx_bvp), bvp(idx_bvp), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    title('Segnale BVP Originale');
    xlabel('Tempo (s)'); ylabel('Ampiezza');
    grid on;
    
    % BVP Filtrato
    subplot(3,1,2);
    plot(tempo_bvp(idx_bvp), bvp_pulito(idx_bvp), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    title('Segnale BVP Filtrato (Chebyshev II passabanda 0.5-5 Hz)');
    xlabel('Tempo (s)'); ylabel('Ampiezza');
    grid on;
    
    % Derivata BVP
    subplot(3,1,3);
    plot(tempo_bvp(idx_bvp), bvp_d1(idx_bvp), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
    title('Prima Derivata BVP');
    xlabel('Tempo (s)'); ylabel('Ampiezza');
    grid on;
    
    figure('Name', ['Segnali EDA - ' soggetti(s).name], 'Position', [150 150 1000 800]);
    
    % EDA Originale
    subplot(4,1,1);
    plot(tempo_eda_orig(idx_eda_orig), eda(idx_eda_orig), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    title('Segnale EDA Originale (4 Hz)');
    xlabel('Tempo (s)'); ylabel('Ampiezza (\muS)');
    grid on;
    
    % EDA Sovracampionato
    subplot(4,1,2);
    plot(tempo_eda_ric(idx_eda_ric), eda_ricampionato(idx_eda_ric), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
    title('Segnale EDA Sovracampionato a 64 Hz');
    xlabel('Tempo (s)'); ylabel('Ampiezza (\muS)');
    grid on;
    
    % EDA Pulito (Filtrato Gaussiano)
    subplot(4,1,3);
    plot(tempo_eda_ric(idx_eda_ric), eda_pulito(idx_eda_ric), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    title('Segnale EDA Filtrato (Filtro Gaussiano, Rimozione Artefatti)');
    xlabel('Tempo (s)'); ylabel('Ampiezza (\muS)');
    grid on;
    
    % EDA Fasico (SCR)
    subplot(4,1,4);
    plot(tempo_eda_ric(idx_eda_ric), fasico(idx_eda_ric), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
    title('Componente Fasica EDA (Picchi SCR, Filtrata con Finestra Bartlett)');
    xlabel('Tempo (s)'); ylabel('Ampiezza (\muS)');
    grid on;
    
    fprintf('Grafici generati con successo per %s.\n', soggetti(s).name);
    break;
end
