% =========================================================================
% Rilevazione dello Stress (Empatica E4) - Baseline del Paper
% Esame di Interfacce Uomo-Macchina
% Versione il più simile possibile al paper.
% =========================================================================

clear; clc; close all;

%% Configurazione Iniziale
basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
soggetti = dir([basePath 'subject_*']);

% Frequenze originali di campionamento di Empatica E4
fs_eda_orig = 4;
fs_temp_orig = 4;
fs_acc_orig = 32;
fs_bvp = 64;

% Parametri per il windowing (segmentazione)
% Finestre da 60 secondi
dimensioneFinestra = 60;

featureTotali = [];
labelTotali = [];

%% Loop di Elaborazione Dati (Preprocessing e Feature Extraction)

% Estraggo il numero di soggetti e lo scrivo in console
numSoggettiDaElaborare = min(length(soggetti), 29); 
fprintf('Inizio dell''estrazione delle feature per %d soggetti (Baseline Paper)...\n', numSoggettiDaElaborare);

% Ciclo per tutti i soggetti ed elaboro i dati
for s = 1:numSoggettiDaElaborare

    % Apro la cartella del singolo soggetto
    cartellaSoggetto = fullfile(basePath, soggetti(s).name, filesep);
    
    % Prendo i file csv con i dati che ci interessano del soggetto
    file_eda = fullfile(cartellaSoggetto, 'EDA.csv');
    file_bvp = fullfile(cartellaSoggetto, 'BVP.csv');
    file_temp = fullfile(cartellaSoggetto, 'TEMP.csv');
    file_acc = fullfile(cartellaSoggetto, 'ACC.csv');
    
    % Verifico che i file non siano nulli (in più)
    if ~isfile(file_eda) || ~isfile(file_bvp) || ~isfile(file_temp) || ~isfile(file_acc)
        continue; 
    end
    
    % Comunico inizio dell'elaborazione in console del soggetto
    fprintf('Elaborazione del %s...\n', soggetti(s).name);
    
    % Carico i dati come matrice (alla fine i file CSV sono essenzialmente matrici con diverse righe e colonne)
    eda = readmatrix(file_eda);
    bvp = readmatrix(file_bvp);
    temp = readmatrix(file_temp);
    acc = readmatrix(file_acc);
    
    % Per sicurezza, se un file come EDA contenesse dati di inizializzazione come timestamp e frequenze di campionamento, le facciamo saltare.
    if length(eda) > 2 && eda(2) == fs_eda_orig
        eda = eda(3:end); 
    end
    if length(bvp) > 2 && bvp(2) == fs_bvp
        bvp = bvp(3:end); 
    end
    if length(temp) > 2 && temp(2) == fs_temp_orig
        temp = temp(3:end); 
    end
    if size(acc, 1) > 2 && acc(2, 1) == fs_acc_orig
        acc = acc(3:end, :); 
    end
    
    % Rimozione di eventuali NaN
    eda(isnan(eda)) = 0;
    bvp(isnan(bvp)) = 0;
    temp(isnan(temp)) = 0;
    acc(isnan(acc)) = 0;
    
    % Preprocessamento EDA
    % Sovracampionamento da 4Hz a 64Hz per allineare i segnali
    fs_eda = fs_bvp; 
    eda_ricampionato = resample(eda, fs_eda, fs_eda_orig);
    
    % Filtro Passa-Basso (0.5-1 Hz) di secondo ordine
    [b_eda, a_eda] = butter(2, 1/(fs_eda/2), 'low');
    eda_pulito = filtfilt(b_eda, a_eda, eda_ricampionato);
    
    % Decomposizione Tonica/Fasica (0.05 Hz Passa-Basso) con filtro di secondo ordine
    [b_tonico, a_tonico] = butter(2, 0.05/(fs_eda/2), 'low');
    tonico = filtfilt(b_tonico, a_tonico, eda_pulito);
    fasico = eda_pulito - tonico;
    
    % Preprocessamento BVP
    % Filtro Passa-Banda (0.5-5 Hz) di quarto ordine
    [b_bvp, a_bvp] = butter(4, [0.5 5]/(fs_bvp/2), 'bandpass');
    bvp_pulito = filtfilt(b_bvp, a_bvp, bvp);
    
    % Derivate del segnale BVP
    bvp_d1 = diff([bvp_pulito(1); bvp_pulito]) * fs_bvp;
    bvp_d2 = diff([bvp_d1(1); bvp_d1]) * fs_bvp;
    
    % Preprocessamento TEMP e ACC
    % Ricampionamento a fs_bvp (64 Hz)
    temp_ricampionata = resample(temp, fs_bvp, fs_temp_orig);
    
    % ACC magnitudine e ricampionamento
    acc_magnitudo = sqrt(acc(:,1).^2 + acc(:,2).^2 + acc(:,3).^2);
    acc_ricampionato = resample(acc_magnitudo, fs_bvp, fs_acc_orig);

    % Filtro mediano per rimuovere artefatti e spike impulsivi
    temp_pulita = medfilt1(temp_ricampionata, 15); 
    acc_pulito = medfilt1(acc_ricampionato, 15);

    % Ulteriore filtro Passa-Basso sulla TEMP 
    [b_temp, a_temp] = butter(2, 0.1/(fs_bvp/2), 'low');
    temp_pulita = filtfilt(b_temp, a_temp, temp_pulita);
    
    % Segmentazione (Finestre da 60 secondi senza overlay per il baseline)
    durata = floor(length(eda_ricampionato) / fs_eda);

    % Nel paper, il passo è 60 secondi (non c'è overlap)
    passo = 60;
    numFinestre = floor((durata - dimensioneFinestra) / passo) + 1;
    
    % Ciclo per ogni segmento/finestra
    for w = 1:numFinestre
        % Traccio il tempo
        tempoInizioSec = (w-1)*passo;
        indici = tempoInizioSec*fs_bvp + 1 : (tempoInizioSec + dimensioneFinestra)*fs_bvp;
        
        % Verifico che l'indice non superi la lunghezza
        if max(indici) > length(bvp_pulito) || max(indici) > length(temp_pulita) || max(indici) > length(acc_pulito)
            break; 
        end
        
        fin_eda = eda_pulito(indici);
        fin_fasica = fasico(indici);
        fin_bvp = bvp_pulito(indici);
        fin_bvp_d1 = bvp_d1(indici);
        fin_bvp_d2 = bvp_d2(indici);
        fin_temp = temp_pulita(indici);
        fin_acc = acc_pulito(indici);
        
        % Estrazione Feature BVP
        % Rilevamento picchi sistolici (distanza minima 0.4 secondi)
        [~, picchi_bvp] = findpeaks(fin_bvp, 'MinPeakDistance', 0.4*fs_bvp);
        
        % Inizializzo tutto a 0
        f_bvp_media_ppi = 0; f_bvp_devstd_ppi = 0; f_bvp_media_hr = 0; 
        f_bvp_devstd_hr = 0; f_bvp_sd2 = 0; f_bvp_rmssd = 0;
        
        if length(picchi_bvp) > 2
             % Intervalli Picco-A-Picco (secondi)
            ppi = diff(picchi_bvp) / fs_bvp;
            
            % Scarto intervalli anomali: baseline paper usa esattamente il range 0.5-1.2
            ppi_validi = ppi(ppi >= 0.5 & ppi <= 1.2);
            
            % Se ci sono troppi artefatti (>15%), scarto il segmento. Baseline usa soglia del 15% (quindi <0.85 validi)
            if length(ppi_validi) < 0.85 * length(ppi)
                continue;
            end
            
            if length(ppi_validi) > 2
                f_bvp_media_ppi = mean(ppi_validi);
                f_bvp_devstd_ppi = std(ppi_validi);
                f_bvp_media_hr = 60 / f_bvp_media_ppi;
                f_bvp_devstd_hr = std(60 ./ ppi_validi);

                % RMSSD Feature HRV
                f_bvp_rmssd = sqrt(mean(diff(ppi_validi).^2));
                f_bvp_sd2 = std(ppi_validi(1:end-1) + ppi_validi(2:end)) / sqrt(2);
            end
        else
            continue; 
        end
        
        % Feature BVP
        f_bvp_devstd = std(fin_bvp);
        f_bvp_d1_devstd = std(fin_bvp_d1);
        f_bvp_d2_devstd = std(fin_bvp_d2);
        
        % Estrazione Feature EDA
        f_eda_media = mean(fin_eda);
        f_eda_devstd = std(fin_eda);
        
        % Rilevamento picchi EDA sulla fasica
        [picchi_eda, ~] = findpeaks(fin_fasica, 'MinPeakDistance', fs_eda);
        picchi_eda = picchi_eda(picchi_eda > 0.01);
        
        f_eda_picchi = length(picchi_eda);
        if f_eda_picchi > 0
            f_eda_ampiezza = mean(picchi_eda); 
        else
            f_eda_ampiezza = 0;
        end
        
        % Estrazione Feature Temperatura e Accelerometro
        f_temp_media = mean(fin_temp);
        f_temp_devstd = std(fin_temp);
        f_temp_pendenza_fit = polyfit(1:length(fin_temp), fin_temp', 1);
        f_temp_pendenza = f_temp_pendenza_fit(1);
        
        f_acc_media = mean(fin_acc);
        f_acc_devstd = std(fin_acc);
        
        % Vettore di Feature finale
        vettore_feature =[f_eda_media, f_eda_devstd, f_eda_picchi, f_eda_ampiezza, ...
                       f_bvp_devstd, f_bvp_d1_devstd, f_bvp_d2_devstd, ...
                       f_bvp_media_ppi, f_bvp_devstd_ppi, f_bvp_media_hr, f_bvp_devstd_hr, f_bvp_sd2, f_bvp_rmssd, ...
                       f_temp_media, f_temp_devstd, f_temp_pendenza, f_acc_media, f_acc_devstd];
                   
        featureTotali = [featureTotali; vettore_feature];
        
        % Etichettatura (Labeling)
        tempoMedio = (tempoInizioSec + dimensioneFinestra/2) / 60; 
        tempoCorrente = tempoMedio;
        eStress = false;
        
        % Task temporizzati basati sul paper (Rest/riposo di 3 minuti, poi alternati)
        if (tempoCorrente > 3 && tempoCorrente <= 13) || ...  % Task 1 
           (tempoCorrente > 15 && tempoCorrente <= 20) || ... % Task 2 
           (tempoCorrente > 22 && tempoCorrente <= 25) || ... % Task 3 
           (tempoCorrente > 27 && tempoCorrente <= 32) || ... % Task 4 
           (tempoCorrente > 34 && tempoCorrente <= 35)        % Task 5 
            eStress = true;
        end
        
        if eStress
            % 1 = Stress
            labelTotali = [labelTotali; 1];
        else
            % 0 = Baseline (Rest)
            labelTotali = [labelTotali; 0];
        end
    end
end

if isempty(featureTotali)
    error('Nessuna feature estratta! Controllare i percorsi dei file e i dati.');
end

%% Machine Learning e Classificazione Baseline
fprintf('\nAddestramento dei classificatori baseline (10-Fold Cross-Validation)...\n');

% Standardizzazione Z-score
featureTotali = zscore(featureTotali);

% Preparazione CV (10 Fold baseline)
cv = cvpartition(labelTotali, 'KFold', 10);
predizioniRF = zeros(size(labelTotali));
predizioniSVM = zeros(size(labelTotali));

% Random Forest standard e SVM lineare/standard come da paper (senza ottimizzazioni sperimentali come RUSBoost)
for i = 1:cv.NumTestSets
    indiciTrain = cv.training(i);
    indiciTest = cv.test(i);
    
    Xtrain = featureTotali(indiciTrain, :);
    Ytrain = labelTotali(indiciTrain, :);
    Xtest = featureTotali(indiciTest, :);
    
    % Random Forest - Bagging classico, simile al default del paper
    t = templateTree('MaxNumSplits', 100);
    modelloRF = fitcensemble(Xtrain, Ytrain, 'Method', 'Bag', 'NumLearningCycles', 50, 'Learners', t);
    predizioniRF(indiciTest) = predict(modelloRF, Xtest);
    
    % SVM standard - kernel RBF standard senza auto-scale avanzati o cubici manuali
    modelloSVM = fitcsvm(Xtrain, Ytrain, 'KernelFunction', 'rbf');
    predizioniSVM(indiciTest) = predict(modelloSVM, Xtest);
end

%% Valutazione e Metriche di Performance
accRF = sum(predizioniRF == labelTotali) / length(labelTotali);
accSVM = sum(predizioniSVM == labelTotali) / length(labelTotali);

% Metriche per Random Forest
matriceConfRF = confusionmat(labelTotali, predizioniRF);
tn = matriceConfRF(1,1); fp = matriceConfRF(1,2);
fn = matriceConfRF(2,1); tp = matriceConfRF(2,2);

precision = tp / (tp + fp);
recall = tp / (tp + fn);
f1_score = 2 * (precision * recall) / (precision + recall);

fprintf('\n=== RISULTATI MODELLI (BASELINE) ===\n');
fprintf('Distribuzione Classi: %d Rest (0), %d Stress (1)\n', sum(labelTotali==0), sum(labelTotali==1));
fprintf('Accuratezza Random Forest (Bagging): %.2f%%\n', accRF * 100);
fprintf('Accuratezza SVM (RBF standard): %.2f%%\n', accSVM * 100);

fprintf('\n=== METRICHE DETTAGLIATE (Random Forest) ===\n');
fprintf('Precision (Stress): %.2f\n', precision);
fprintf('Recall (Stress): %.2f\n', recall);
fprintf('F1-Score: %.2f\n', f1_score);

% Matrici di Confusione
figure('Name', 'Analisi delle Performance (Baseline Paper)', 'Color', 'w', 'Position',[100 100 900 400]);

subplot(1,2,1);
confusionchart(labelTotali, predizioniRF, 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
title('Matrice Confusione - RF Baseline');

subplot(1,2,2);
confusionchart(labelTotali, predizioniSVM, 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
title('Matrice Confusione - SVM Baseline');

save('processed_data_baseline.mat', 'featureTotali', 'labelTotali');
fprintf('\nDati salvati in "processed_data_baseline.mat".\n');
