% =========================================================================
% Rilevazione dello Stress (Empatica E4) di G.C.
% Esame di Interfacce Uomo-Macchina
% =========================================================================

clear; clc; close all;

%% Configurazione Iniziale
basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
subjects = dir([basePath 'subject_*']);

% Frequenze di campionamento note per Empatica E4
fs_eda_orig = 4;
fs_bvp = 64; 

% Parametri per il windowing (segmentazione)
% Finestre da 60 secondi
windowSize = 60;

allFeatures = [];
allLabels =[];

%% Loop di Elaborazione Dati (Preprocessing e Feature Extraction)
numSubjectsToProcess = min(length(subjects), 29); 
fprintf('Inizio dell''estrazione delle feature per %d soggetti...\n', numSubjectsToProcess);

for s = 1:numSubjectsToProcess
    subjectDir = fullfile(basePath, subjects(s).name, filesep);
    
    eda_file = fullfile(subjectDir, 'EDA.csv');
    bvp_file = fullfile(subjectDir, 'BVP.csv');
    temp_file = fullfile(subjectDir, 'TEMP.csv');
    acc_file = fullfile(subjectDir, 'ACC.csv');
    
    if ~isfile(eda_file) || ~isfile(bvp_file) || ~isfile(temp_file) || ~isfile(acc_file)
        continue; 
    end
    
    fprintf('Elaborazione del %s...\n', subjects(s).name);
    
    % Caricamento Dati
    eda = readmatrix(eda_file);
    bvp = readmatrix(bvp_file);
    temp = readmatrix(temp_file);
    acc = readmatrix(acc_file);
    
    % Controllo di sicurezza (intestazioni)
    if length(eda) > 2 && eda(2) == fs_eda_orig
        eda = eda(3:end); 
    end
    if length(bvp) > 2 && bvp(2) == fs_bvp
        bvp = bvp(3:end); 
    end
    fs_temp_orig = 4;
    if length(temp) > 2 && temp(2) == fs_temp_orig
        temp = temp(3:end); 
    end
    fs_acc_orig = 32;
    if size(acc, 1) > 2 && acc(2, 1) == fs_acc_orig
        acc = acc(3:end, :); 
    end
    
    % Rimozione eventuali NaN
    eda(isnan(eda)) = 0;
    bvp(isnan(bvp)) = 0;
    temp(isnan(temp)) = 0;
    acc(isnan(acc)) = 0;
    
    % Preprocessing EDA
    % Sovracampionamento da 4Hz a 64Hz per allineare i segnali
    fs_eda = fs_bvp; 
    eda_resampled = resample(eda, fs_eda, fs_eda_orig);
    
    % Filtraggio Low-pass (passa-basso) (0.5-1 Hz) per rimuovere rumore ad alta frequenza
    [b_eda, a_eda] = butter(2, 1/(fs_eda/2), 'low');
    eda_clean = filtfilt(b_eda, a_eda, eda_resampled);
    
    % Decomposizione Tonica/Fasica (0.05 Hz Low-pass)
    [b_tonic, a_tonic] = butter(2, 0.05/(fs_eda/2), 'low');
    tonic = filtfilt(b_tonic, a_tonic, eda_clean);
    phasic = eda_clean - tonic;
    
    % Preprocessing BVP
    % Filtro Band-pass (passa-banda) 0.5-5 Hz
    [b_bvp, a_bvp] = butter(4, [0.5 5]/(fs_bvp/2), 'bandpass');
    bvp_clean = filtfilt(b_bvp, a_bvp, bvp);
    
    % Derivate del segnale BVP
    bvp_d1 = diff([bvp_clean(1); bvp_clean]) * fs_bvp;
    bvp_d2 = diff([bvp_d1(1); bvp_d1]) * fs_bvp;
    
    % Preprocessing TEMP e ACC
    % Resample a fs_bvp (64 Hz)
    temp_resampled = resample(temp, fs_bvp, fs_temp_orig);
    
    % ACC magnitude e resampling
    acc_mag = sqrt(acc(:,1).^2 + acc(:,2).^2 + acc(:,3).^2);
    acc_resampled = resample(acc_mag, fs_bvp, fs_acc_orig);

    % Filtro mediano (Non Lineare) per rimuovere artefatti spike impulsivi
    temp_clean = medfilt1(temp_resampled, 15); 
    acc_clean = medfilt1(acc_resampled, 15);

    % Ulteriore filtro Low-Pass sulla TEMP essendo un segnale a lentissima variazione
    [b_temp, a_temp] = butter(2, 0.1/(fs_bvp/2), 'low');
    temp_clean = filtfilt(b_temp, a_temp, temp_clean);
    
    % Segmentazione (Finestre da 60 secondi con overlap di 30)
    duration = floor(length(eda_resampled) / fs_eda);

    % Secondi tra una finestra e l'altra
    stepSize = 30;
    numWindows = floor((duration - windowSize) / stepSize) + 1;
    
    for w = 1:numWindows
        startTimeSec = (w-1)*stepSize;
        idx = startTimeSec*fs_bvp + 1 : (startTimeSec + windowSize)*fs_bvp;
        
        if max(idx) > length(bvp_clean) || max(idx) > length(temp_clean) || max(idx) > length(acc_clean)
            break; 
        end
        
        win_eda = eda_clean(idx);
        win_phasic = phasic(idx);
        win_bvp = bvp_clean(idx);
        win_bvp_d1 = bvp_d1(idx);
        win_bvp_d2 = bvp_d2(idx);
        win_temp = temp_clean(idx);
        win_acc = acc_clean(idx);
        
        % Estrazione Feature BVP
        % Rilevamento picchi sistolici (distanza minima 0.4s)
        [~, bvp_locs] = findpeaks(win_bvp, 'MinPeakDistance', 0.4*fs_bvp);
        
        % Inizializzazione a 0
        f_bvp_mean_ppi = 0; f_bvp_std_ppi = 0; f_bvp_mean_hr = 0; 
        f_bvp_std_hr = 0; f_bvp_sd2 = 0;
        
        if length(bvp_locs) > 2
             % Intervalli Peak-to-Peak (picco-picco) (sec)
            ppi = diff(bvp_locs) / fs_bvp;
            
            % Scarto intervalli anomali (fuori dal range 400ms-1500ms)
            valid_ppi = ppi(ppi >= 0.4 & ppi <= 1.5);
            
            % Se ci sono troppi artefatti (>20% anomalo), scarto il segmento
            if length(valid_ppi) < 0.80 * length(ppi)
                % Ho ancora i traumi Infor in Baan C di questa roba, fantastico tirocinio
                continue;
            end
            
            if length(valid_ppi) > 2
                f_bvp_mean_ppi = mean(valid_ppi);
                f_bvp_std_ppi = std(valid_ppi);
                f_bvp_mean_hr = 60 / f_bvp_mean_ppi;
                f_bvp_std_hr = std(60 ./ valid_ppi);

                % RMSSD (Root Mean Square of Successive Differences) Feature HRV standard
                f_bvp_rmssd = sqrt(mean(diff(valid_ppi).^2));
                f_bvp_sd2 = std(valid_ppi(1:end-1) + valid_ppi(2:end)) / sqrt(2);
            end
        else
            continue; 
        end
        
        % Feature Statistiche BVP
        f_bvp_std = std(win_bvp);
        f_bvp_d1_std = std(win_bvp_d1);
        f_bvp_d2_std = std(win_bvp_d2);
        
        % Estrazione Feature EDA
        f_eda_mean = mean(win_eda);
        f_eda_std = std(win_eda);
        
        % Rilevamento picchi EDA (Skin Conductance Responses SCR) sulla fasica
        [pks_eda, ~] = findpeaks(win_phasic, 'MinPeakDistance', fs_eda);
        pks_eda = pks_eda(pks_eda > 0.01);
        
        f_eda_peaks = length(pks_eda);
        if f_eda_peaks > 0
            f_eda_amp = mean(pks_eda); 
        else
            f_eda_amp = 0;
        end
        
        % Estrazione Feature Temperatura e Accelerometro
        f_temp_mean = mean(win_temp);
        f_temp_std = std(win_temp);
        f_temp_slope_fit = polyfit(1:length(win_temp), win_temp', 1);
        f_temp_slope = f_temp_slope_fit(1);
        
        f_acc_mean = mean(win_acc);
        f_acc_std = std(win_acc);
        
        % Vettore di Feature finale
        feat_vector =[f_eda_mean, f_eda_std, f_eda_peaks, f_eda_amp, ...
                       f_bvp_std, f_bvp_d1_std, f_bvp_d2_std, ...
                       f_bvp_mean_ppi, f_bvp_std_ppi, f_bvp_mean_hr, f_bvp_std_hr, f_bvp_sd2, f_bvp_rmssd, ...
                       f_temp_mean, f_temp_std, f_temp_slope, f_acc_mean, f_acc_std];
                   
        allFeatures = [allFeatures; feat_vector];
        
        % Etichettatura (Labeling)
        % Centro della finestra in minuti
        midTime = (startTimeSec + windowSize/2) / 60; 
        currentTime = midTime;
        isStress = false;
        
        % Task temporizzati basati sul paper (Rest/riposo di 3 minuti, poi Task e Rest alternati)
        if (currentTime > 3 && currentTime <= 13) || ...  % Task 1 
           (currentTime > 15 && currentTime <= 20) || ... % Task 2 
           (currentTime > 22 && currentTime <= 25) || ... % Task 3 
           (currentTime > 27 && currentTime <= 32) || ... % Task 4 
           (currentTime > 34 && currentTime <= 35)        % Task 5 
            isStress = true;
        end
        
        if isStress
            % 1 = Stress
            allLabels = [allLabels; 1];
        else
            % 0 = Baseline (Rest)
            allLabels = [allLabels; 0];
        end
    end
end

if isempty(allFeatures)
    error('Nessuna feature estratta! Controlla i percorsi dei file e i dati.');
end

%% Machine Learning e Classificazione
fprintf('\nAddestramento dei classificatori con 10-Fold Cross-Validation...\n');

% Standardizzazione Z-score (Per comparare le scale)
allFeatures = zscore(allFeatures);

% Verifica presenza di sufficienti esempi di entrambe le classi per il K-Fold
unique_labels = unique(allLabels);
if length(unique_labels) < 2
    error('Errore: I dati contengono solo 1 classe. Impossibile addestrare il classificatore.');
end

% Preparazione CV (10 Fold)
cv = cvpartition(allLabels, 'KFold', 10);
rfPred = zeros(size(allLabels));
svmPred = zeros(size(allLabels));

% Setup Alberi per l'Ensamble Random Forest
t = templateTree('MaxNumSplits', 500, 'MinLeafSize', 5);

for i = 1:cv.NumTestSets
    trainIdx = cv.training(i);
    testIdx = cv.test(i);
    
    Xtrain = allFeatures(trainIdx, :);
    Ytrain = allLabels(trainIdx, :);
    Xtest = allFeatures(testIdx, :);
    
    % Random Forest (RUSBoost per bilanciare il dataset in automatico)
    rfModel = fitcensemble(Xtrain, Ytrain, 'Method', 'RUSBoost', ...
        'NumLearningCycles', 100, 'Learners', t, 'LearnRate', 0.1);
    rfPred(testIdx) = predict(rfModel, Xtest);
    
    % SVM (Kernel RBF Radial Basis Function, migliora i risultati per segnali fisiologici)
    svmModel = fitcsvm(Xtrain, Ytrain, 'KernelFunction', 'rbf', ...
        'BoxConstraint', 10, 'Standardize', false, 'KernelScale', 'auto');
    svmPred(testIdx) = predict(svmModel, Xtest);
end

%% Valutazione e Metriche di Performance
rfAcc = sum(rfPred == allLabels) / length(allLabels);
svmAcc = sum(svmPred == allLabels) / length(allLabels);

% Calcolo Metriche Dettagliate per Random Forest
confMatRF = confusionmat(allLabels, rfPred);
tn = confMatRF(1,1); fp = confMatRF(1,2);
fn = confMatRF(2,1); tp = confMatRF(2,2);

precision = tp / (tp + fp);
recall = tp / (tp + fn);
f1_score = 2 * (precision * recall) / (precision + recall);

fprintf('\n=== RISULTATI MODELLI ===\n');
fprintf('Distribuzione Classi: %d Rest (0), %d Stress (1)\n', sum(allLabels==0), sum(allLabels==1));
fprintf('Accuratezza Random Forest (RUSBoost): %.2f%%\n', rfAcc * 100);
fprintf('Accuratezza SVM (Cubica): %.2f%%\n', svmAcc * 100);

fprintf('\n=== METRICHE DETTAGLIATE (Random Forest) ===\n');
fprintf('Precision (Stress): %.2f\n', precision);
fprintf('Recall (Stress): %.2f\n', recall);
fprintf('F1-Score: %.2f\n', f1_score);

% Matrici di Confusione
figure('Name', 'Analisi delle Performance', 'Color', 'w', 'Position',[100 100 900 400]);

subplot(1,2,1);
confusionchart(allLabels, rfPred, 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
title('Matrice Confusione - RF (RUSBoost)');

subplot(1,2,2);
confusionchart(allLabels, svmPred, 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
title('Matrice Confusione - SVM (Cubica)');

save('processed_data.mat', 'allFeatures', 'allLabels');
fprintf('\nDati salvati in "processed_data.mat".\n');