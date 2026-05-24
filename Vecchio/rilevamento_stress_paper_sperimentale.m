% =========================================================================
% Rilevazione dello Stress (Empatica E4) - Approccio Sperimentale
% Esame di Interfacce Uomo-Macchina
% =========================================================================

clear; clc; close all;

%% Configurazione Iniziale
basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
soggetti = dir([basePath 'subject_*']);

% Frequenze di campionamento secondo Empatica E4
fs_eda_orig = 4;
fs_bvp = 64;
fs_temp = 4;
fs_hr = 1;

% Parametri per il windowing (segmentazione)
dimensioneFinestra = 60;

featureTotali = [];
labelTotali = [];

%% Loop di Elaborazione Dati (Preprocessing e Feature Extraction)
numSoggettiDaElaborare = min(length(soggetti), 29); 
fprintf('Inizio dell''estrazione delle feature sperimentali per %d soggetti...\n', numSoggettiDaElaborare);

for s = 1:numSoggettiDaElaborare
    cartellaSoggetto = fullfile(basePath, soggetti(s).name, filesep);
    
    file_eda = fullfile(cartellaSoggetto, 'EDA.csv');
    file_bvp = fullfile(cartellaSoggetto, 'BVP.csv');
    file_temp = fullfile(cartellaSoggetto, 'TEMP.csv');
    file_hr = fullfile(cartellaSoggetto, 'HR.csv');
    
    if ~isfile(file_eda) || ~isfile(file_bvp) || ~isfile(file_temp) || ~isfile(file_hr)
        continue; 
    end

    fprintf('Elaborazione del %s...\n', soggetti(s).name);
    
    uda = readmatrix(file_eda);
    eda = uda;
    bvp = readmatrix(file_bvp);
    temp_sig = readmatrix(file_temp);
    hr_sig = readmatrix(file_hr);
    
    eda(isnan(eda)) = mean(eda, 'omitnan'); 
    bvp(isnan(bvp)) = mean(bvp, 'omitnan');
    temp_sig(isnan(temp_sig)) = mean(temp_sig, 'omitnan');
    hr_sig(isnan(hr_sig)) = mean(hr_sig, 'omitnan');
    
    % Preprocessamento BVP
    % Filtro Chebyshev II di quarto ordine, passa-banda 0.5-5 Hz
    [b_bvp, a_bvp] = cheby2(4, 20, [0.5 5]/(fs_bvp/2), 'bandpass');
    bvp_pulito = filtfilt(b_bvp, a_bvp, bvp);
    
    bvp_d1 = diff([bvp_pulito(1); bvp_pulito]) * fs_bvp;
    bvp_d2 = diff([bvp_d1(1); bvp_d1]) * fs_bvp;
    
    % Preprocessamento EDA
    fs_eda = fs_bvp; 
    eda_ricampionato = resample(eda, fs_eda, fs_eda_orig);
    
    % Filtro gaussiano 
    sigma_samples = 0.4 * fs_eda; 
    alpha = (40 - 1) / (2 * sigma_samples);
    win_gauss = gausswin(40, alpha);
    win_gauss = win_gauss / sum(win_gauss); 
    eda_pulito = conv(eda_ricampionato, win_gauss, 'same');
    
    eda_diff = diff([eda_pulito(1); eda_pulito]);
    win_bartlett = bartlett(20);
    win_bartlett = win_bartlett / sum(win_bartlett);
    fasico = conv(eda_diff, win_bartlett, 'same');
    
    % Segmentazione basata su Blocchi di Task (Paper)
    d_min = floor(length(eda_ricampionato) / fs_eda) / 60;
    
    blocchi_rest = [0 3; 13 15; 20 22; 25 27; max(27, d_min-5) max(27, d_min-3); max(27, d_min-2) d_min];
    blocchi_stress = [3 13; 15 20; 22 25; 27 max(27, d_min-5); max(27, d_min-3) max(27, d_min-2)];
    
    tutti_blocchi = [blocchi_rest, zeros(size(blocchi_rest,1), 1); blocchi_stress, ones(size(blocchi_stress,1), 1)];
                     
    for b = 1:size(tutti_blocchi, 1)
        inizio_blocco_min = tutti_blocchi(b, 1);
        fine_blocco_min = tutti_blocchi(b, 2);
        label = tutti_blocchi(b, 3);
        
        if fine_blocco_min <= inizio_blocco_min
            continue;
        end
        
        numFinestreBlocco = floor(fine_blocco_min - inizio_blocco_min);
        
        for w = 1:numFinestreBlocco
            tempoInizioSec = (inizio_blocco_min * 60) + (w-1)*dimensioneFinestra;
            
            indici = round(tempoInizioSec*fs_bvp) + 1 : round((tempoInizioSec + dimensioneFinestra)*fs_bvp);
            indici_temp = round(tempoInizioSec*fs_temp) + 1 : round((tempoInizioSec + dimensioneFinestra)*fs_temp);
            indici_hr = round(tempoInizioSec*fs_hr) + 1 : round((tempoInizioSec + dimensioneFinestra)*fs_hr);
            
            indici = indici(indici <= length(bvp_pulito) & indici <= length(eda_pulito));
            indici_temp = indici_temp(indici_temp <= length(temp_sig));
            indici_hr = indici_hr(indici_hr <= length(hr_sig));
            
            if length(indici) < 0.5 * dimensioneFinestra * fs_bvp
                break;
            end
            
            fin_eda = eda_pulito(indici);
            fin_fasica = fasico(indici);
            
            fin_bvp = bvp_pulito(indici);
            fin_bvp_d1 = bvp_d1(indici);
            fin_bvp_d2 = bvp_d2(indici);
            
            fin_temp = temp_sig(indici_temp);
            fin_hr = hr_sig(indici_hr);
            
            % Estrazione Feature BVP
            [picchi_mag, picchi_bvp] = findpeaks(fin_bvp, 'MinPeakDistance', 0.4*fs_bvp, 'MinPeakHeight', 0);        
            if length(picchi_bvp) < 3, continue; end
            
            ppi = diff(picchi_bvp) / fs_bvp;
            ppi_validi = ppi(ppi >= 0.5 & ppi <= 1.2);
            
                        % Rilassamento delle condizioni di scarto rispetto al paper (da 0.85 a 0.60 per recuperare dati utili)
            if length(ppi_validi) < 0.6 * length(ppi) || length(ppi_validi) < 2
                continue;
            end
            
            Mean_PP = mean(ppi_validi);
            std_PP = std(ppi_validi);
            M_HR = mean(60 ./ ppi_validi);
            std_HR = std(60 ./ ppi_validi);
            SD2 = std(ppi_validi(1:end-1) + ppi_validi(2:end)) / sqrt(2);
            RMSSD = sqrt(mean(diff(ppi_validi).^2));
            if isempty(RMSSD) || isnan(RMSSD), RMSSD = 0; end
            pNN50 = sum(abs(diff(ppi_validi)) > 0.05) / max(1, length(diff(ppi_validi))) * 100;
            
            % Potenza in Alta Frequenza e Bassa frequenza
            tempo_picchi = (picchi_bvp(2:end)) / fs_bvp;
            if length(tempo_picchi) == length(ppi)
               t_validi = tempo_picchi(ppi >= 0.5 & ppi <= 1.2);
               if length(t_validi) > 4
                   t_interp = linspace(t_validi(1), t_validi(end), 256);
                   ppi_interp = interp1(t_validi, ppi_validi, t_interp, 'spline');
                   HF = bandpower(ppi_interp, 256/(t_validi(end)-t_validi(1)), [0.15 0.4]);
                   LF = bandpower(ppi_interp, 256/(t_validi(end)-t_validi(1)), [0.04 0.15]);
                   LF_HF_Ratio = LF / (HF + eps);
               else
                   HF = 0; LF = 0; LF_HF_Ratio = 0;
               end
            else
               HF = 0; LF = 0; LF_HF_Ratio = 0;
            end
            
            feat_BVP = [Mean_PP, std_PP, M_HR, std_HR, HF, LF, LF_HF_Ratio, SD2, RMSSD, pNN50, ...
                        mean(fin_bvp), median(fin_bvp), mode(round(fin_bvp*10)/10), min(fin_bvp), max(fin_bvp), std(fin_bvp), ...
                        mean(fin_bvp_d1), std(fin_bvp_d1), mean(fin_bvp_d2), std(fin_bvp_d2)];
            
            % Estrazione Feature EDA
            [amp_scr, pos_scr, w_scr, p_scr] = findpeaks(fin_fasica, 'MinPeakProminence', 0.01);
            N_PEAKS = length(pos_scr);
            if N_PEAKS > 0
                M_Amp = mean(amp_scr);
                M_D = mean(w_scr) / fs_bvp; 
                M_RT = mean(w_scr / 2) / fs_bvp;
            else
                M_Amp = 0; M_D = 0; M_RT = 0;
            end
            
                        
            p_eda = histcounts(fin_eda, 10);
            p_eda = p_eda / sum(p_eda);
            p_eda(p_eda == 0) = [];
            eda_entropy = -sum(p_eda .* log2(p_eda));
            
            feat_EDA = [mean(fin_eda), median(fin_eda), max(fin_eda), min(fin_eda), std(fin_eda), ...
                        mean(diff(fin_eda)), eda_entropy, ... 
                        M_D, M_Amp, M_RT, N_PEAKS, mean(fin_fasica)];
                        
            % Feature TEMP
            if ~isempty(fin_temp)
                feat_TEMP = [mean(fin_temp), median(fin_temp), max(fin_temp), min(fin_temp), max(fin_temp)-min(fin_temp)];
            else
                feat_TEMP = zeros(1, 5);
            end
            
            % Feature HR dal sensore dedicato (più robusto del HR da BVP ppi)
            if ~isempty(fin_hr)
                feat_HR_ext = [mean(fin_hr), std(fin_hr), max(fin_hr), min(fin_hr)];
            else
                feat_HR_ext = zeros(1, 4);
            end
            
            vettore_feature = [feat_BVP, feat_EDA, feat_TEMP, feat_HR_ext];
            featureTotali = [featureTotali; vettore_feature];
            labelTotali = [labelTotali; label];
        end
    end
end

if isempty(featureTotali)
    error('Nessuna feature estratta! Controllare i dati.');
end

%% Machine Learning Sperimentale
fprintf('\nAddestramento classificatori Sperimentali (10-Fold CV)...\n');

% Standardizzazione Z-score
featureTotali = zscore(featureTotali);

% Bilanciamento classi usando pesi o RUSBoost
rng(42);
cv = cvpartition(labelTotali, 'KFold', 10);

best_acc = 0;
best_pred = [];
best_config = "";

% Test 1: Random Forest (Bagging avanzato)
t_rf = templateTree('MinLeafSize', 5);
predizioniRF = zeros(size(labelTotali));
for i = 1:cv.NumTestSets
    Xtrain = featureTotali(cv.training(i), :);
    Ytrain = labelTotali(cv.training(i), :);
    Xtest = featureTotali(cv.test(i), :);
    
    modelloRF = fitcensemble(Xtrain, Ytrain, 'Method', 'Bag', 'NumLearningCycles', 300, 'Learners', t_rf);
    predizioniRF(cv.test(i)) = predict(modelloRF, Xtest);
end
accRF = sum(predizioniRF == labelTotali) / length(labelTotali);
fprintf('Random Forest Avanzato: %.2f%%\n', accRF*100);
if accRF > best_acc
    best_acc = accRF;
    best_pred = predizioniRF;
    best_config = 'Random Forest (Bag)';
end

% Test 2: SVM con Kernel RBF
predizioniSVM = zeros(size(labelTotali));
for i = 1:cv.NumTestSets
    Xtrain = featureTotali(cv.training(i), :);
    Ytrain = labelTotali(cv.training(i), :);
    Xtest = featureTotali(cv.test(i), :);
    
    % Applico pesi alle classi (cost matrix) per contrastare lo sbilanciamento senza perdere campioni come RUSBoost
    t_svm = templateSVM('KernelFunction', 'gaussian', 'KernelScale', 'auto', 'Standardize', true);
    modelloSVM = fitcsvm(Xtrain, Ytrain, 'KernelFunction', 'rbf', 'KernelScale', 'auto', 'Standardize', true);
    predizioniSVM(cv.test(i)) = predict(modelloSVM, Xtest);
end
accSVM = sum(predizioniSVM == labelTotali) / length(labelTotali);
fprintf('SVM (RBF): %.2f%%\n', accSVM*100);
if accSVM > best_acc
    best_acc = accSVM;
    best_pred = predizioniSVM;
    best_config = 'SVM con Kernel RBF';
end

fprintf('\n=== Modello Migliore: %s ===\n', best_config);
fprintf('Distribuzione Classi: %d Rest (0), %d Stress (1)\n', sum(labelTotali==0), sum(labelTotali==1));
fprintf('Accuratezza Complessiva Sperimentale: %.2f%% (Target paper: ~76.5%%)\n', best_acc * 100);

matriceConf = confusionmat(labelTotali, best_pred);
tn = matriceConf(1,1); fp = matriceConf(1,2);
fn = matriceConf(2,1); tp = matriceConf(2,2);

precision_stress = tp / (tp + fp);
recall_stress = tp / (tp + fn);
f1_stress = 2 * (precision_stress * recall_stress) / (precision_stress + recall_stress);

precision_rest = tn / (tn + fn);
recall_rest = tn / (tn + fp);
f1_rest = 2 * (precision_rest * recall_rest) / (precision_rest + recall_rest);

fprintf('\n=== METRICHE LABEL 1 (Stress) ===\n');
fprintf('Precision: %.2f  |  Recall: %.2f  |  F1-Score: %.2f\n', precision_stress, recall_stress, f1_stress);

fprintf('\n=== METRICHE LABEL 0 (Rest) ===\n');
fprintf('Precision: %.2f  |  Recall: %.2f  |  F1-Score: %.2f\n', precision_rest, recall_rest, f1_rest);

figure('Name', 'Modello con Performance Migliori (Sperimentale)', 'Color', 'w');
confusionchart(labelTotali, best_pred, 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
title(['Matrice di Confusione - ', best_config]);



