## Esame di Interfacce Uomo-Macchina

Paper scelto: PPG and EDA dataset collected with Empatica E4 for stress assessment (dataset_ppg_eda_stress_datainbrief2024.pdf).

Paper di supporto:
A Method for Stress Detection Using Empatica E4 Bracelet and
Machine-Learning Techniques (Empatica_dataset_stress_Sensors2024.pdf).

Link al dataset: https://data.mendeley.com/datasets/kb42z77m2g/2 (Scaricato nella cartella Dataset)
Il Dataset non dovrebbe finire su git come da .gitignore (si spera).

Esecuzione: 
Matlab, con:
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Wavelet Toolbox
- Optimization Toolbox
- Deep Learning Toolbox (Forse sono un po' ambizioso, GC)

## Installazione ed Esecuzione

### SETUP:
- Scaricare il Dataset da https://data.mendeley.com/datasets/kb42z77m2g/2
- Creare una cartella "Dataset" ed Estrarlo al suo interno, alla fine la struttura dovrebbe essere la seguente (si, doppia cartella empatica una dentro l'altra):
(Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/...).

Come IDE il progetto prevede che sia installato Matlab con i pacchetti menzionati all'inizio in Esecuzione, ossia:
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Wavelet Toolbox
- Optimization Toolbox

Inoltre, non fa male VSCode, con estensione Matlab installata (non Matlab Unofficial, quella è deprecata), così da poter compilare con il triangolo in alto a destra.
In alternativa, direttamente dalla cartella del progetto, si può eseguire il codice e creare il modello con il comando:
- matlab -batch "main_stress_detection"

Dovrebbe creare il modello, generare in automatico una matrice di confusione e dei log in console, oltre ad un file processed_data.mat.

## Presentazione Powerpoint
C'è una versione Work-In-Progress della presentazione direttamente nella root del progetto, si chiama:
- Progetto Interfacce Uomo-Macchina.pptx.

Vi chiedo di fare attenzione a quello che modificate, salvate sempre e fate molte commit (così i rollback sono più semplici, in caso di errori).



