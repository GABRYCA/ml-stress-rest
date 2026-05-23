basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
soggetti = dir([basePath 'subject_*']);
fs_bvp = 64; dimensioneFinestra = 60;
numR = 0; numS = 0;
for s=1:length(soggetti)
    file_bvp = fullfile(basePath, soggetti(s).name, 'BVP.csv');
    if ~isfile(file_bvp), continue; end
    bvp = readmatrix(file_bvp);
    if length(bvp) > 2 && bvp(2) == fs_bvp, bvp = bvp(3:end); end
    d = length(bvp)/fs_bvp/60;
    
    tasks_r = [0 3; 13 15; 20 22; 25 27; d-5 d-3; d-2 d];
    tasks_s = [3 13; 15 20; 22 25; 27 d-5; d-3 d-2];
    
    for i=1:size(tasks_r,1)
        dur = tasks_r(i,2) - tasks_r(i,1);
        numR = numR + floor(dur);
    end
    for i=1:size(tasks_s,1)
        dur = tasks_s(i,2) - tasks_s(i,1);
        numS = numS + floor(dur);
    end
end
fprintf('Rest segments: %d\nStress segments: %d\nTotal: %d\n', numR, numS, numR+numS);
