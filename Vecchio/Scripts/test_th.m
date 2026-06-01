basePath = 'Dataset/EmpaticaE4Stress/EmpaticaE4Stress/Data_29_subjects/Subjects/';
soggetti = dir([basePath 'subject_*']);
fs_bvp = 64; 
% Test thresholds
for th = [0.70, 0.75, 0.80, 0.82, 0.83]
    numR=0; numS=0;
    for s=1:29
        f = fullfile(basePath, soggetti(s).name, 'BVP.csv');
        if ~isfile(f), continue; end
        bvp = readmatrix(f);
        if length(bvp)>2 && bvp(2)==fs_bvp, bvp=bvp(3:end); end
        [b,a] = cheby2(4,20,[0.5 5]/(fs_bvp/2), 'bandpass');
        bvp_pulito = filtfilt(b,a,bvp);
        d_min = floor(length(bvp_pulito)/fs_bvp/60);
        tasks = [0 3 0; 13 15 0; 20 22 0; 25 27 0; max(27,d_min-3) max(27,d_min-1) 0; ...
                 3 13 1; 15 20 1; 22 25 1; 27 max(27,d_min-3) 1; max(27,d_min-1) d_min 1];
        
        for bidx=1:size(tasks,1)
            t_in = tasks(bidx,1); t_end = tasks(bidx,2); lab = tasks(bidx,3);
            if t_end<=t_in, continue; end
            for w=1:floor(t_end-t_in)
                inds = round((t_in*60 + (w-1)*60)*fs_bvp)+1 : round((t_in*60 + w*60)*fs_bvp);
                inds = inds(inds <= length(bvp_pulito));
                if length(inds)<30*fs_bvp, break; end
                seg = bvp_pulito(inds);
                [~,pl] = findpeaks(seg,'MinPeakDistance',0.4*fs_bvp);
                if length(pl)<3, continue; end
                ppi = diff(pl)/fs_bvp;
                val = ppi(ppi>=0.5 & ppi<=1.2);
                if length(val) >= th*length(ppi) && length(val)>=2
                    if lab==0, numR=numR+1; else, numS=numS+1; end
                end
            end
        end
    end
    fprintf('Th=%.2f -> Rest: %d, Stress: %d, Total: %d\n', th, numR, numS, numR+numS);
end
