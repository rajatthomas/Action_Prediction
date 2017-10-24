function data_comb = get_data(run_num, movie_type)
% helper function for hierarchical_intactness.m

if ~ischar(run_num)
    run_num = num2str(run_num);
end

Amovies = [23, 4, 31, 10, 20, 36, 25, 8];
Bmovies = [11, 28, 21, 1, 5, 40, 9, 27];
Cmovies = [7, 30, 12, 24, 29, 16, 32, 3];

Dmovies = [1, 5, 8, 9, 10, 11, 12, 16, 20]; % some are time warped others super intact
Emovies = [8, 11, 20, 5, 1, 9, 16, 12, 10];

if strcmp(movie_type, 'intact')
    Dmovie_select = Dmovies;
end

if strcmp(movie_type, 'scrambled')
    Dmovie_select = Dmovies+20;
end

if strcmp(movie_type, 'intact') || strcmp(movie_type, 'scrambled')
    runID = ['A' run_num];
    load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);
    to_select = intersect(Dmovie_select, Amovies);
    cfg = [];
    cfg.trials = ismember(data.trialinfo, to_select);
    data_A = ft_selectdata(cfg, data);
    
    runID = ['B' run_num];
    load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);
    to_select = intersect(Dmovie_select, Bmovies);
    cfg = [];
    cfg.trials = ismember(data.trialinfo, to_select);
    data_B = ft_selectdata(cfg, data);
    
    runID = ['C' run_num];
    load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);
    to_select = intersect(Dmovie_select, Cmovies);
    cfg = [];
    cfg.trials = ismember(data.trialinfo, to_select);
    data_C = ft_selectdata(cfg, data);
    
    data_comb = ft_appenddata([],data_A, data_B, data_C);
end

if strcmp(movie_type, 'super_intact') || strcmp(movie_type, 'time_warped')
    
    runID = ['D' run_num];
    load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);
    to_select = [1, 9, 10, 12, 20]; % super intact movies
    if strcmp(movie_type, 'time_warped')
        to_select = setdiff(Dmovies, to_select);
    end
    cfg = [];
    cfg.trials = ismember(data.trialinfo, to_select);
    data_D = ft_selectdata(cfg, data);
    
    runID = ['E' run_num];
    load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);
    to_select = [8, 11, 5, 16]; % super intact movies
    if strcmp(movie_type, 'time_warped')
        to_select = setdiff(Emovies, to_select);
    end
    cfg = [];
    cfg.trials = ismember(data.trialinfo, to_select);
    data_E = ft_selectdata(cfg, data);
    
    data_comb = ft_appenddata([], data_D, data_E);
end

end