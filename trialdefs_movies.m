function [trl_new] = trialdefs_movies(cfg)
% trialdefs_movie.m
% Input: cfg.dataset = mat file containing ECoG data

% Output: trl
% trl -> [begin end offset 'MovieName']
% dimension of trl -> #'S  4' trials X 4

% read the header information and the events from the data
%load(cfg.dataset);
% loads the following
% 1. hdr (header)
% 2. raw (data from all channels)

% List of movie specific markers

%allmovieFileName  = ['M01';'M02';'M03';    'M04';  'M05';  'M06';  'M07';  'M08';  'M10';  'M11';  'M12';  'M13'; ...
%        'M14';  'M15';  'M16';  'M17';  'M18';  'M19';  'M21';  'M22';  ...
%        'S01';  'S02';  'S03';  'S04';  'S05';  'S06';  'S07';  'S08';  'S10';  'S11';  'S12';  'S13'; ...
%        'S14';  'S15';  'S16';  'S17';  'S18';  'S19';  'S21';  'S22'];


% search for "Stimulus" events
if (strcmp(cfg.runID(1),'A'))
    
    if(cfg.patientNum == 2)
        movie_begin_triggers = [15093, 97382, 167106, 214825, 265542, 308345, 355980, 444503];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [15630, 97919, 167643, 215362, 266079, 308882, 356517, 445040];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [17450, 99739, 169463, 217182, 267899, 310702, 358337, 446860];
    end
    
    % Plot channel raw(63,:) and detect the triggers manually. <Only for beginning of movie>    
    % movieFileName = ['S03';'M04';'S12';'M11';'M22';'S17';'S05';'M08'];
    movieTags  = [23, 4, 31, 10, 20, 36, 25, 8]; % check allmovieFileName for indices 
    % Used mat2fieldtrip.m and AP_ECoG_erp_trialdefs.m to create the
    % following (changed start_index = start_index + 1 -> start_index to
    % include first trigger
end

if (strcmp(cfg.runID(1),'B'))
    
    if(cfg.patientNum == 2)
        movie_begin_triggers = [16501, 62217, 178545, 233424, 288427, 374979, 418782, 491574];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [38420, 84136, 200464, 255343, 310346, 396898, 440701, 513493];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [17060, 62776, 179104, 233983, 288986, 375538, 419341, 492133];
    end
    
    %movieFileName = ['M12';'S08';'S01';'M13';'M05';'S22';'M10';'S07'];
    movieTags = [11, 28, 21, 12, 5, 40, 9, 27];
end

if (strcmp(cfg.runID(1),'C'))
    
    if(cfg.patientNum == 2)
        movie_begin_triggers = [14705, 102509, 153226, 210106, 281834, 355627, 403280, 458273];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [19180, 106984, 157701, 214581, 286309, 360102, 407755, 462748];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [15440, 103244, 153961, 210841, 282569, 356362, 404015, 459008];
    end
    
    %movieFileName = ['M07';'S11';'M01';'S04';'S10';'M17';'S13';'M03'];
    movieTags = [7, 30, 1, 24, 29, 16, 32, 3];
end


% Run D and E for patient 3 are with warped timeshift movies
if (strcmp(cfg.runID(1),'D'))
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [17507, 77705, 175990, 292540, 366350, 417065, 477175, 533590, 589425 ];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [16412, 76610, 174895, 291445, 365255, 415970, 476080, 532495, 588330];
    end
    
    %movieFileName = ['M01';'M05';'M08';'M09';'M10';'M11';'M12';'M16', 'M20'];
    movieTags = [1, 5, 8, 9, 10, 11, 12, 16, 20];
end


if (strcmp(cfg.runID(1),'E'))
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [16985, 132315, 178180, 242025, 329780, 409750, 501740, 551385, 614150];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [17730, 133060, 178925, 242770, 330525, 410495, 502485, 552130, 614895];
    end
    
    %movieFileName = ['M08';'M11';'M20';'M05';'M01';'M09';'M16';'M12', 'M10'];
    movieTags = [8, 11, 20, 5, 1, 9, 16, 12, 10];
end

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * cfg.Fs);
posttrig =  round(cfg.trialdef.post * cfg.Fs);

n_movies = length(movieTags);
trl_new = zeros(n_movies,4);

for i = 1:n_movies
    trl_new(i,1) = movie_begin_triggers(i) + pretrig;
    trl_new(i,2) = movie_begin_triggers(i) + posttrig;
    trl_new(i,3) = pretrig;
    trl_new(i,4) = movieTags(i); 
end