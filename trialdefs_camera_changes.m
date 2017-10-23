function [trl_new] = trialdefs_camera_changes(cfg)
% AP_ECoG_trialdefs.m
% Input: cfg.dataset = mat file containing ECoG data

% Output: trl
% trl -> [begin end offset 'MovieName']
% dimension of trl -> #'S  4' trials X 4

% read the header information and the events from the data
%load(cfg.dataset);
% loads the following
% 1. hdr (header)
% 2. raw (data from all channels)


% search for "Stimulus" events
if (strcmp(cfg.runID(1),'A'))
    numCameraChanges = [40,42,34,28,20,35,39,39];
    % These movie numbers are given based on its position in the original
    % eeg Runs. All M movies are < 20 and all S movies > 20. All movies are
    % not included in the ECoG runs.
    
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
    camera_trigger_file = 'runA_trl.mat';
end

if (strcmp(cfg.runID(1),'B'))
    numCameraChanges = [34,39,27,27,39,20,38,39];
    
    if(cfg.patientNum == 2)
        movie_begin_triggers = [16501, 62217, 178545, 233424, 288427, 374979, 418782, 491574];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [38420, 84136, 200464, 255343, 310346, 396898, 440701, 513493];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [17060, 62776, 179104, 233983, 288986, 375538, 419341, 492133];
    end
    
    
    %movieFileName = ['M12';'S08';'S01';'M01';'M05';'S22';'M10';'S07'];
    movieTags = [11, 28, 21, 1, 5, 40, 9, 27];
    camera_trigger_file = 'runB_trl.mat';
end

if (strcmp(cfg.runID(1),'C'))
    numCameraChanges = [39,28,27,42,38,35,27,40];
    
    if(cfg.patientNum == 2)
        movie_begin_triggers = [14705, 102509, 153226, 210106, 281834, 355627, 403280, 458273];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '1')
        movie_begin_triggers = [19180, 106984, 157701, 214581, 286309, 360102, 407755, 462748];
    end
    
    if(cfg.patientNum == 3 && cfg.runID(2) == '2')
        movie_begin_triggers = [15440, 103244, 153961, 210841, 282569, 356362, 404015, 459008];
    end
    
    
    movieFileName = ['M07';'S11';'M01';'S04';'S10';'M17';'S13';'M03'];
    movieTags = [7, 30, 1, 24, 29, 16, 32, 3];
    camera_trigger_file = 'runC_trl.mat'; % load trl from when we had all triggers.
end


movie_begin_triggers = movie_begin_triggers - 100; % correction because above corresponds to middle of trigger
load(camera_trigger_file); % load trl
trigger_indices = trl(:,1)-trl(:,3);
triggers = movie_begin_triggers(1) + trigger_indices - trigger_indices(1);

% Used numCameraChanges to test
% For example in the case of runID == 'C', 
% triggers(numCameraChanges(1)+1) = movie_begin_triggers(2) % nearly equal

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * cfg.Fs);
posttrig =  round(cfg.trialdef.post * cfg.Fs);

trl_new = trl;

trl_new(:,1) = triggers + pretrig;
trl_new(:,2) = triggers + posttrig;
trl_new(:,3) = pretrig; % offset
% trl_new(:,4) are the labels of the movies which remain the same as trl 


