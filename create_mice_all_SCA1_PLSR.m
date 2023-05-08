% create_mice_all_SCA1_PLSR.m
% Sarah West
% 5/5/23

clear all;

% Create the experiment name.
parameters.experiment_name ='SCA1 PLSR';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\'];  

mice_all(1).name = '007';
mice_all(1).number_of_sources = 14;
mice_all(1).date = '07-01-2021';

mice_all(2).name = '819';
mice_all(2).number_of_sources = 20;
mice_all(2).date = '07-02-2021';

save([parameters.dir_exper 'mice_all.mat'], 'mice_all');