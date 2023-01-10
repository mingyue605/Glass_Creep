clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS AND PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kob-Andersen 80-20 systems, 10,000 particles
natoms = 1e4;

% Analyze: same initial configuration, 14 runs (sim121), deltaN = 1e4,
%          uniaxial compression with stress = 0.5 in the xx direction

iniconf = 'sim125'; % initial configuration, can choose from 121 to 130

nruns = 15; % from 0 to 13

sim_name = randperm(30,5);

stress = 0.5; % magnitude of the applied stress
stress_tensor = stress*[1 0 0; % [xx xy xz
                        0 0 0; %  xy yy yz
                        0 0 0]; % xz yz zz]

deltaN = 1e4; % number of steps between recorded frames (i.e., cycles)
N = 1e7; % total number of steps
ncycles = N/deltaN; % number of cycles

% the outputs from the simulation are:
% 0 to ncycles files named: #NL.xyz (total number of files = ncycles+1)
% 1 to ncycles files named: #L.xyz (total number of files = ncycles)
% boxsizeL.txt - size of the simulation box from 1L to ncyclesL (total ncycles)
% boxsizeNL.txt - size of the simulation box from 0NL to ncyclesNL (total ncycles+1)
% the reference undeformed configuration is 0NL.xyz

% Inputs for the reference configuration: cycle and L/NL
cycle_ref = 0; % undeformed configuration
type_ref = 'NL'; % can be L or NL (if cycle_ref = 0, it can only be NL)

% Inputs for the reference configuration: cycle and L/NL
cycles = (10:10:1000); % undeformed configuration
type_curr = 'NL'; % can be L or NL (if cycle_ref = 0, it can only be NL)

% path to simulation data
%path_to_runs = strcat('/Volumes/LaCie/Work/compression_creep/deltaN1E4/simulation_data/stress0',num2str(stress*10),'/',iniconf);
path_to_runs = strcat('/Volumes/Mingyue_Plu/deltaN1E4/simulation_data/force',num2str(stress*10),'/',iniconf);

%%% END OF INPUTS AND PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the features of the reference configuration

% load the coordinates of the reference configuration
name_ref = strcat(num2str(cycle_ref),type_ref,'.xyz');
coords_ref = importdata(strcat(path_to_runs,'/run1/',name_ref),' ',2);
atom_type = coords_ref.data(:,1);
coords_ref = coords_ref.data(:,2:4);

% load the simulation box size of the reference configuration
Lbox = importdata(strcat(path_to_runs,'/run1/','boxsize',type_ref,'.txt'));

if strcmp(type_ref,'NL')
    row = cycle_ref + 1;
else
    row = cycle_ref;
end
Lbox_ref = Lbox(row,:);
%%
% Calculate the features of the reference configuration if the file doesn't
% already exists
features_filename_out = strcat('features_',iniconf,'_',num2str(cycle_ref),type_ref,'.csv');
if isfile(features_filename_out)==0

    % short range and medium range order features according to Jain et al.
    disp('Calculating SRO & MRO features')
    F_sromro = sromro_features(coords_ref,atom_type,Lbox_ref); % x 60

    % radial symmetry functions according to Liu et al.
    disp('Calculating radial symmetry function features')
    F_radsymfun = radsymfun_features(coords_ref,atom_type,Lbox_ref); % x 100

    % calculate the orientation matrix features (the rotation matrix of the Voronoi
    % cell principal axis and the stress tensor) | x 9
    disp('Calculating orientation matrix features')
    F_orimat = orimat_features(coords_ref,Lbox_ref,stress_tensor);

    % write out the features to a csv file
    F = [F_sromro F_radsymfun F_orimat];
    F_labels = [cellstr(strcat(repmat('SRO',12,1),num2str((1:12)')));
        cellstr(strcat(repmat('MRO',48,1),num2str((1:48)')));
        cellstr(strcat(repmat('SYM',100,1),num2str((1:100)')));
        cellstr(strcat(repmat('ROT',9,1),num2str((1:9)')))];
    T = array2table(F);
    T.Properties.VariableNames = F_labels;
    writetable(T,features_filename_out);

end


%% Calculate the targets
for iiii=1:numel(cycles)

    cycle_curr = cycles(iiii);

    % Calculate the targets between the reference and current configurations

    dglobal_x_all = zeros(natoms,nruns);
    dglobal_y_all = zeros(natoms,nruns);
    dglobal_z_all = zeros(natoms,nruns);
    
    dglobal_total_all = zeros(natoms,nruns);
    
    d2minlocal_total_all = zeros(natoms,nruns);

    out_oneBig_all = zeros(natoms,nruns);
    out_allBig_all = zeros(natoms,nruns);
    
    Gaussian_fit = zeros()

    for ii=1:nruns % calculate the targets for each independent run

        disp(strcat('Cycle:',{' '},num2str(cycle_curr),'; Run:',{' '},num2str(ii)))

        % load the coordinates of the current configuration
        name_curr = strcat(num2str(cycle_curr),type_curr,'.xyz');
        coords_curr = importdata(strcat(path_to_runs,'/run',num2str(ii),'/',name_curr),' ',2);
        coords_curr = coords_curr.data(:,2:4);

        % load the simulation box size of the current configuration
        Lbox = importdata(strcat(path_to_runs,'/run',num2str(ii),'/boxsize',type_curr,'.txt'));
        if strcmp(type_curr,'NL')
            row = cycle_curr + 1;
        else
            row = cycle_curr;
        end
        Lbox_curr = Lbox(row,:);

        % Calculate the targets between the reference and current

        targets_filename_out = strcat('targets_',iniconf,'_from',num2str(cycle_ref),type_ref,'_to',num2str(cycle_curr),type_curr,'_run',num2str(ii),'.csv');

        % Calculate non-affine displacements with respect to global deformations
        disp('Calculating non-affine global displacements')
        [dglobal_components, dglobal_total] = global_non_affine(coords_ref,coords_curr,Lbox_ref,Lbox_curr);

        % Calculate local non-affine displacements according to Falk & Lagner
        disp('Calculating non-affine local displacements')
        [~,dsquare_min_norm] = local_non_affine(coords_ref,coords_curr,Lbox_ref,Lbox_curr);

        % Calculate atoms that rearranged
        disp('Calculating atoms that rearranged')
        [out_oneBig,out_allBig]= outlier_detector(dglobal_components);

        % write targets to a csv file
        targets = [dglobal_components dglobal_total dsquare_min_norm out_oneBig out_allBig];
        targets_labels = {'dglobal_x','dglobal_y','dglobal_z','dglobal_total','d2local_total','out_oneBig','out_allBig'};
        T = array2table(targets);
        T.Properties.VariableNames = targets_labels;
        writetable(T,targets_filename_out);

        % Store the targets in global variables
        dglobal_x_all(:,ii) = dglobal_components(:,1);
        dglobal_y_all(:,ii) = dglobal_components(:,2);
        dglobal_z_all(:,ii) = dglobal_components(:,3);
        dglobal_total_all(:,ii) = dglobal_total;
        d2minlocal_total_all(:,ii) = dsquare_min_norm;
        out_oneBig_all(:,ii) = out_oneBig;
        out_allBig_all(:,ii) = out_allBig;

    end

    % calculate and save the averages and standard deviations over
    % all the independent runs

    targets_summary_filename_out = strcat('targets_',iniconf,'_from',num2str(cycle_ref),type_ref,'_to',num2str(cycle_curr),type_curr,'.csv');
    if isfile(targets_summary_filename_out)==0

        % write targets to a csv file
        targets = [mean(dglobal_x_all,2), mean(dglobal_y_all,2), mean(dglobal_z_all,2),...
            mean(dglobal_total_all,2), mean(d2minlocal_total_all,2),...
            mean(out_oneBig_all,2), mean(out_allBig_all,2),...
            std(dglobal_x_all,0,2), std(dglobal_y_all,0,2), std(dglobal_z_all,0,2),...
            std(dglobal_total_all,0,2), std(d2minlocal_total_all,0,2),...
            std(out_oneBig_all,0,2), std(out_allBig_all,0,2)];

        targets_labels = {'avg_dglobal_x','avg_dglobal_y','avg_dglobal_z',...
            'avg_dglobal_total','avg_d2local_total',...
            'avg_out_oneBig','avg_out_allBig',...
            'std_dglobal_x','std_dglobal_y','std_dglobal_z',...
            'std_dglobal_total','std_d2local_total',...
            'std_out_oneBig','std_out_allBig'};

        T = array2table(targets);
        T.Properties.VariableNames = targets_labels;
        writetable(T,targets_summary_filename_out);

    end

end


