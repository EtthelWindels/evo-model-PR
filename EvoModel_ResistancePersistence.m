%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to simulate bacterial evolution of resistance and persistence under intermittent antibiotic treatment

% Authors: Etthel M. Windels & Lloyd Cool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
%tic
nsim = 1000;


%%% TIME

time_step = 30; % [minutes]
cycles = 10;

total_time = 24*60*cycles;
time_seq = 0:time_step:total_time;
nr_timesteps = length(time_seq);
times_treatment = time_seq(mod(time_seq,24*60)==0 & (time_seq)>=24*60);
times_treatment = num2cell(times_treatment(:));
times_growth = time_seq(mod(time_seq,24*60)==5*60 & (time_seq)>=24*60);
times_growth = num2cell(times_growth(:));


%%% MIC AND PERSISTENCE SWITCHING RATES

log_MIC = 1:6;                              % log2(MIC) for different resistance mutant classes
MIC = 2.^log_MIC;
log_a = [-4.5 -3.5 -2.5 -1.5 -0.5];         % log10(forward persistence switching rate) for different persistence mutant classes
R_classes = length(log_MIC);
P_classes = length(log_a);
nr_classes = R_classes*P_classes;
b_g = repelem(log10(0.1),nr_classes);       % log10(backward persistence switching rate) during growth
b_t = repelem(log10(0.1),nr_classes);       % log10(backward persistence switching rate) during treatment
log_MIC_matrix = ones(length(log_a),length(log_MIC)).*log_MIC;
MIC_matrix = ones(P_classes,R_classes).*MIC;
MIC_vec = reshape(MIC_matrix,1,[]);
log_a_matrix = ones(R_classes,P_classes).*log_a;
phenotypes = cell(1,2);
phenotypes{1,1} = reshape(log_a_matrix.',1,[]);
phenotypes{1,2} = reshape(log_MIC_matrix,1,[]);


%%% MUTATION

p_mut = 1e-4;                                   % overall mutation rate
mut_prob = DFEmaker(R_classes,P_classes,p_mut); % mutation probabilities to different classes
mut_prob_matrix = zeros(nr_classes,nr_classes);

i = 1;
for r = 1:R_classes
    for p = 1:P_classes
        int = mut_prob{r,p};
        mut_prob_matrix(i,:) = reshape(int.',1,[]);
        i = i+1;
    end
end


%%% GROWTH

Gmax = ones(1,nr_classes)*0.6*time_step/60;  % maximum growth rate (log10(cell count change per hour))
Gmin = ones(1,nr_classes)*-10*time_step/60;   % minimum growth rate (log10(cell count change per hour))
k_ab = 1;                                    % Hill coefficient for antibiotic dose-response curve
k_nutr = 1;                                  % Hill coefficient for nutrient dose-response curve


%%% TREATMENT

AB_levels = [3.125 6.25 12.5 25 50 100 200 400];        % antibiotic treatment concentrations
nutr_levels = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];  % nutrient concentrations during treatment
AB_low = 0;


%%% INITIAL POPULATION

N0 = 1000;  % initial number of normal cells
P0 = 1;     % initial number of persister cells
K = 1e6;    % carrying capacity


%%% INITIALIZATION OF DATA ARRAY

storage_time_points = cell(1,6,length(AB_levels),length(nutr_levels));
N_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');    % number of normal cells per class per time point
P_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');    % number of persister cells per class per time point
D_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');    % number of dead cells per class per time point

AB_index = 1;
for AB = AB_levels
    nutr_index = 1;
    for nutr = nutr_levels
        storage_time_points{1,1,AB_index,nutr_index} = time_seq;
        storage_time_points{1,2,AB_index,nutr_index} = N_all_sim;
        storage_time_points{1,3,AB_index,nutr_index} = P_all_sim;
        storage_time_points{1,4,AB_index,nutr_index} = D_all_sim;
        storage_time_points{1,5,AB_index,nutr_index} = AB;
        storage_time_points{1,6,AB_index,nutr_index} = nutr;
        nutr_index = nutr_index+1;
    end
    AB_index = AB_index+1;
end


%%% EVOLUTION

AB_index = 1;
for AB_cond = AB_levels
    %disp(['AB_cond:',num2str(AB_cond)])

    nutr_index = 1;
    for nutr_cond = nutr_levels
        %disp(['nutr_cond:',num2str(nutr_cond)])
        N_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');
        P_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');
        D_all_sim = zeros(nr_classes,nr_timesteps+1,nsim,'single');

        parfor sim_nr = 1:nsim
            N_one_sim = zeros(nr_classes,nr_timesteps+1,'single');
            N_one_sim(1,1) = N0;
            D_one_sim = zeros(nr_classes,nr_timesteps+1,'single');
            P_one_sim = zeros(nr_classes,nr_timesteps+1,'single');
            P_one_sim(1,1) = P0;
            gA = zeros(nr_classes,1); % growth rate per class
            nutr_t = 1;
            AB_t = AB_low;
            b = b_g;
            time_index = 1;
            rng(sim_nr) % seed

            for t = time_seq
                [nutr_t,AB_t,N_one_sim(:,time_index),P_one_sim(:,time_index),D_one_sim(:,time_index),b] = update_nutr_AB(t,nutr_t,AB_t,AB_cond,AB_low,nutr_cond,times_treatment,times_growth,N_one_sim(:,time_index),P_one_sim(:,time_index),D_one_sim(:,time_index),b,b_g,b_t);
                %disp(['AB_t:',num2str(AB_t)])
                %disp(['nutr_t:',num2str(nutr_t)])
                total_viable_cells = sum(N_one_sim(:,time_index))+sum(P_one_sim(:,time_index));

                if sum(total_viable_cells)>0

                    % Growth and death
                    gA_log = Gmax-((Gmax-Gmin).*(AB_t./MIC_vec).^k_ab)./((AB_t./MIC_vec).^k_ab-Gmin./Gmax);
                    pos = gA_log>0;
                    neg = gA_log<=0;
                    gA(pos) = log2(10.^(gA_log(pos))).*nutr_t;
                    gA(neg) = 1-10.^(gA_log(neg)*(1/(1+(0.5/nutr_t)^k_nutr)));
                    change_in_cell_count = mybinornd(N_one_sim(:,time_index).',abs(gA));
                    change_in_cell_count(isnan(change_in_cell_count)) = 0;
                    D_one_sim(neg,time_index+1) = D_one_sim(neg,time_index)+change_in_cell_count(neg).';
                    N_one_sim(:,time_index+1) = N_one_sim(:,time_index);
                    N_one_sim(neg,time_index+1) = N_one_sim(neg,time_index+1)-change_in_cell_count(neg).';

                    % Mutation
                    change_in_cell_count_sign = change_in_cell_count.*sign(gA_log);
                    prob_index = 1;
                    for cells = change_in_cell_count_sign
                        if cells > 0
                            probability = mut_prob_matrix(prob_index,:)/sum(mut_prob_matrix(prob_index,:));
                            N_one_sim(:,time_index+1) = N_one_sim(:,time_index+1)+mnrnd(cells,probability).';
                        end
                        prob_index = prob_index+1;
                    end

                    % Phenotype switching
                    forward_rates = 10.^(phenotypes{1})*(1-nutr_t)*time_step/60;
                    backward_rates = 10.^(b)*nutr_t*time_step/60;
                    cells_to_persister = mybinornd(N_one_sim(:,time_index+1).',forward_rates);
                    cells_to_normal = mybinornd(P_one_sim(:,time_index).',backward_rates);
                    N_one_sim(:,time_index+1) = N_one_sim(:,time_index+1)-cells_to_persister.'+cells_to_normal.';
                    P_one_sim(:,time_index+1) = P_one_sim(:,time_index)+cells_to_persister.'-cells_to_normal.';

                    % Update nutrient concentration
                    nutr_t = max(nutr_t-sum(change_in_cell_count_sign(change_in_cell_count_sign>0))/K,0);
                    if nutr_t<0
                        nutr_t = 0;
                    end
                end
                time_index = time_index+1;
            end
            N_all_sim(:,:,sim_nr) = N_one_sim;
            P_all_sim(:,:,sim_nr) = P_one_sim;
            D_all_sim(:,:,sim_nr) = D_one_sim;
        end
        storage_time_points{1,2,AB_index,nutr_index} = N_all_sim;
        storage_time_points{1,3,AB_index,nutr_index} = P_all_sim;
        storage_time_points{1,4,AB_index,nutr_index} = D_all_sim;
        nutr_index = nutr_index+1;
    end
    AB_index  = AB_index+1;
end
%toc


%%% DATA STORAGE

datafr = ones(length(nutr_levels)*length(AB_levels),5);
datafr_full = ones(length(nutr_levels)*length(AB_levels)*nsim,5);

i = 1;
k = 1;
for AB_index = AB_levels
    j = 1;
    for b = nutr_levels
        total = storage_time_points{1,2,k,j}+storage_time_points{1,3,k,j};
        percentage_extinct = 1-sum(sum(total(:,nr_timesteps+1,:))>0)/nsim;
        mic_av= nanmean((MIC_vec*squeeze(total(:,nr_timesteps+1,:)))./(squeeze(sum(total(:,nr_timesteps+1,:))).'));
        pers_av = nanmean((phenotypes{1}*squeeze(total(:,nr_timesteps+1,:)))./squeeze(sum(total(:,nr_timesteps+1,:))).');
        aB = storage_time_points{1,5,k,j};
        nutr = storage_time_points{1,6,k,j};
        datafr(i,1) = aB;
        datafr(i,2) = nutr;
        datafr(i,3) = percentage_extinct;
        datafr(i,4) = mic_av;
        datafr(i,5) = pers_av;

        mic = (MIC_vec*squeeze(total(:,nr_timesteps+1,:)))./(squeeze(sum(total(:,nr_timesteps+1,:))).');
        pers = (phenotypes{1}*squeeze(total(:,nr_timesteps+1,:)))./squeeze(sum(total(:,nr_timesteps+1,:))).';
        datafr_full((nsim*(i-1)+1):(nsim*i),1) = 1:nsim;
        datafr_full((nsim*(i-1)+1):(nsim*i),2) = aB;
        datafr_full((nsim*(i-1)+1):(nsim*i),3) = nutr;
        datafr_full((nsim*(i-1)+1):(nsim*i),4) = mic;
        datafr_full((nsim*(i-1)+1):(nsim*i),5) = pers;
        i = i+1;
        j = j+1;
    end
k = k+1;
end

writematrix(datafr,'PRevo_out_average.csv')
writematrix(datafr_full,'PRevo_out_full.csv')


%%% MUTATION PROBABILITIES

function dfe = DFEmaker(R_classes,P_classes,p_mut)
    % Calculate mutation probabilities from each mutant class to all other classes.
    dfe = cell(R_classes,P_classes);
    tmp1 = zeros(R_classes,P_classes);
    tmp2 = repelem(0,P_classes);
    for t = (0:P_classes-1)
        tmp2(t+1) = 10*0.1^t;
    end
    tmp1(1,1:P_classes) = tmp2;
    row = 1;
    while row < R_classes
        row = row + 1;
        tmp1(row,:) = tmp1(row-1,:)*0.1;
    end
    tmp1 = tmp1/(sum(tmp1,'all')-tmp1(1,1))*p_mut;
    tmp1(1,1) = 1-p_mut;
    for i = 1:R_classes
        for j = 1:P_classes
            mat = zeros(R_classes,P_classes);
            if i==R_classes && j==P_classes
                mat(i:R_classes,j:P_classes) = 1;
            else
                mat(i:R_classes,j:P_classes) = tmp1(1:(R_classes-i+1),1:(P_classes-j+1));
                mat(i,j) = 1-(sum(mat,'all')-mat(i,j));
            end
            dfe{i,j} = mat;
        end
    end
end


%%% GROWTH MEDIUM CHANGES

function [nutr_t,AB_t,N_cells,P_cells,D_cells,b] = update_nutr_AB(t,nutr_t,AB_t,AB_cond,AB_low,nutr_cond,times_treatment,times_growth,N_cells,P_cells,D_cells,b,b_g,b_t)
    % Update population sizes of all classes at start of treatment phase and start of growth phase, based on corresponding dilution factors and medium volume changes.
    switch t
        case times_treatment
            nutr_t = nutr_cond;
            dilution = 1-nutr_cond;
            AB_t = AB_cond;
            b = b_t;

            total_cells = sum(N_cells)+sum(P_cells)+sum(D_cells);
            total_transferred = round(total_cells*200/500*dilution); % medium volume reduces from 500 to 200 ul
            transferred_per_type = mnrnd(total_transferred,[sum(N_cells) sum(P_cells) sum(D_cells)]/total_cells).';
            if sum(N_cells)~=0
                prob_n = (N_cells/sum(N_cells)).';
                N_cells(prob_n>0) = mnrnd(transferred_per_type(1),prob_n(prob_n>0));
            end
            if sum(P_cells)~=0
                prob_p = (P_cells/sum(P_cells)).';
                P_cells(prob_p>0) = mnrnd(transferred_per_type(2),prob_p(prob_p>0));
            end
            if sum(D_cells)~=0
                prob_d = (D_cells/sum(D_cells)).';
                D_cells(prob_d>0) = mnrnd(transferred_per_type(3),prob_d(prob_d>0));
            end
            N_cells(isnan(N_cells)) = 0;
            P_cells(isnan(P_cells)) = 0;
            D_cells(isnan(D_cells)) = 0;

        case times_growth
            nutr_t = 1;
            AB_t = AB_low;
            b = b_g;

            total_cells = sum(N_cells)+sum(P_cells)+sum(D_cells);
            total_transferred = round(total_cells*0.01*500/200); % 1:100 dilution after treatment cycle; medium volume increases from 200 to 500 ul
            transferred_per_type = mnrnd(total_transferred,[sum(N_cells) sum(P_cells) sum(D_cells)]/total_cells).';
            if sum(N_cells)~=0
                prob_n = (N_cells/sum(N_cells)).';
                N_cells(prob_n>0) = mnrnd(transferred_per_type(1),prob_n(prob_n>0));
            end
            if sum(P_cells)~=0
                prob_p = (P_cells/sum(P_cells)).';
                P_cells(prob_p>0) = mnrnd(transferred_per_type(2),prob_p(prob_p>0));
            end
            if sum(D_cells)~=0
                prob_d = (D_cells/sum(D_cells)).';
                D_cells(prob_d>0) = mnrnd(transferred_per_type(3),prob_d(prob_d>0));
            end
            N_cells(isnan(N_cells)) = 0;
            P_cells(isnan(P_cells)) = 0;
            D_cells(isnan(D_cells)) = 0;
        otherwise
    end
end



%%% BINOMIAL DISTRIBUTION

function ret = mybinornd(n,p)
  % Given vectors of binomial parameters, generate draws from corresponding binomial distributions.
    ret = zeros(1,length(n));
    for k = 1:length(n)
        ret(k) = binornd(n(k),p(k));
    end
    check = ret>n;
    if sum(check) >= 1
        ret(check) = n(check);
    end
end