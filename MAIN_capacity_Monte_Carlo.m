%% Estimates the potential variability in scanner capacity enhancement with DRB using preliminary data in sequence variation and MC model
% Mikael Brix, 6.12.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
% Define file paths
mskFile = '/Users/mikaeljuntunen/Documents_Mikael/Studies/KTM/Masters thesis/DATA/sekvenssikestot_KOOSTE_TODELLISET_AJAT.xlsx';
examinationsFile = '/Users/mikaeljuntunen/Documents_Mikael/Studies/KTM/Masters thesis/DATA/Kuntaliittokoodit_ja_kestot_2022_TODELLISET_AJAT.xlsx';
acceleration_testData = '/Users/mikaeljuntunen/Documents_Mikael/Studies/KTM/Masters thesis/DATA/Kiihdytys.mat';
saveLoc = '/Users/mikaeljuntunen/Documents_Mikael/Studies/KTM/Masters thesis/Results/1000_MC_iterations';

fontsize = 12;
fontname = 'Times new roman';
% Simulation parameters
useGlobal = 1; % 1 = use global value for drb acceleration (obtained as a mean of test subjects), if 0, sequence specific speed-up estimates are used
positioning_time_minimum = 13;%7.5;10;15
num_iter = 1000; % number of monte carlo iteration

%% Initialize speed-up estimates
if useGlobal
    % Using global acceleration factor
    speed_up = struct;
    speed_up.t2_sag_tse.mean = 0.597;
    speed_up.t2_sag_tse.median = 0.608;
    speed_up.t2_sag_tse.sd = 0.153;

    speed_up.stir = speed_up.t2_sag_tse;
    speed_up.pd_sag_tse = speed_up.t2_sag_tse;

    speed_up.pd_cor_tse = speed_up.t2_sag_tse;
    speed_up.t1_sag_tse = speed_up.t2_sag_tse;
    speed_up.pd_tra_tse = speed_up.t2_sag_tse;
    speed_up.t2_tra_tse = speed_up.t2_sag_tse;
    speed_up.t1_tra_tse = speed_up.t2_sag_tse;
    speed_up.t1_cor_tse = speed_up.t2_sag_tse;

    speed_up.t2_tra_blade = speed_up.t2_sag_tse;
    speed_up.t2_cor_tse = speed_up.t2_sag_tse;

else
    % Using sequence specific acceleration factors
    speed_up = struct;
    speed_up.t2_sag_tse.mean = 0.697;
    speed_up.t2_sag_tse.median = 0.748;
    speed_up.t2_sag_tse.sd = 0.066;

    speed_up.stir.mean = 0.547;
    speed_up.stir.median = 0.617;
    speed_up.stir.sd = 0.19;

    speed_up.pd_sag_tse.mean = 0.592;
    speed_up.pd_sag_tse.median = 0.592;
    speed_up.pd_sag_tse.sd = 0.2;

    speed_up.pd_cor_tse.mean = 0.756;
    speed_up.pd_cor_tse.median = 0.756;
    speed_up.pd_cor_tse.sd = 0.0;

    speed_up.t1_sag_tse.mean = 0.401;
    speed_up.t1_sag_tse.median = 0.375;
    speed_up.t1_sag_tse.sd = 0.042;

    speed_up.pd_tra_tse.mean = 0.536;
    speed_up.pd_tra_tse.median = 0.545;
    speed_up.pd_tra_tse.sd = 0.018;

    % Added without statistics
    speed_up.t2_tra_tse.mean = 0.6;
    speed_up.t2_tra_tse.median = 0.6;
    speed_up.t2_tra_tse.sd = 0.018;


    speed_up.t1_tra_tse.mean = 0.401;
    speed_up.t1_tra_tse.median = 0.375;
    speed_up.t1_tra_tse.sd = 0.042;
    % Added part ends

    speed_up.t1_cor_tse.mean = 0.617;
    speed_up.t1_cor_tse.median = 0.664;
    speed_up.t1_cor_tse.sd = 0.19;

    speed_up.t2_tra_blade.mean = 0.612;
    speed_up.t2_tra_blade.median = 0.600;
    speed_up.t2_tra_blade.sd = 0.089;

    speed_up.t2_cor_tse.mean = 0.71;
    speed_up.t2_cor_tse.median = 0.768;
    speed_up.t2_cor_tse.sd = 0.101;
end
DRB_sequences = fieldnames(speed_up);
%% Optional, visualize speed-up distribution
speed_ups = random('normal', speed_up.t2_sag_tse.mean,  speed_up.t2_sag_tse.sd,1000,1); 
speed_ups(speed_ups>1) = [];
speed_ups(speed_ups<0) = [];
figure;
histogram(speed_ups*100,25,'FaceColor','k');
set(gca,'fontsize',fontsize+2,'fontname',fontname)
xlabel('Speed-up (%)','fontsize',fontsize+4,'fontname',fontname);
ylabel('','fontsize',fontsize+4,'fontname',fontname);
% % Anderson-Darling test: returns a test decision for the null hypothesis that the data in vector x is from a population with a normal distribution, using the Anderson-Darling test. 
% [h,p] = adtest(Kiihdytys); % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, or 0 otherwise.
% skew = skewness(Kiihdytys) = -1.3142 --> right-skewed
% kurt = kurtosis(Kiihdytys) = 5.0733 > 3 --> more outlier prone
load(acceleration_testData);
Kiihdytys(Kiihdytys<20 | isnan(Kiihdytys)) = [];
[transdat,lambda] = boxcox(Kiihdytys);% transdat=(data^lambda−1)/lambda 
% Kiihdytys = (lambda*transdat+1).^(1/lambda)
% figure;subplot(1,3,1);hist(Kiihdytys,10);xlim([0,100]);title('Original (%)');
% subplot(1,3,2);hist(speed_ups*100,10);xlim([0,100]);title('Normal distribution (%)');
% subplot(1,3,3);conv_data = (lambda*speed_ups*100+1).^(1/lambda);hist(conv_data*mean(speed_ups*100)/mean(conv_data),10);xlim([0,100]);title('Transformed Normal distribution');
% Determine the cumulative density function
CDF = cumsum(hist(Kiihdytys,20));
% ChatGPT suggests that this is either beta or gamma-distributed
pd = fitdist(Kiihdytys,'gamma');
xrange = 1:100;
figure;
histogram(Kiihdytys,100,'Normalization','Probability');
xlim([0,100]);hold on;
plot(xrange,pdf(pd,xrange));


pd = fitdist(Kiihdytys/100,'beta');
xrange = (1:100)/100;
figure;
histogram(Kiihdytys/100,20,'Normalization','Probability');
xlim([0,1]);hold on;
plot(xrange,pdf(pd,xrange)*5/100);

%% Import msk data
msk_times = importMSKfile(mskFile, 'MSK-protokollat todelliset ajat');
msk_times.Aika = msk_times.Aika*60*24;
msk_times_corrected = msk_times(1,:);
% Remove lyhenne with numbers - these indicate unwanted protocols
loop = 1;
for ii = 1:size(msk_times,1)
    idx = regexp(string(msk_times.Lyhenne(ii)),'\d*','Match');
    if isempty(idx)
        msk_times_corrected(loop,:) = msk_times(ii,:);
        loop = loop + 1;
    end
end
%% Import neuro data
neuro_times = importNeurofile(mskFile, 'Neuro todelliset ajat');
neuro_times.Aika = neuro_times.Aika*60*24;
neuro_times_corrected = neuro_times(1,:);
% Remove lyhenne with numbers - these indicate unwanted protocols
loop = 1;
for ii = 1:size(neuro_times,1)
    idx = regexp(string(neuro_times.Lyhenne(ii)),'\d*','Match');
    if isempty(idx)
        neuro_times_corrected(loop,:) = neuro_times(ii,:);
        loop = loop + 1;
    end
end
%% Import body

%% Combine
times_total = cat(1,msk_times_corrected,neuro_times_corrected);
%%
examination_numbers = importTutkimuskestot(examinationsFile, 'Statistics');
results = [];
%% Run the Monte carlo
slots = [20,30,45,60];
slot_identifiers = {'A','B','C','D'};
times_DRB = [];
for ii = 2:size(examination_numbers,1)
    fprintf('Processing examination (%s) %d of %d\n',examination_numbers.CPTCode(ii),ii,size(examination_numbers,1));
    % Identify the current appointment slot - check the fourth CPT code
    % element
    tempname = char(examination_numbers.CPTCode(ii));
    idx_slot = find(~cellfun(@isempty,strfind(slot_identifiers,tempname(4))));
    
    idxs = find(~cellfun(@isempty,strfind(string(times_total.Tunniste),examination_numbers.CPTCode(ii))));
    sequence_table = times_total(idxs,:);
    results{ii-1,1} = examination_numbers.CPTCode(ii);
    results{ii-1,2} = examination_numbers.KuvaustenLkm(ii);
    for iter = 1:num_iter 
        DRB_total_time = 0;
        routine_total_time = 0;
        for ii_sequence = 1:size(sequence_table,1)
            if isempty(strfind(string(sequence_table.Sekvenssi(ii_sequence)),'(')) & isempty(strfind(string(sequence_table.Sekvenssi(ii_sequence)),'Kokonaisaika'))
                sequence_name = string(sequence_table.Sekvenssi(ii_sequence));
                sequence_name = lower(sequence_name);
                sequence_name = strrep(sequence_name,' ','_'); % Replace space with null.
                idx = [];
                for ii2 = 1:numel(DRB_sequences)
                    if strfind(sequence_name,DRB_sequences{ii2})
                        idx = ii2;
                    end
                end
                if ~isempty(idx)
                    % Estimate the potential speed-up with monte carlo
                    speed_up_estimate = random('normal', speed_up.(DRB_sequences{idx}).mean,  speed_up.(DRB_sequences{idx}).sd); 
                    % Added on 17.12.2023
                    if speed_up_estimate > 0.9
                        speed_up_estimate = 0.9;
                    end
                    if speed_up_estimate < 0
                        speed_up_estimate = 0;
                    end
                else
                    speed_up_estimate = 0; % no speed up
                end
               % fprintf('Sequence: %s, speed-up: %4.1f\n',sequence_name,speed_up_estimate);
    %             disp(speed_up_estimate);
                routine_total_time = routine_total_time + sequence_table.Aika(ii_sequence);
                DRB_total_time = DRB_total_time + sequence_table.Aika(ii_sequence)*(1-speed_up_estimate);
            end
        end
        times_DRB(ii-1,iter) = DRB_total_time;
        % Criteria for moving to a shorter appointment slot
        time_reduction = routine_total_time - DRB_total_time;
        if idx_slot > 1 % No shorter slot than 20 min
            positioning_time = slots(idx_slot-1) - DRB_total_time;
            if time_reduction >= (slots(idx_slot)-slots(idx_slot-1)) | positioning_time >= positioning_time_minimum
                results{ii-1,iter+2} = examination_numbers.KuvaustenLkm(ii)*1.5;
            else
                results{ii-1,iter+2} = examination_numbers.KuvaustenLkm(ii)*1;
            end
        else
            results{ii-1,iter+2} = examination_numbers.KuvaustenLkm(ii)*1;
        end
    end
end
% %% Visualize confidence interval plot
% x = 1:size(results,1);
% y = cell2mat(results(1:end,3:end))/sum(cell2mat(results(1:end,2)));
% avg_data = mean(y,2)';
% std_data = std(y,0,2)';
% fill([x, flip(x)], [avg_data+std_data, flip(avg_data-std_data)], [0.8 0.8 0.8])
% hold on
% plot(x, avg_data, 'k','linewidth',2)
%% Visualize confidence interval plot of cumulative sum
cumSum_noDRB = cumsum(cell2mat(results(1:end,2)),1);
x = 1:size(results,1);
% y = cumsum(cell2mat(results(1:end,3:end))/sum(cell2mat(results(1:end,2)),1))*100;
%% Original cumsum method
y = cumsum(cell2mat(results(1:end,3:end)))./cumSum_noDRB*100;
avg_data = mean(y,2)';
std_data = std(y,0,2)';
%% Adapted method
% Alternative method
y = cell2mat(results(1:end,3:end));
yCorr = y;
for ii = 1:size(yCorr)
    if ii < size(yCorr)
        yCorr(ii,:)= sum(y(1:ii,:),1) + sum(cell2mat(results(ii+1:end,2)));
    else
        yCorr(ii,:)= sum(y(1:ii,:),1);
    end
end
yCorr = yCorr/sum(cell2mat(results(1:end,2)),1)*100;
avg_data = mean(yCorr,2)';
std_data = std(yCorr,0,2)';
figure;
fill([x, flip(x)], [avg_data+std_data, flip(avg_data-std_data)], [0.8 0.8 0.8])
hold on
plot(x, avg_data, 'k','linewidth',2)
xlim([0,140]);
ylim([100,140]);

% ylim([100,160]);
legend('',sprintf('%d%s ± %1.1f%s',round(avg_data(end)),'%',round(std_data(end),1),'%'),'location','northeast');
set(gca,'fontsize',fontsize+2,'fontname',fontname)
xlabel('Number of examination codes included','fontsize',fontsize+4,'fontname',fontname);
ylabel('Scanner capacity compared to routine (%)','fontsize',fontsize+4,'fontname',fontname);
if useGlobal
    exportgraphics(gcf,sprintf('%s/Scanner_capacity_fixed_DRB_acceleration_positioning_minimum_%dmin.png',saveLoc,positioning_time_minimum),'Resolution',300);
else
    exportgraphics(gcf,sprintf('%s/Scanner_capacity_sequence_specific_DRB_acceleration_positioning_minimum_%dmin.png',saveLoc,positioning_time_minimum),'Resolution',300);
end