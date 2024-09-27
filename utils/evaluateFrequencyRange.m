function [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel]= evaluateFrequencyRange(powerspectra_image,alphabeta,gamma,minrange)


set_pval=0.05;
individual_goodness_result = nan(3,1);
nanboolean=~isnan(powerspectra_image(:,1));
startrow= find(nanboolean ==1,1);
endrow = startrow+sum(nanboolean)-1;

for superficialchannel = startrow: endrow-minrange
    for deepchannel = superficialchannel+minrange:endrow
        powspec_window = powerspectra_image(superficialchannel:deepchannel,:);

        maxpow =squeeze(max(powspec_window,[],1));
        relpow = squeeze(powspec_window)./maxpow;
        lowband = mean(relpow(:,alphabeta),2);
        highband = mean(relpow(:,gamma),2);



        %extra normalization toggle
        %lowband = lowband-min(lowband);
        %lowband = lowband/max(lowband);
        %highband = highband-min(highband);
        %highband = highband/max(highband);

        [lowband_b,bint,r,rint,stats_low] = regress(lowband,horzcat([1:length(lowband)]',ones(length(lowband),1)));
        lowband_rsquared=stats_low(1);
        lowband_pval=stats_low(3);

        [highband_b,bint,r,rint,stats_high] = regress(highband,horzcat([1:length(highband)]',ones(length(highband),1)));
        highband_rsquared=stats_high(1);
        highband_pval=stats_high(3);
        reward = 0.04*(deepchannel-superficialchannel+1)+0.72;

        beta_peak_locations = find((lowband ==max(lowband)));
        gamma_peak_locations = find((highband ==max(highband)));

        if mean(beta_peak_locations)> length(lowband)/2
            beta_peak_index = max(find(lowband ==max(lowband)));
        else
            beta_peak_index = min(find(lowband ==max(lowband)));
        end

        if mean(gamma_peak_locations)> length(highband)/2
            gamma_peak_index= max(find(highband ==max(highband)));
        else
            gamma_peak_index= min(find(highband ==max(highband)));
        end

        if beta_peak_index==1
            beta_peak_max_check = [superficialchannel==startrow & lowband(beta_peak_index)>lowband(beta_peak_index+1)];
        elseif beta_peak_index==length(lowband)
            beta_peak_max_check =[deepchannel ==endrow & lowband(beta_peak_index)>lowband(beta_peak_index-1)];
        else
            beta_peak_max_check = [lowband(beta_peak_index)>lowband(beta_peak_index+1) & lowband(beta_peak_index)>lowband(beta_peak_index-1)];
        end

        if gamma_peak_index==1
            gamma_peak_max_check = [superficialchannel==startrow & highband(gamma_peak_index)>highband(gamma_peak_index+1)];
        elseif gamma_peak_index==length(highband)
            gamma_peak_max_check =[deepchannel ==endrow & highband(gamma_peak_index)>highband(gamma_peak_index-1)];
        else
            gamma_peak_max_check = [highband(gamma_peak_index)>highband(gamma_peak_index+1) & highband(gamma_peak_index)>highband(gamma_peak_index-1)];
        end

        if gamma_peak_max_check & beta_peak_max_check
        else
            goodness=nan;
            individual_goodness_result = horzcat(individual_goodness_result, [goodness; superficialchannel;deepchannel]);
            continue
        end


        if lowband_pval<set_pval & highband_pval<set_pval
            if lowband_b(1)>0 & highband_b(1)<0
                goodness = highband_rsquared*lowband_rsquared*reward;
            elseif lowband_b(1)<0 & highband_b(1)>0
                goodness = -highband_rsquared*lowband_rsquared*reward;
            else
                goodness=nan;
                continue
            end
            individual_goodness_result = horzcat(individual_goodness_result, [goodness; superficialchannel;deepchannel]);
        else
            continue
        end
    end
end


a = max((individual_goodness_result(1,:)));
b = min((individual_goodness_result(1,:)));
if a>=abs(b)

    goodness_k_value = a ;
    superficialchannel = individual_goodness_result(2,find(individual_goodness_result(1,:) ==a,1));
    deepchannel=individual_goodness_result(3,find(individual_goodness_result(1,:) ==a,1));
    upright_boolean=1;



else
    goodness_k_value = b ;
    superficialchannel = individual_goodness_result(2,find(individual_goodness_result(1,:) ==b,1));
    deepchannel=individual_goodness_result(3,find(individual_goodness_result(1,:) ==b,1));
    upright_boolean=0;
end
powspec_window = squeeze(powerspectra_image(superficialchannel:deepchannel,:));

maxpow =squeeze(max(powspec_window,[],1));
relpow = squeeze(powspec_window)./maxpow;
lowband = mean(relpow(:,alphabeta),2);
highband = mean(relpow(:,gamma),2);


beta_peak_locations = find(lowband ==max(lowband));
gamma_peak_locations = find(highband ==max(highband));

if mean(beta_peak_locations)> length(lowband)/2
    beta_peak_index = max(find(lowband ==max(lowband)));

else
    beta_peak_index = min(find(lowband ==max(lowband)));
end
if mean(gamma_peak_locations)> length(highband)/2
    gamma_peak_index= max(find(highband ==max(highband)));
else

    gamma_peak_index= min(find(highband ==max(highband)));
end

betamaxchannel = superficialchannel+beta_peak_index-1;
gammamaxchannel = superficialchannel+gamma_peak_index-1;


band_diff=highband-lowband;
band_diff_inverted = lowband-highband;
crossoverchannels=[];
for channel_index = 1:length(lowband)-2
    if goodness_k_value>0

        %gamma is greater than beta, and then it crosses over to
        %beta being greater than gamma
        if highband(channel_index)>lowband(channel_index) && lowband(channel_index+1)>highband(channel_index+1)
            if  abs(band_diff(channel_index))<= abs(band_diff(channel_index+1))
                crossoverchannels= [crossoverchannels channel_index];
            else
                crossoverchannels=[crossoverchannels channel_index+1];
            end


        elseif highband(channel_index)>lowband(channel_index) &&highband(channel_index+1)==lowband(channel_index+1) && lowband(channel_index+2)>highband(channel_index+2)
            crossoverchannels=[crossoverchannels channel_index+1];
        end


    elseif goodness_k_value<0
        if highband(channel_index)<lowband(channel_index) && lowband(channel_index+1)<highband(channel_index+1)
            if  abs(band_diff(channel_index))<= abs(band_diff(channel_index+1))
                crossoverchannels=[crossoverchannels channel_index];
            else
                crossoverchannels=[crossoverchannels channel_index+1];
            end


        elseif highband(channel_index)<lowband(channel_index) &&highband(channel_index+1)==lowband(channel_index+1) && lowband(channel_index+2)<highband(channel_index+2)
            crossoverchannels=[crossoverchannels channel_index+1];
        end
    end
end
%now, all valid cross over points have been collected

number_of_crosses = length(crossoverchannels);
crossover_ratings = horzcat(crossoverchannels',nan(number_of_crosses,1));
if number_of_crosses ==0
    final_crossoverchannel =nan;
elseif number_of_crosses==1
    final_crossoverchannel =crossoverchannels+superficialchannel-1;
else
    if goodness_k_value>0

        for index = 1: number_of_crosses
            crossover_choice = crossoverchannels(index);
            % calculates difference between highband and lowband on
            % each side of the crossover
            crossover_ratings(index,2) = sum(band_diff(1:crossover_choice))-sum(band_diff(crossover_choice:length(lowband)));
        end


    elseif goodness_k_value<0
        for index = 1: number_of_crosses
            crossover_choice = crossoverchannels(index);
            crossover_ratings(index,2) = sum(band_diff_inverted(1:crossover_choice))-sum(band_diff_inverted(crossover_choice:length(lowband)));
        end
    end
    %return index of ideal crossover point
    [~,indexval] = max(crossover_ratings(:,2));

    final_crossoverchannel = crossover_ratings(indexval,1)+superficialchannel-1;
end
end














