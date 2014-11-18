function [features] = extract_features(data,fs,options)
%EXTRACT_FEATURES Summary of this function goes here
%   Detailed explanation goes here

data = clean_data(data);
features = [];
[r,~] = size(data);

freq_bands = [0.5,4;4,7;7,13;13,15;14,30;30,45;65,120];
freq_bands = freq_bands/200;
band_names = {'delta','theta','alpha','low_beta','high_beta','beta','alpha'};


%%%%%%%%%%%%
%%%%%%%%%%%%

% % % % if options.univariate %%% univariate features
% % % %     for channel = 1:r
% % % %         x = data(channel,:); 
% % % %         if options.power%% power of the signal
% % % %             features.(['ch',int2str(channel)]).power = sum(x.^2);
% % % %         end
% % % %         if options.psd %% power spectral density...
% % % %             [fd,~] = get_fft(x,fs);
% % % %             PSD = (sum(fd.^2))/length(fd);
% % % %             features.(['ch',int2str(channel)]).psd = PSD;       
% % % %         end
% % % %     end
% % % % end
    
if options.bivariate %%% bivariate features
%     n = r*(r-1)/2;
%     cc = zeros(1,n);
%     cc_lag = zeros(1,n);
%     ni = zeros(1,n);
%     de = zeros(1,n);
%     for q = 1:length(band_names)
%        features.(band_names{q}).coh = zeros(1,n);
%        features.(band_names{q}).ed = zeros(1,n);
%        features.(band_names{q}).coh = zeros(1,n);
%     end
%     ch = 0;
f=1;
    for channel1 = 2:4 %%% iterate through each row,column pair
    for channel2 = 1:(channel1-1) %%% and get signals for comparison
    
        
        fprintf('%d%d\n',channel1,channel2);
        x_a = data(channel1,:);
        x_b = data(channel2,:);              
        
        if options.cc %%% get max cross_correlation and lag of max cross corellation
            [acorr,l] = xcorr(x_a,x_b);
            [max_cross_corellation,x] = max(acorr);
            lag = l(x);
            features(f) = max_cross_corellation;
            f=f+1;
            features(f) = lag;
            f=f+1;
        end
%         fprintf('cc')            
        if options.ni || options.de 
        %%% nonlinear interdependence or dynamic entrainment ... this way 
        %%% you just do time delay embedding once... saves time...

        % get time delay embedded matrices X_a and X_b
        T = 10;
        d = 10;
        K = 5;
        max_offset = T*(d-1);
        N = length(x_a)-max_offset;
        X_a = zeros(d,N);
        X_b = zeros(d,N);

        for ind = 0:(d-1)
            offset = ind*T;
            X_a((ind+1),:) = x_a((offset+1):(end-max_offset+offset));
            X_b((ind+1),:) = x_b((offset+1):(end-max_offset+offset));
        end
        [~,N] = size(X_a);
        if options.ni
            %%find K nearest neighbor of each element of X_a , X_b in state space
            distance_a_a = ones(N,N) ;
            distance_b_b = zeros(N,N) ;

            for ind_1 = 1:N
            for ind_2 = 1:N
                if ind_1 == ind_2
                    distance_a_a(ind_1,ind_2) = 999;
                    distance_b_b(ind_1,ind_2) = 999;
                else
                    distance_a_a(ind_1,ind_2) = sqrt(sum((X_a(:,ind_1)-X_a(:,ind_2)).^2));
                    distance_b_b(ind_1,ind_2) = sqrt(sum((X_b(:,ind_1)-X_b(:,ind_2)).^2));
                end
            end
            end
            

            for ind = 1:d
                [~,min_a_a] = sort(distance_a_a);
                [~,min_b_b] = sort(distance_b_b);
            end

            min_k_a_a = min_a_a(1:K,:); % get indices of the min k
            min_k_b_b = min_b_b(1:K,:); % distances 

% % % % % %             min_k_distance_a_a = distance_a_a(min_k_a_a);
% % % % % %             sum(min_k_distance_a_a);
% % % % % %             min_k_distance_b_b = distance_b_b(min_k_b_b);
% % % % % %             sum(min_k_distance_b_b);
            %%% sum distances to k nearest neighbors in state space...
            R_a_a = sum(distance_a_a(min_k_a_a))/K;
            R_a_b = sum(distance_a_a(min_k_b_b))/K;
            R_b_b = sum(distance_b_b(min_k_b_b))/K;
            R_b_a = sum(distance_b_b(min_k_a_a))/K;
            
            S_a_b = sum(R_a_a./R_a_b)/N;
            S_b_a = sum(R_b_b./R_b_a)/N;
            non_linear_interdependence = (S_a_b + S_b_a)/2;
%             fprintf('%d\n',non_linear_interdependence)
            features(f) = non_linear_interdependence;
            f=f+1;
%             fprintf('ni')
        end   
            
        if options.de %%% dynamical entrainment

            dt = 20;
            
%             t_j_a = min_a_a(1,1:(N-dt)); %
%             t = 1:(N-dt);
%             max(t_j_a+dt)
%             lambda_X_a_dt = X_a(:,(dt+1):N)-X_a(:,(t_j_a+dt));
%             lambda_X_a = X_a(:,t) - X_a(:,t_j_a);
%             t_j_b = min_b_b(1,1:N); %
%             lambda_X_b_dt = X_b(:,(dt+1):N-1)-X_b(:,(t_j_b+dt));
%             lambda_X_b = X_b(:,t) - X_b(:,t_j_b);

            t_j_a = 2:(N-dt); %
            t = 1:(N-dt-1);
            lambda_X_a_dt = X_a(:,t+dt)-X_a(:,(t_j_a+dt));
            lambda_X_a = X_a(:,t) - X_a(:,t_j_a);
            t_j_b = 2:(N-dt); %
            lambda_X_b_dt = X_b(:,t+dt)-X_b(:,(t_j_b+dt));
            lambda_X_b = X_b(:,t) - X_b(:,t_j_b);


            STL_max_a = sum(    log2(  sqrt( sum( lambda_X_a_dt .^2 ))/ sqrt( sum( lambda_X_a .^2 )) )  ) / ((N-dt)*dt);
            STL_max_b = sum(    log2(  sqrt( sum( lambda_X_b_dt .^2 ))/ sqrt( sum( lambda_X_b .^2 )) )  ) / ((N-dt)*dt);

            dynamical_entrainment = abs(STL_max_a-STL_max_b);
            features(f) = dynamical_entrainment;
            f=f+1;
%             fprintf('de')
        end
        end
        
        %%%%%
        %%%%% Wavelet-based synchrony measures
        %%%%%
        %%%%% 
        
    
            
        if options.pls || options.ed || options.coh
        %%% phase locking synchrony of i,j, entropy and overall synchrony...    
        %%% get phi_a and phi_b using hilbert/ wavelet
        %%% transform, or whatever... separate into signals
        %%% for eeg frequency bands....
        order = 1;
        
        for band = 1:length(freq_bands(:,1))
            
            lowFreq = freq_bands(band,1);
            hiFreq = freq_bands(band,2);
            [b,a] = butter(order,[lowFreq,hiFreq]/(fs/2),'bandpass');
            y_a = filter(b,a,x_a);
            y_b = filter(b,a,x_b);

            phi_a = unwrap(angle(hilbert(y_a)));
            phi_b = unwrap(angle(hilbert(y_b)));
            phi_difference = phi_a - phi_b;

            %%%
            if options.pls %%% phase locking synchrony....
            y = sum(exp((phi_difference)*1i))/length(phi_a);
            phase_locking_synchrony = sqrt(sum((y.*conj(y)).^2));
            features(f) = phase_locking_synchrony;
            f=f+1;
            end
            
            if options.ed%%% entropy difference
            M = 10;
            pm = zeros(1,M);
            last = length(phi_difference);
            phi_min = min(phi_difference);
            phi_max = max(phi_difference);
            for ind = 1:M
                qwer = ((1/M)*(phi_max-phi_min)*ind+phi_min);
                abcdef = sum( phi_difference > qwer);
                pm = (last - abcdef)/length(phi_difference);
                if pm == 0
                    pm = 0.0000000000000000001;
                end
                last = abcdef;
            end
            entropy_difference = ((log(M)+ sum(pm.*log(pm))))/log(M);
            features(f) = entropy_difference;
            f=f+1;
            end
            
            if options.coh %% coherence
            coherence = mscohere(phi_a,phi_b);
            features(f) = sum(coherence);
            f=f+1;
            end
            
        end
        end
    end
    end
%     features.cc = cc;
%     features.ni = ni;
%     features.de = de;
end
% [~,c]=size(features);
% for q= 1:c
%     current = features(:,q);
%     current = current-mean(current);
%     current = current/std(current);
%     features(:,q) = current(:,1);
% end
end
