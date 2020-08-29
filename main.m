function [label,time_fall]=main(pathname)

%% find relevant path
find_slash=strfind(pathname,'\');
path_of_folder=pathname(1:find_slash(end-1));

%% load the model
load('model_RUSBOOST.mat');

%% main loop
correct_path=[path_of_folder,'*.csv*'];
fileList = dir(correct_path);
cur_file_name=pathname(find_slash(end-1)+1:end-1);

%% the record signals
% find the file number from file list
file_names={fileList.name};
filenum=find(strcmp(file_names(1,:),cur_file_name));

%% main code
% initialization of the parameters:
filetype='csv';
window=120; %[sec]
overlap=60; %[sec]
amount_par=4; %the amount of parameters
signaltime=60*20; %sec, it is going to change according to the signal's real length
f_A=50; %[Hz] 
f_G=50; %[Hz]
f_P= round(0.99); %[Hz]
f_T= 1; %[sec]
group=[];
feature_mat=[];
labels=[];
count=0;
prev_mat=0;
signal_orig_time=ones(1,signaltime*f_A);
isfirst=1;
A=[];G=[];P=[];T=[];
wind_mat={};
for row=1:window-overlap:length(signal_orig_time)/f_A-window
    feature_vec=[];
    if isfirst==1
       [vec_A,vec_G,vec_P,vec_T,~,A,G,P,T,signal_orig_time] = main_load_data(correct_path,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time);
        isfirst=0;
    else
        [vec_A,vec_G,vec_P,vec_T,~,~,~,~,~,~] = main_load_data(correct_path,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time);    
    end    
    %creating the features for one window:
    wind_mat={vec_A,vec_G,vec_P,vec_T}; % cell because the vectors doesn't have the same length
    if row==1
        prev_mat=wind_mat;
    end
    feature_vec=Extract_Features(wind_mat,prev_mat,f_A,f_G,f_P,f_T);
    count=count+1;
    feature_mat(count,:)=feature_vec;
end
prev_mat=wind_mat;

X=feature_mat;
Y=labels;



%% Normalization
%load the minimum and maximum values
% min_max_path=[path_of_folder,'\min_max_norm.mat'];
load('min_max_norm.mat');

min_x=min_max_norm(1,:);
max_x=min_max_norm(2,:);

X_norm=(X-repmat(min_x,length(X),1))./repmat((max_x-min_x),length(X),1);

%% Features discretization
% bin_num=10;
% X_discrete=X_norm;
% % implementing equal frequency discretization
% for r=1:size(X_norm,2)
%     bin_thresholds=tsprctile(X_norm(:,r),0:bin_num:100);
%     bin_value=bin_thresholds(1:10)+diff(bin_thresholds)/2;
%     for rr=1:10
%         X_discrete(bin_thresholds(rr)<X_norm(:,r) & X_norm(:,r)<=bin_thresholds(rr+1),r)=bin_value(rr);
%     end
% end

%% predicting using the model
load('best_features.mat'); %load the best features index
[Y_pred,score] = predict(model_RUSBOOST, X_norm(:,best_features)); %Predict labels using classification tree

Y_pred(1:2)=0; %there is no fall in the first 2 minutes

% initiation
time_fall=0;
label=0;

falls_vec=find(Y_pred==1);
if ~isempty(falls_vec)
    score_falls=score(falls_vec,2);
    diff_falls_vec=diff(falls_vec);
    if ~isempty(diff_falls_vec)
        adjacent_falls=find(diff_falls_vec==1);
        if ~isempty(adjacent_falls)
            new_adjacent_falls=[adjacent_falls;adjacent_falls(end)+1];
            scores_ind=new_adjacent_falls;
            score_adjacent=score_falls(scores_ind);
            mean_vec=[];
            for i=1:length(scores_ind)-1
                mean_vec(i)=mean(score_adjacent(i:i+1));
            end     
            [~,ind_max]=max(mean_vec);
            local_time=adjacent_falls(ind_max);
            time_fall=falls_vec(local_time);
            label=1;
        end
    end
end
 
end