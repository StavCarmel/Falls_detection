
function model_RUSBOOST=Retrain(pathname)

%% fixing the file names in the folder
new_path=[pathname,'*.c*'];
fixed=fixing_names(pathname,new_path);

%% main loop
correct_path=[pathname,'*.csv*'];
cd =correct_path; %?
fileList = dir(correct_path);

% parameters for the loop:
numberOfFiles = length(fileList);
filetype='csv';
window=120; %[sec]
overlap=60; %[sec]
amount_par=4; %the amount of parameters
signaltime=60*20; %sec, it is going to change according to the signal's real length
f_A=50; %[Hz] 
f_G=50; %[Hz]
f_P= round(0.99); %[Hz]
f_T= 1; %[sec]
continue_cell={};
group=[];
feature_mat=[];
labels=[];
count_group=1;
count=0;
prev_mat=0;
am_feat=10; % the amount of features
wind_mat=[];
for filenum=1:amount_par:numberOfFiles
    signal_orig_time=ones(1,signaltime*f_A);
    filenum
    isfirst=1;
    A=[];G=[];P=[];T=[];
    for row=1:window-overlap:length(signal_orig_time)/f_A-window
        feature_vec=[];
        if isfirst==1
            [if_continue,vec_A,vec_G,vec_P,vec_T,time_A,time_G,time_P,time_T,label,A,G,P,T,signal_orig_time,group_name] = Retrain_load_data(correct_path,filetype,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time);
            isfirst=0;
            if if_continue
                continue_cell=[continue_cell,fileList(filenum).name];
                break
            end
        else
            [~,vec_A,vec_G,vec_P,vec_T,~,~,~,~,label,~,~,~,~,~,group_name] = Retrain_load_data(correct_path,filetype,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time);    
        end              
        %creating the features for one window:
        wind_mat={vec_A,vec_G,vec_P,vec_T}; % cell because the vectors doesn't have the same length
        if filenum==1
            group_prev=group_name;
            count_group=count_group;
        else
            if strcmp(group_name,group_prev)
                count_group=count_group;
            else
                count_group=count_group+1;
            end
        end
        if row==1
            prev_mat=wind_mat;
        end
        feature_vec=Extract_Features(wind_mat,prev_mat,f_A,f_G,f_P,f_T);
        count=count+1;
        feature_mat(count,:)=feature_vec;
        labels(count,1)=label;  
        group=[group;count_group];
        group_prev=group_name;
    end
    prev_mat=wind_mat;
end

X=feature_mat;
Y=labels;

%% Normalization
X_norm=(X-repmat(min(X),length(X),1))./repmat((max(X)-min(X)),length(X),1);
%% save the max and the min values of the best 5 features 
min_max_norm=[min(X);max(X)];
save('min_max_norm','min_max_norm')

% %% Features discretization
% bin_dist=20;
% X_discrete=X_norm;
% % implementing equal frequency discretization
% for r=1:size(X_norm,2)
%     bin_thresholds=tsprctile(X_norm(:,r),0:bin_dist:100);
%     bin_value=bin_thresholds(1:end-1)+diff(bin_thresholds)/2;
%     for rr=1:length(bin_value)
%         X_discrete(bin_thresholds(rr)<X_norm(:,r) & X_norm(:,r)<=bin_thresholds(rr+1),r)=bin_value(rr);
%     end
% end


%% training & Test sets
% person=round(0.7*group(end));
% ind_1=find(group==person);
% %taking the 0.7 precent index in which it changes from one person to another (we assume that all the patients aree listed one after another and not mixed)
% ind_patient=ind_1(end);
% ind_training=1:ind_patient;
% ind_test=(ind_patient+1:size(X_norm,1));
% 
% X_training=X_norm(ind_training,:);
% X_test=X_norm(ind_test,:);
% Y_training=Y(ind_training,1);
% Y_test=Y(ind_test,1);

% %checking the amount of every class, that every group has almost the same from every class
% tabulate(Y)
% tabulate(Y_training)
% tabulate(Y_test)

% after checking the model on train and test set, we want to create a model
% from all the data:
X_training=X_norm;
Y_training=Y;
%%  checking the correlation for every pair of features
feat_feat_cor=corr(X_training,X_training,'Type','Spearman');
feat_feat_cor2=feat_feat_cor-eye(am_feat);
[max_corr,~]=max(max(abs(feat_feat_cor2)))
[ind_2x,ind_2y]=find(feat_feat_cor2>0.8);

%% Select the best feature
feat_label_cor=corr(X_training,Y_training,'Type','Spearman');
[~,best_ind]=max(abs(feat_label_cor));
% [ind, weights]=relieff(X_training,Y_training,3)
% best_ind=ind(1);

%%  Select n best features using SFS
n_feat=5;
% SFS with CFS- Spearman correlation:

max_ind=best_ind;

for n_features=2:n_feat
    fprintf('\nLooking for the best %d features combination\n',n_features)
    CFS=zeros(size(X_training,2)-length(max_ind),1);
    available_features=setdiff((1:size(X_training,2)),max_ind);
    for r=available_features
        CFS(r)=calculate_CFS(feat_label_cor,feat_feat_cor,[max_ind,r]);
        disp(['CFS criterion of feature(s) ',num2str(max_ind),' with feature #',num2str(r),' = ',num2str(CFS(r))])
    end
    [~,ind]=max(CFS);
    max_ind=[max_ind,ind];
    disp(['Max CFS Spearman correlation for feature combination: ',num2str(max_ind)])    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

best_features=max_ind;
save('best_features','best_features');
disp(['The best 5 features are numbers: ',num2str(best_features)])
disp('------------------------------------------')
save('best_features','best_features')



 %% exploring the features
%  gplotmatrix(X_training(:,best_features),[],Y_training)

%% creating model
%% RUSBOOST
num_feat=size(X_training(:,best_features),2);
t=templateTree('MaxNumSplits',num_feat);
model_RUSBOOST=fitcensemble(X_training(:,best_features),Y,'Method','RUSBoost','NumLearningCycles',1500,'Learners',t,'LearnRate',0.03,'nprint',300);

%% the model on the training set
% model_RUSBOOST=fitcensemble(X_training(:,best_features),Y_training,'Method','RUSBoost','NumLearningCycles',1500,'Learners',t,'LearnRate',0.03,'nprint',300);
% 
% % training
% [Y_pred_train,score_train] = predict(model_RUSBOOST, X_training(:,best_features)); %Predict labels using classification tree
% Tree_perf = classperf(Y_training, Y_pred_train ,'Positive',1,'Negative',0); %Evaluate classifier performance
% confusionchart(Y_training,Y_pred_train);
% 
% 
% % test
% [Y_pred_test,score_RUSBOOST] = predict(model_RUSBOOST, X_test(:,best_features)); %Predict labels using classification tree
% Tree_perf_test = classperf(Y_test, Y_pred_test ,'Positive',1,'Negative',0); %Evaluate classifier performance
% confusionchart(Y_test,Y_pred_test);
% f1_score_test=(2*Tree_perf_test.Sensitivity*Tree_perf_test.PositivePredictiveValue)/(Tree_perf_test.Sensitivity+Tree_perf_test.PositivePredictiveValue)
% [x_ROC,y_ROC,threshold,AUC1]=perfcurve(Y_test,score_RUSBOOST(:,2),1);
% [x_PRC,y_PRC,threshold,AUC2]=perfcurve(Y_test,score_RUSBOOST(:,2),1,'XCrit','tpr','YCrit','ppv');
% 
% Sensitivity=Tree_perf_test.Sensitivity;
% PPV=Tree_perf_test.PositivePredictiveValue;
% Specificity=Tree_perf_test.Specificity;
% 
% f1_score=(2*x_PRC.*y_PRC)./(x_PRC+y_PRC);
% max(f1_score)
% [~,max_f]=max(f1_score)
% New_Sensitivity=x_PRC(max_f);
% new_PPV=y_PRC(max_f);
% New_Specifity=1-x_ROC(max_f);
% 
% figure;plot(x_ROC,y_ROC);line([0 1],[0 1],'color','r');xlim([0 1]);ylim([0 1]);title('ROC curve');xlabel('FPR'),ylabel('TPR');
% figure;plot(x_PRC,y_PRC);yline(1/10,'r');title('PRC curve');xlabel('Recall'),ylabel('Presicion');
% 
% figure
% instance=1:length(Y_test);
% ind=find(Y_test==1);
% plot(instance(ind),score_RUSBOOST(ind,2),'bx')
% hold all
% ind=find(Y_test==0);
% plot(instance(ind),score_RUSBOOST(ind,2),'g*')
% ind=find(Y_test~=Y_pred_test);
% plot(instance(ind),score_RUSBOOST(ind,2),'ro','MarkerSize',8)
% line([instance(1),instance(end)],yline(2.26),'color',[1 0.75 0])
% legend({'Fall','No fall','Misclassified','Threshold'})
% ylabel('Score')
% xlabel('Instance')
% ylim([0 4])
% 
% thershold_for_deterioration=threshold(max_f);
% Ymodel=double(score_RUSBOOST(:,2)>thershold_for_deterioration);
% figure
% ind=find(Y_test==1);
% plot(instance(ind),score_RUSBOOST(ind,2),'bx')
% hold all
% ind=find(Y_test==0);
% plot(instance(ind),score_RUSBOOST(ind,2),'g*')
% ind=find(Ymodel~=Y_test);
% plot(instance(ind),score_RUSBOOST(ind,2),'ro','MarkerSize',8)
% line([instance(1),instance(end)],[thershold_for_deterioration thershold_for_deterioration],'color',[1 0.75 0])
% legend({'Fall','No fall','Misclassified','Threshold'})
% ylabel('Score')
% xlabel('Instance')



% confusionmat(Ytest,Ymodel)
% %% SVM
% c=25; %the cost;
% cost_mat=[0,1;c,0];
% model_SVM=fitcsvm(X_training,Y_training,'Cost',cost_mat);
% [Y_pred_SVM,score_SVM] = predict(model_SVM, X_test(:,best_features));%Predict labels using classification tree
% Tree_perf_SVM = classperf(Y_test, Y_pred_SVM ,'Positive',1,'Negative',0);
% confusionchart(Y_test,Y_pred_SVM)
% 
% [x_ROC,y_ROC,threshold,AUC]=perfcurve(Y_test,score_SVM(:,2),1);
% [x_PRC,y_PRC,threshold,AUC]=perfcurve(Y_test,score_SVM(:,2),1,'XCrit','tpr','YCrit','ppv');
% 
% 
% figure;plot(x_ROC,y_ROC);line([0 1],[0 1],'color','r');xlim([0 1]);ylim([0 1])
% figure;plot(x_PRC,y_PRC);line([0 1],[1 0],'color','r');xlim([0 1]);ylim([0 1])
% 
% %% Bagging
% model_BAG = fitcensemble(X_training,Y_training,'Method','bag', ...
% 'NumLearningCycles',1000,'Learners',t,'nprint',500,'Cost',[0 1;15 0]);
% 
% % model_BAG = fitcensemble(X_training,Y_training,'Method','bag', ...
% % 'NumLearningCycles',10,'Learners',t,'nprint',500);
% 
% [Y_pred_BAG,score_BAG] = predict(model_BAG, X_test(:,best_features));%Predict labels using classification tree
% Tree_perf_BAG = classperf(Y_test, Y_pred_BAG ,'Positive',1,'Negative',0);
% confusionchart(Y_test,Y_pred_BAG)
% 
% [x_ROC,y_ROC,threshold,AUC]=perfcurve(Y_test,score_BAG(:,2),1);
% [x_PRC,y_PRC,threshold,AUC]=perfcurve(Y_test,score_BAG(:,2),1,'XCrit','tpr','YCrit','ppv');
%  f1_score=(2*x_PRC.*y_PRC)./(x_PRC+y_PRC);
% 
% figure;plot(x_ROC,y_ROC);line([0 1],[0 1],'color','r');xlim([0 1]);ylim([0 1])
% figure;plot(x_PRC,y_PRC);line([0 1],[1 0],'color','r');xlim([0 1]);ylim([0 1])

end


