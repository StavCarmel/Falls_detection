function [if_continue,vec_A,vec_G,vec_P,vec_T,time_A,time_G,time_P,time_T,label,A,G,P,T,signal_orig_time,group_name] = Retrain_load_data(folder,filetype,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time)
% num 0/1 -isfirst     : means that if it is the first window of this session than don't load again and use the tables from before. less running time
% num 0/1 -if_continue : if there is somethnig wrong with the recordnig
% files then the code will pass those recording and continue to load the
% next batch of files

%% initialization
if_continue=0;
% cd =folder;
% pathname = folder;
fileList = dir(folder);
numberOfFiles = length(fileList);

path_of_folder=folder(1:end-6);
% calculating how many indicies we need to take from every signal in order
% to get one window
ind_A=f_A*window;
ind_G=f_G*window;
ind_P=f_P*window;
ind_T=window/f_T;

%% associating the file names to its corresponding record
i=filenum;
fileName1 = fileList(i).name;
fileName2 = fileList(i+1).name;
fileName3 = fileList(i+2).name;
fileName4 = fileList(i+3).name;

files={fileName1,fileName2,fileName3,fileName4};
for file=1:4
    if contains(files{file},'acc') 
       fileA=files{file};
    elseif contains(files{file},'gyr') 
        fileG=files{file};
    elseif contains(files{file},'bar') 
        fileP=files{file};
    else
        fileT=files{file};
    end
end
    % labeling

str_time_fall=fileName1(12:13);
if ~strcmp(str_time_fall,'00')
    time_fall=str2num(str_time_fall);
    sec_time_fall=time_fall;
    % looking at minutes
    range_falling=[sec_time_fall,sec_time_fall+((window-overlap)/60)];
    left=(row-1)/60; %limits of the window, row -1 because it starts from row=1 row= 61
    right=(row-1)/60+2;
    if left<=range_falling(1)&& right>=range_falling(2)
       label=1;
    else
        label=0;
    end
else
    label=0; %if the record doesn't include fall
end 

% the sensor name of this gruop in order to separate to train and test later on
group_name=fileA(1:4);

%% loading the files as tables
if isfirst==1 % it means that this is the first window for this session than load the files and create a table
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames') 
    A=readtable(strcat(path_of_folder,fileA),'ReadVariableNames',true);  
    G=readtable(strcat(path_of_folder,fileG),'ReadVariableNames',true);
    P=readtable(strcat(path_of_folder,fileP),'ReadVariableNames',true); 
    T=readtable(strcat(path_of_folder,fileT),'ReadVariableNames',true); 

    %% checking that there are no missing titles (students' data collection problems)
    A_names=A.Properties.VariableNames;
    G_names=G.Properties.VariableNames;
    P_names=P.Properties.VariableNames;
    T_names=T.Properties.VariableNames;

    if sum(contains(A_names,'elapsed_s_'))+sum(contains(G_names,'elapsed_s_'))+sum(contains(P_names,'elapsed_s_'))+sum(contains(T_names,'elapsed_s_'))<4
        if_continue=1;
    end

    if sum(contains(A_names,'temperature_C_'))+ sum(contains(G_names,'temperature_C_'))+sum(contains(P_names,'temperature_C_'))+sum(contains(T_names,'temperature_C_'))==0 
        if_continue=1;
    end



    %% checking that all the recording has at least 90% of the 20 minutes recording
    if ~if_continue
        if length(A.elapsed_s_)<0.90*50*20*60 || length(G.elapsed_s_)<0.90*50*20*60 ||length(T.elapsed_s_)<0.90*20*60 ||length(P.elapsed_s_)<0.90*0.99*20*60 
            warning(['One or more of the recordings related to ',fileA, ' aren''t long enough'])
            if_continue=1;

        end
    end

    %% checking the consecutive missing samples in the records
    if ~if_continue
        time_A=A.elapsed_s_;
        time_G=G.elapsed_s_;
        time_P=P.elapsed_s_;
        time_T=T.elapsed_s_;

        diff_vec_A=diff(time_A);
        diff_vec_G=diff(time_G);
        diff_vec_P=diff(time_P);
        diff_vec_T=diff(time_T);

        gap1=50*(1/f_A);
        gap2=50*(1/f_G);
        gap3=5*(1/f_P);
        gap4=5*(1/f_T);

        if sum(diff_vec_A>=gap1)>0
            warning(['There are more than 50 consecutive missing samples in file ',fileA])
            if_continue=1;
        end

        if sum(diff_vec_G>=gap2)>0
            warning(['There are more than 50 consecutive missing samples in file ',fileG])
            if_continue=1;
        end

        if sum(diff_vec_P>=gap3)>0
            warning(['There are more than 5 consecutive missing samples in file ',fileP])
            if_continue=1;
        end

        if sum(diff_vec_T>=gap4)>0
            warning(['There are more than 5 consecutive missing samples in file ',fileT])
            if_continue=1;
        end
        
        if ~if_continue
            signal_orig_time=(A.elapsed_s_); % the amount of samples in the acc vector
        end
    end

    %% setting empty values if there was a problem and the code will pass this specific value, to avoid error 
    if if_continue
        vec_A=[];
        vec_G=[];
        vec_P=[];
        vec_T=[];
        time_A=[];
        time_G=[];
        time_P=[];
        time_T=[];
        label=[];
        group_name=[];
        return
    end

    %% checking if resampling is needed

    need_samp_A=sum(diff_vec_A>1/f_A);
    if need_samp_A>0
       [a,time_A]=resample(([A.x_axis_g_,A.y_axis_g_,A.z_axis_g_]),time_A,f_A,3,1);
       A_mat=[time_A,a];
       A=array2table(A_mat,'VariableNames',{'elapsed_s_','x_axis_g_','y_axis_g_','z_axis_g_'});
    end

    need_samp_G=sum(diff_vec_G>1/f_G);
    if need_samp_G>0
       [g,time_G]=resample(([G.x_axis_deg_s_,G.y_axis_deg_s_,G.z_axis_deg_s_]),time_G,f_G,3,1);
       G_mat=[time_G,g];
       G=array2table(G_mat,'VariableNames',{'elapsed_s_','x_axis_deg_s_','y_axis_deg_s_','z_axis_deg_s_'});
    end

    need_samp_P=sum(diff_vec_P>1/f_P);
    if need_samp_P>0
       [p,time_P]=resample((P.pressure_Pa_),time_P,f_P,3,1);
       P_mat=[time_P,p];
       P=array2table(P_mat,'VariableNames',{'elapsed_s_','pressure_Pa_'});
    end

    need_samp_T=sum(diff_vec_T>1/f_T);
    if need_samp_T>0
       [t,time_T]=resample((T.temperature_C_),time_T,f_T,3,1);
       T_mat=[time_T,t];
       T=array2table(T_mat,'VariableNames',{'elapsed_s_','temperature_C_'});
    end

end

%% taking the window samples and creating the vectors
A_window=A((row-1)*f_A+1:(row-1)*f_A+ind_A,:);
G_window=G((row-1)*f_G+1:(row-1)*f_G+ind_G,:);   
P_window=P((row-1)*f_P+1:(row-1)*f_P+ind_P,:);   
T_window=T((row-1)*f_T+1:(row-1)*f_T+ind_T,:);

vec_A=mean([A_window.x_axis_g_,A_window.y_axis_g_,A_window.z_axis_g_],2);
vec_G=mean([G_window.x_axis_deg_s_,G_window.y_axis_deg_s_,G_window.z_axis_deg_s_],2);
vec_P=P_window.pressure_Pa_;
vec_T=T_window.temperature_C_;

%% we need the times only when isfirst=1; so put zero to avoid error
if ~isfirst
time_A=0;
time_G=0;
time_P=0;
time_T=0;
end
end

