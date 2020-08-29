function [vec_A,vec_G,vec_P,vec_T,label,A,G,P,T,signal_orig_time] = main_load_data(correct_path,window,filenum,row,f_A,f_G,f_P,f_T,isfirst,A,G,P,T,overlap,signal_orig_time)
% num 0/1 -isfirst     : means that if it is the first window of this session than don't load again and use the tables from before. less running time
% num 0/1 -if_continue : if there is somethnig wrong with the recordnig
% files then the code will pass those recording and continue to load the
% next batch of files

%% initialization
fileList = dir(correct_path);
numberOfFiles = length(fileList);
path_of_folder=correct_path(1:end-6);

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

    %% checking if resampling is needed
    time_A=A.elapsed_s_;
    time_G=G.elapsed_s_;
    time_P=P.elapsed_s_;
    time_T=T.elapsed_s_;

    diff_vec_A=diff(time_A);
    diff_vec_G=diff(time_G);
    diff_vec_P=diff(time_P);
    diff_vec_T=diff(time_T);
    
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
    
    signal_orig_time=(A.elapsed_s_);% the amount of samples in the acc vector- after resampling
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

end

