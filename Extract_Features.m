function feature_vec = Extract_Features(wind_cell,prev_cell,fA,fG,fP,fT)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% wind_mat is a matrix containing the vectors of every signal (Acc, Gyr,Bar, Tmp)
% the cells dim is n X 4 , n is the number of samples in the window (rows in the matrix) and
% 4  is the number of signals (columns).
% building the matrix:
% ACC_vec - sqrt(x^2+y^2+z^2)
% Gyr_vec -sqrt(x^2+y^2+z^2)
% Bar_Vec- the pressure column from the excell file
% Tmp_Vec- the temperature column from the excell file
% put every signal by this order in the columns of the cell;

%fA is the sampling frequency of Acc- 50
%fG is the sampling frequency of Gyr- 50
%fP is the sampling frequency of Bar- round(0.99)
%fT is the sampling frequency of Tmp- 1

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

feature_vec=zeros(1,10);

A_sig=wind_cell{:,1};
G_sig=wind_cell{:,2};
P_sig=wind_cell{:,3};
T_sig=wind_cell{:,4};

A_prev=prev_cell{:,1};
G_prev=prev_cell{:,2};
P_prev=prev_cell{:,3};
T_prev=prev_cell{:,4};

[P_P,f_P]=pwelch(P_sig-mean(P_sig),[],[],[],fP);
[P_G,f_G]=pwelch(G_sig-mean(G_sig),[],[],[],fG);
[P_A,f_A]=pwelch(A_sig-mean(A_sig),[],[],[],fA);
[P_Aprev,f_Aprev]=pwelch(A_prev-mean(A_prev),[],[],[],fA);

%% exploring the signal
% figure();
% plot([1:length(A_sig)],A_sig)
% hold on
% plot([1:length(A_prev)],A_prev)
% title('Acc time')
% legend('curr','prev')
% 
% figure();
% plot(f_A,log10(P_A))
% hold on
% plot(f_Aprev,log10(P_Aprev))
% title('Acc frequency')
% legend('curr','prev')


%% 1. The difference between the PtP of the Acc and the std of the rest of the window
[max_A,idx_max_A]=max(A_sig);
[min_A,idx_min_A]=min(A_sig);

if idx_max_A>idx_min_A
    Right=idx_max_A;
    Left=idx_min_A;
else
    Right=idx_min_A;
    Left=idx_max_A;
end
    
PtP_A=max_A-min_A;
vec_std=[A_sig(1:Left);A_sig(Right:end)];
std_restof_A=std(vec_std);
diff_PtP_std_A=abs(PtP_A-std_restof_A);

feature_vec(1,1)=diff_PtP_std_A;

%% 2. the ratio between the (max_std-min_std) and the time between the windows
wind=10*fA; %[sec]
max_std=0;
idx_max=0;
min_std=inf;
idx_min=0;
for idx=1:wind:length(A_sig)-wind
    cur_wind=A_sig(idx:idx+wind);
    cur_std=std(cur_wind);
    if cur_std>max_std
        max_std=cur_std;
        idx_max=wind+idx;
    elseif cur_std<min_std
        min_std=cur_std;
        idx_min=wind+idx;
    end
end

if idx_max>idx_min
    Right=idx_max;
    Left=idx_min;
else
    Right=idx_min;
    Left=idx_max;
end

dist_max_min_std=(max_std-min_std)/(Right-Left);
feature_vec(1,2)=dist_max_min_std;

%% 3 correlation with template

load('template_values.mat')
norm_A=(A_sig-min(A_sig))./(max(A_sig)-min(A_sig));
norm_template=(template_values-min(template_values))./(max(template_values)-min(template_values));
[match,lags] = xcov(norm_A,norm_template);
best_match=max(match);
feature_vec(1,3)=best_match;

%% 4 number of local extrema 
feature_vec(1,4)=sum(diff(sign(diff(A_sig)))~=0);

%% 5 areau under curve of Acc
feature_vec(1,5)=trapz(A_sig);

%% 6 difference in std of Acc compared with previous section
feature_vec(1,6)=abs(std(A_sig)-std(A_prev)); 

%% 7 difference between areau under curve of power of Acc between frequencies 10 and 25 (there we saw a difference between fall and no fall)
feature_vec(1,7)=trapz(P_A(10<f_A & f_A<25))-trapz(P_Aprev(10<f_Aprev & f_Aprev<25));

% 8 number of values from current Acc window that are higher than the previous window's mean
feature_vec(1,8)=sum(abs(A_sig)>abs(max(A_prev)));

% 9 difference between the max of power of Acc between frequencies 16 and 18 (there we saw a big peak both in windows containing fall and no fall)
feature_vec(1,9)=max(P_A(16<f_A & f_A<18))-max(P_Aprev(16<f_Aprev & f_Aprev<18));

% 10 correlation between the flipped Acc sig and the curr Acc sig, only if the std is high 
if std(A_sig)>0.1
    [r,~] = corr(A_prev,fliplr(A_sig));
feature_vec(1,10)=r;
else
    feature_vec(1,10)=0;
end

% features with high correlation to others
% max_diff_A=max(diff(P_A));
% feature_vec(1,3)=max_diff_A;
% feature_vec(1,6)=max(A_sig)-min(A_sig);  % dynamic range capped by 5e4
% feature_vec(1,10)=median(P_A(10<f_A & f_A<25));

end

