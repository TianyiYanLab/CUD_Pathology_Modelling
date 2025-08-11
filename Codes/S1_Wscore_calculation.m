%Load subject lists of CUD
CUD_list='D:\OneDrive\1_Preprocessing\SUDMEX_CONN\SUDMexCONN_Case.txt';
fileID = fopen(CUD_list, 'r');
sub_list = textscan(fileID, '%s');  
fclose(fileID); 
CUD = sub_list{1}; 

%Load Non-SUD subject lists from the HCP-YA and SUDMEX-CONN
HCP_list='D:\OneDrive\2_subs_filtering\HCP_subject_list_NonSUD_v3_302.txt';
HCP = load(HCP_list); 

HC_list='D:\OneDrive\1_Preprocessing\SUDMEX_CONN\SUDMexCONN_HC.txt';
fileID2 = fopen(HC_list, 'r');
sub_list2 = textscan(fileID2, '%s');  
fclose(fileID2); 
HC = sub_list2{1}; 


%Prepare the covariates data for Cases
csv ='D:\fMRI_Data\SUDMEX_CONN\5123331\Connectome_demographics_050219.xlsx';
x_CUD = extract_subject_info(CUD_list,csv,{'age','sex'},'rmNaN');

%Prepare the covariates data for HCs
[sex,age] = CBIG_parse_delimited_txtfile('D:\OneDrive\Data_docs\HCP_docs\restricted_S1200.csv', {'Gender'}, ...
{'Age_in_Yrs'}, 'Subject', CBIG_text2cell(HCP_list), ',');
sex = strcmp(sex, 'F')+1; % Male:1,Female:2
x_HCP = [sex,age];

csv_HC = 'D:\fMRI_Data\SUDMEX_CONN\5123331\Connectome_demographics_050219.xlsx';
x_HC = extract_subject_info(HC_list,csv_HC,{'age','sex'},'rmNaN');
x_HC=[x_HC;x_HCP];


%Prepare the cortical thickness data for Cases 
tmp=load('D:\OneDrive\1_Preprocessing\thickness136.mat');
thickness136 = tmp.thickness;
y_CUD = thickness136(1:length(CUD), :); %corresponding to the CUD_list!

%Prepare the cortical thickness data for HCP'HCs
tmp=load('D:\OneDrive\5_FC_NBS_HAU\thickness416.mat');
y_HCP = tmp.thickness(end-length(HCP)+1:end, :);

%Prepare the cortical thickness data for SUDMEX-CONN'HCs
tmp=load('D:\OneDrive\1_Preprocessing\thickness136.mat');
y_HC = tmp.thickness(end-length(HC)+1:end, :);
y_HC=[y_HC;y_HCP];

W_score=parcelwise_Wscore(y_HC,y_CUD,x_HC,x_CUD);
