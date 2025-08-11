function S9_PLS_input_preparation(pheno_mat,version,epicenter_version)
% Inputs:
%    - pheno_mat：A #Subjects x #behaviours matrix
%    - version：PLS analysis index (e.g., v80_v2)
%    - epicenter_version：Epicenter map index (e.g., 4modes_dataHansenSSC)
% Outputs:
%    Input files for the S8_Behavioral_PLS_analysis.ipynb

    wdir=['D:\OneDrive\5_TMS_response\data_cud'];   
    
    rows_with_NaN =any(isnan(pheno_mat), 2);
    phenos=pheno_mat(~rows_with_NaN ,:);

    subs53_list='D:\OneDrive\SUDMEX_TMS\SUDMexTMS.txt';
    fileID = fopen(subs53_list, 'r');
    sub_list = textscan(fileID, '%s');  
    fclose(fileID); 
    subs53 = sub_list{1}; 

    subs = subs53(~rows_with_NaN,:);
    writecell(subs, [wdir '\subject_list_' version '_' epicenter_version '.txt']);

    %regress age and sex from the phenotypes 
    age_sex = extract_subject_info([wdir '\subject_list_' version '_' epicenter_version '.txt'], ...
        'D:\fMRI_Data\SUDMEX_TMS\10409461\Demographic_data.xlsx',{'q1_age','q1_sex'},'rmNaN');
    [phenos_hat, ~] = CBIG_regress_X_from_y_train(phenos, zscore(age_sex));
    save([wdir '\phenos_TMS_' version '_' epicenter_version '.mat'],'phenos_hat')
    
    tmp=load('D:\OneDrive\5_TMS_response\CT_Wscore53_S400_7Net_vHCP2.mat');
    Wscore74_s400_7Net =tmp.Wscore_S400_7Net;

    [~, idx] = ismember(subs, subs53);
    Wscores_S4007Net = Wscore74_s400_7Net(idx,:);
    savemat=([wdir '\ct_Wscore_' version '_' epicenter_version '_s400_7Net.mat']);
    save(savemat,'Wscores_S4007Net')

    epiMaps_dir=['D:\OneDrive\5_TMS_response\Epicenter\SC\indi_epicenterMap_all53_' epicenter_version];
    epiMaps=[];
    for i=1:length(idx)
        thisind=idx(i);
        tmp=load([epiMaps_dir '\epicenterMap_avgranking_S4007Net_' num2str(thisind) 'th.mat']);
        epiMaps(i,:)=tmp.my_array;
    end
    save([wdir '\epicenterMap_avgranking_' version '_' epicenter_version '_S4007Net.mat'],'epiMaps')
