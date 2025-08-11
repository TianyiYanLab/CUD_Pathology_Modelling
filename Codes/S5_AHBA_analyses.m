% ---------------------------------------------------------------------------------------
% Refrence: Morgan et al. PNAS, 2019.
% https://github.com/SarahMorgan/Morphometric_Similarity_SZ
% Adapted by Ziteng Han, 2024.
% ---------------------------------------------------------------------------------------

%% input data
nregs=400; % number of regions

% AHBA data, output of abagen.ipynb; Predictors
X = readtable(['D:\OneDrive\ABHA\AHBA_expression_data_S400Net7.csv'], ...
    'VariableNamingRule','preserve');   
geneData = X(:, 2:end);
geneMatrix = table2array(geneData);
columnNames  =  X.Properties.VariableNames;
geneNames  =  columnNames(2:end);   %Gene names

%epicenter map; Response variable
Map_dir =('D:\OneDrive\5_TMS_response\Epicenter\SC') ;
map=load([Map_dir '\epicenterMap_TMS_avgranking_S4007Net_4modes_dataHansenSSC.mat']);
Y = map.my_array';  

% z-score normalized
X = zscore(geneMatrix);  % Predictor variables (gene expression)
Y = zscore(Y);     % Response variable 

%% Setp1: perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
dim=15;

plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on

%%% plot correlation of PLS component 1 with t-statistic:
figure
plot(XS(:,1),map.my_array','r.')
[R,p]=corr(XS(:,1),map.my_array') 
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('epicenter ranking','FontSize',14);
grid on

%% Step2: permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included

rep=1000;
for dim=1:10
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    temp=cumsum(100*PCTVAR(2,1:dim));
    Rsquared = temp(dim);
        for j=1:rep
            %j
            order=randperm(size(Y,1));
            Yp=Y(order,:);
    
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
    
            temp=cumsum(100*PCTVAR(2,1:dim));
            Rsq(j) = temp(dim);
        end
    dim
    R(dim)=Rsquared
    p(dim)=length(find(Rsq>=Rsquared))/rep
end
figure
plot(1:dim, p,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on

%% Step3: bootstrap to get the gene list

genes=geneNames; % this needs to be imported first
geneindex=1:length(genes);

%number of bootstrap iterations:
bootnum=1000;

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],map.my_array');

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);

%print out results
csvwrite('PLS1_ROIscores.csv',XS(:,1));
csvwrite('PLS2_ROIscores.csv',XS(:,2));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

% print out results
% later use first column of these csv files for pasting into Metescape (for
% bootstrapped ordered list of genes) 
fid1 = fopen('PLS1_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1)

fid2 = fopen('PLS2_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2)