clc;clear;
close all;
% profile on;
%%
addpath('/home/s1785969/RDS/MATLAB/DAtoolbox');
addpath('/home/s1785969/RDS/MATLAB/BlockMatchingAlgoMPEG');
addpath('/home/s1785969/RDS/MATLAB/');
addpath('/home/s1785969/RDS/MATLAB/DA_Image_REG');
addpath('/home/s1785969/RDS/MATLAB/DicomGUI');
addpath('/home/s1785969/RDS/MATLAB/DicomGUI/Scripting');
addpath('/home/s1785969/RDS/MATLAB/ImgOutputs/blar_hyper_run_ML_3D/newrun/run11');
cd /home/s1785969/RDS/MATLAB/ImgOutputs/blar_hyper_run_ML_3D/newrun/run11
%%
mypath='/home/s1785969/RDS/MATLAB/ImageDB/Bladder';
patfolder=dir(mypath);
patfolder(1:3,:) = []; %Removing .,.., directory path along with the blar 1 since it is validation data
bigfolderlen=length(patfolder);
subfolder=cell(bigfolderlen,1);

for folderi=1:bigfolderlen  
    subfolder{folderi,1}=dir(strcat(mypath,'/',patfolder(folderi).name));
    Names=extractfield(subfolder{folderi,1},'name');
    Names=Names';
    index=find(contains(Names,'.'));
    subfolder{folderi,1}(index)=[];
    Names(index)=[];
end

%% Gradient Search Parameters
lambda1=-0.5;%x initial point
lambda2=0.5;%y initial point
dl1=0.0001;% x-step size
dl2=0.0001;% y-step size
alpha=0.001;% learning rate only after identifying proper learning
% rate
beta=0.9;% momentum hyperparameter
epochnumber=1; % No. of Epochs
NumIter=1000; % No. of Iterations
Grad=[Inf,Inf];% Gradient inital
Tol=10^-16; % Tolerance
NormGrad=Inf; % Norm of Gradient
Precision=10^-16;
StepsizeDiff=1;

%Initialisation
DSCL11=zeros(NumIter,epochnumber);
DSCL12=zeros(NumIter,epochnumber);
DSCL21=zeros(NumIter,epochnumber);
DSCL22=zeros(NumIter,epochnumber);
GradNum=zeros(NumIter,epochnumber*2);
StepNum=zeros(NumIter,epochnumber*2);
vel=zeros(NumIter,epochnumber*2);
GradNom=zeros(NumIter,epochnumber);
Mlambda1=zeros(NumIter,epochnumber);
Mlambda2=zeros(NumIter,epochnumber);
Malpha=zeros(NumIter,epochnumber);
%%
% LogAlpha=logspace(-8,1,NumIter);
% LogAlpha=LogAlpha';
runerrorflag=1;
for epoch=1:epochnumber

PlanScanName=cell(NumIter,1);
TretScanName=cell(NumIter,1);
PatName=cell(NumIter,1);
PathFolder=cell(NumIter,1);

Li=1; % Gradient Stepping Count

% while Li<=NumIter && StepsizeDiff>Precision && NormGrad>Tol
% while Li<=NumIter && NormGrad>Tol
while Li<=NumIter
%     if runerrorflag==1
%         load('1-157.mat');
% %         Li=Li+1;
%         runerrorflag=0;
%     end

bigfolderind=randi([1 bigfolderlen]);
subfolderlen=length(subfolder{bigfolderind,1});
subfolderind=randi([2 subfolderlen]);

PlanScanName{Li,1}=subfolder{bigfolderind,1}(1,1).name;
TretScanName{Li,1}=subfolder{bigfolderind,1}(subfolderind,1).name;
PatName{Li,1}=patfolder(bigfolderind).name;
PathFolder{Li,1}=subfolder{bigfolderind,1}(subfolderind,1).folder;


directory_name1=char(strcat(PathFolder{Li,1},'/',PlanScanName{Li,1}));
directory_name2=char(strcat(PathFolder{Li,1},'/',TretScanName{Li,1}));
 
                Lambda1plusdel=lambda1+(dl1/2);
                Lambda1miusdel=lambda1-(dl1/2);
                Lambda2plusdel=lambda2+(dl2/2);
                Lambda2miusdel=lambda2-(dl2/2);
                
[~,~,~,DSCGrad11Lam1,~]=BLMIR_1(directory_name1,directory_name2,Lambda1miusdel,lambda2,1);
[~,~,~,DSCGrad12Lam1,~]=BLMIR_1(directory_name1,directory_name2,Lambda1plusdel,lambda2,1);
[~,~,~,DSCGrad21Lam1,~]=BLMIR_1(directory_name1,directory_name2,lambda1,Lambda2miusdel,1);
[~,~,~,DSCGrad22Lam1,~]=BLMIR_1(directory_name1,directory_name2,lambda1,Lambda2plusdel,1);
                      
        DSCL11(Li,epoch)=1-DSCGrad11Lam1;
        DSCL12(Li,epoch)=1-DSCGrad12Lam1;
        DSCL21(Li,epoch)=1-DSCGrad22Lam1;
        DSCL22(Li,epoch)=1-DSCGrad21Lam1;
        
        %Gradient - Successive Approximation
        GradLam1=(DSCL12(Li,epoch)-DSCL11(Li,epoch))/dl1;
        GradLam2=(DSCL22(Li,epoch)-DSCL21(Li,epoch))/dl2;                

%        Gradient Vector Past value for momentum
        if Li==1
            vel1_P=0;
            vel2_P=0;
            Grad1_P=0;
            Grad2_P=0;
            Step1_P=0;
            Step2_P=0;
        else
            vel1_P=vel(Li-1,epoch*2-1);
            vel2_P=vel(Li-1,epoch*2);
            Step1_P=StepNum(Li-1,epoch*2-1);
            Step2_P=StepNum(Li,epoch*2);
        end
 
        
        vel1C=(1-beta)*GradLam1+(beta*vel2_P);
        vel2C=(1-beta)*GradLam2+(beta*vel2_P);
        
        Step1_C=(alpha*vel1C);
        Step2_C=(alpha*vel2C);
        
        lambda1=lambda1-Step1_C;
        lambda2=lambda2-Step2_C;           
        
        GradNum(Li,epoch*2-1)=GradLam1;
        GradNum(Li,epoch*2)=GradLam2;

        StepNum(Li,epoch*2-1)=Step1_C;
        StepNum(Li,epoch*2)=Step2_C;
        
        vel(Li,epoch*2-1)=vel1C;
        vel(Li,epoch*2)=vel2C;
        
        Grad=[GradLam1,GradLam2];
        NormGrad=norm(Grad);        
        GradNom(Li,epoch)=NormGrad;
        
%         Malpha(Li,epoch)=alpha;
        Mlambda1(Li,epoch)=lambda1;
        Mlambda2(Li,epoch)=lambda2;        
        


StepsizeDiff=norm([Step1_C Step2_C]-[Step1_P Step2_P]); % Convergence Condition check



% save(matfil,'DSCGrad11Lam1','DSCGrad12Lam1','DSCGrad21Lam1','DSCGrad22Lam1'...
%     ,'DSCL11','DSCL12','DSCL21','DSCL22','Grad','StepsizeDiff','GradLam1',...
%     'GradNom','GradNum','StepNum','vel','Lambda1miusdel','Lambda1plusdel','Lambda2miusdel','Lambda2plusdel','Li','Malpha','Mlambda1','Mlambda2','Names'...
%     ,'NormGrad','NumIter','PatName','PathFolder','PlanScanName'...
%     ,'Tol','TretScanName','alpha','bigfolderind','bigfolderlen'...
%     ,'directory_name1','directory_name2','dl1','dl2','epoch','epochnumber','folderi','index','lambda1'...
%     ,'lambda2','matfil','mypath','patfolder','subfolder'...
%     ,'subfolderind','subfolderlen');
matfil=[num2str(epoch),'-',num2str(Li)];
Li=Li+1;
save(matfil);



end

end
%%
%% Functions 
function [Cr,Rel2,Rel1,DSC2,DSC1]=BLMIR_1(directory_name1,directory_name2,lambda1A,lambda2,~)
[IP1a,IT1a,MP1a,MT1a,Z1a,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
ID1x=find(TF1x);
ID2x=find(TF2x);
% MP2=MP1{ID1x,1};
% MT2=MT1{ID2x,1};
MP2a=MP1a{ID1x,1};
MT2a=MT1a{ID2x,1};
% [IP1a,IT1a,MP2a,MT2a,Z1a,ZLa]=rearrangeslice(IP1a,IT1a,MP2a,MT2a,Z1a); % Common no. of slices selection
%     if flag==0
%     [sr,sc,sp]=size(IP1a);
%     IP1a=randi(2000,[sr,sc,sp],'int16');
%     MP2a=false(sr,sc,sp);
%     end
m=500;
n=500;
indSize=3; % IndSize must always be odd for the symmetrical Z-directional slices
mbSize=5; % Block Size
p=5; % Search Parameter
% lambda1A=0.25;% Distance Regularisation
% lambda2=0.005;
costV=2^16+1;
% lambda1A=[0.25 0.5 1 5 10];
% lambda2=[0.25 0.5 1 5 10];
% lambda1A=lambda1A';
% lambda2=lambda2'; % Orientation Regularisation
% LamL=length(lambda2);
LenX=m/mbSize; % Block col index
LenY=n/mbSize; % Block row index
[Cr,Rel2,Rel1,DSC2,DSC1]=BLMIR3DDiceFx_24_inbuilt(Z1a,IP1a,IT1a,MP2a,MT2a,indSize,mbSize,p,lambda1A,lambda2,m,n,LenX,LenY,costV);%1 lambda
end
%%
function [im1,im2,MaskV1,MaskV2,Z1,ZL,SN1,SN2]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,Imgmodal)
%Load data 1 with structure names
myhandle1=Open_DicomGui();
[~,~]=DG_load_data(myhandle1,directory_name1,Imgmodal);
[SN1,~] = DG_get_names(myhandle1);
SN1=SN1';
[im1,Z1,~]=DG_get_image_data(myhandle1);

%Load data 2 with structure names
myhandle2=Open_DicomGui();
[~,~]=DG_load_data(myhandle2,directory_name2,Imgmodal);
[SN2,~] = DG_get_names(myhandle2);
SN2=SN2';
[im2,~,~]=DG_get_image_data(myhandle2);


[~,~,ZL]=size(im1);

SN1L=length(SN1);
SN2L=length(SN2);

MaskV1=cell(SN1L,1);
for contouri=1:SN1L
    [MaskV1{contouri,1},~]=DG_generate_volume_mask(myhandle1, contouri);
end

MaskV2=cell(SN2L,1);
for contouri=1:length(SN2)
    [MaskV2{contouri,1},~]=DG_generate_volume_mask(myhandle2, contouri);
end
Close_DicomGui(myhandle2);
Close_DicomGui(myhandle1);
end
%% Dice score, Reliability score & Block matching for all slices using parfor
function [Cr,ECV1,ECV,DSC2,DSC1]=BLMIR3DDiceFx_24_inbuilt(Z1,IP1,IT1,MP2,MT2,indSize,mbSize,p,lambda1A,lambda2,m,n,LenX,LenY,costV)
ZL=length(Z1); % Total length of slices per scan
CI1=false(m,n,ZL);
MVC1=zeros(ZL,1);
MVT1=zeros(ZL,1);
MVI1=zeros(ZL,1);
MTC1=false(m,n,ZL);
CI2=false(m,n,ZL);
MVC2=zeros(ZL,1);
MVT2=zeros(ZL,1);
MVI2=zeros(ZL,1);
MTC2=false(m,n,ZL);
E=zeros(ZL,1);
E1=zeros(ZL,1);
C=zeros(ZL,1);
    parfor ind=1:ZL
    C(ind,1)=corr2(IP1(:,:,ind),IT1(:,:,ind));   
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1(indz);
    ZLocV=(Zloc-Zloc(indzT==0))*-10; % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
    imgTreatment=IT1(:,:,ind);
    imgPlanning=IP1(:,:,indz);
    % masTreatment=MT2(:,:,ind);
    masPlanning=MP2(:,:,indz);
    mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV); % 1st Pass
    mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda2,m,n,LenX,LenY,costV); % 2nd Pass
    MTE1=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,indL,mbSize);% Mask prediction - 1st Pass
    MTE2=motioncomp3d_1_inbuilt(masPlanning,mV2,m,n,indL,mbSize);% Mask Prediction - 2nd Pass
    MTE1=dustthemask(MTE1);% Mask Correction
    MTE2=dustthemask(MTE2);
    % ITE=motioncomp3d_1(imgPlanning,mV2,m,n,indL,mbSize);
    E(ind,1)=sum(mV1(:,7));
    E1(ind,1)=sum(mV2(:,7));
                MTC1(:,:,ind)=logical(MTE1);
                MT2m1=imresize(MT2(:,:,ind),[m n]);              
                CI1(:,:,ind)=MT2m1&MTC1(:,:,ind);            
                MVI1(ind,1)=length(find(CI1(:,:,ind)));
                MVC1(ind,1)=length(find(MTC1(:,:,ind)));
                MVT1(ind,1)=length(find(MT2m1));                
                MTC2(:,:,ind)=logical(MTE2);
                MT2m2=imresize(MT2(:,:,ind),[m n]);              
                CI2(:,:,ind)=MT2m2&MTC2(:,:,ind);            
                MVI2(ind,1)=length(find(CI2(:,:,ind)));
                MVC2(ind,1)=length(find(MTC2(:,:,ind)));
                MVT2(ind,1)=length(find(MT2m2));
    end         
        VC1=sum(MVC1);
        VT1=sum(MVT1);
        VI1=sum(MVI1);
        DSC1=(2*VI1)/(VT1+VC1);      
        VC2=sum(MVC2);
        VT2=sum(MVT2);
        VI2=sum(MVI2);
        DSC2=(2*VI2)/(VT2+VC2);       
        ECV1=sum(E1)/(ZL*LenX*LenY*costV);
        ECV=sum(E)/(ZL*LenX*LenY*costV);
        Cr=mean(C);
end
%% Block matching - First Pass
function mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV)
% [m,n]=size(imgTreatment);
it=imresize(imgTreatment,[m n]);
mp=int16(zeros(m,n,indL));
    for si=1:indL
    mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
imgTreatment=it;
imgPlanning=mp;
% costV=65537;%Cost Value
% costV=2^16+1;
prow=2*p + 1;
pcol=2*p + 1;
costs1 = ones(prow, prow, indL) * costV;%Cost matrix
dist = zeros(prow, prow, indL);
mbCount = 1;
mbSize_N=mbSize*mbSize;
mV1 = zeros(round(m*n/mbSize^2),7); % 3D motion vector + Z index
blkfactor=1;
for i = 1 : mbSize : m-(blkfactor*mbSize)+1
    for j = 1 : mbSize : n-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
        imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block            
            imgPsubset1=imgPlanning1(refBlkVer:refBlkVer+mbSize-1,:);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate                                
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
%                 costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N);
                costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N)+ dist(m1+p+1,n1+p+1,zi);           
            end
        end
    end        
        [zcostmin,Zminloc]=min(costs1(:));
        [zr,zc,zp]=ind2sub(size(costs1),Zminloc);        
        mV1(mbCount,1) = (zr-p-1);    % row1 co-ordinate for the vector
        mV1(mbCount,2) = (zc-p-1);    % col1 co-ordinate for the vector
        mV1(mbCount,3) = ZLocV(zp);
        mV1(mbCount,4) = indz(zp);
        mV1(mbCount,5) = zp;      
        P1=[m1 n1 ZLocV(zi)];
%         P1=[floor(i+mbSize/2) floor(j+mbSize/2) ZLocV(2)];
        P2=[mV1(mbCount,1) mV1(mbCount,2) mV1(mbCount,3)];
        ThetainSlope=dot(P1,P2)/norm(cross(P1,P2));
                if isnan(ThetainSlope)==1 
                    ThetainSlope=0; 
                end
        mV1(mbCount,6)=atand(ThetainSlope);  
%         mV1(mbCount,7)=costV-zcostmin;
        mV1(mbCount,7)=zcostmin;     % Reliability Score
        mbCount = mbCount + 1;
        costs1 = ones(prow, pcol, indL) * costV;
        dist = zeros(prow, prow, indL);
    end
end 


end
%% Block matching - Second Pass
function mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda2,m,n,LenX,LenY,costV)
it=imresize(imgTreatment,[m n]);
mp=int16(zeros(m,n,indL));
    for si=1:indL
    mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
imgTreatment=it;
imgPlanning=mp;
Vthe=mV1(:,6);
Vthe1=reshape(Vthe,[LenX,LenY]);
% OmegTh=omegdeiff(Vthe1,LenX,LenY);
OmegTh=omegdiif_DrDave(Vthe1,LenX,LenY);
% BlocOd=omegdiif_DrDave(omeg,LenX,LenY)
% costV=65537;%Cost Value
% costV=2^16+1;%Cost Value
prow=2*p + 1;
pcol=2*p + 1;
costs1 = ones(prow, prow, indL) * costV;%Cost matrix
dist = ones(prow, prow, indL);
blkfactor=1;
mbSize_N=mbSize*mbSize;
mV2 = zeros(round(m*n/mbSize^2),7); % 3D motion vector + Z index
mbCount=1;
    for i = 1 : mbSize : m-(blkfactor*mbSize)+1
    for j = 1 : mbSize : n-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
        imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block 
            br = floor((refBlkVer-1)/mbSize)+1; %row vector for Omeg matrix
            imgPsubset1=imgPlanning1(refBlkVer:refBlkVer+mbSize-1,:);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate  
                bc = floor((refBlkHor-1)/mbSize)+1;%col vector
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
                  costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N) + dist(m1+p+1,n1+p+1,zi) + (lambda2*OmegTh(br,bc));
            end
        end
    end        
        [zcostmin,Zminloc]=min(costs1(:));
        [zr,zc,zp]=ind2sub(size(costs1),Zminloc);     
        mV2(mbCount,1) = (zr-p-1);    % row1 co-ordinate for the vector
        mV2(mbCount,2) = (zc-p-1);    % col1 co-ordinate for the vector
        mV2(mbCount,3) = ZLocV(zp);
        mV2(mbCount,4) = indz(zp);
        mV2(mbCount,5) = zp;       
%         mV2(mbCount,7)=costV-zcostmin;
        mV2(mbCount,7)=zcostmin; % Reliability Score 
        mbCount = mbCount + 1;
        costs1 = ones(prow, pcol, indL) * costV;
        dist = zeros(prow, prow, indL);
    end
    end 
end
%% Motion Compensated Image Generation
function imageComp=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,indL,mbSize)
% [sr,sc,~]=size(masPlanning);
mp=zeros(m,n,indL);
    for si=1:indL
    mp(:,:,si)=imresize(masPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
masPlanning=mp;
imageComp=zeros(m,n);
mbCount = 1;
    for mi = 1:mbSize:m-mbSize+1
    for mj = 1:mbSize:n-mbSize+1      
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % this means we are scanning in order       
        dy = mV1(mbCount,1);
        dx = mV1(mbCount,2);
        refBlkVerm = mi + dy;
        refBlkHorm = mj + dx;
        refBlkDeptm = mV1(mbCount,5);
        imageComp(mi:mi+mbSize-1,mj:mj+mbSize-1) = masPlanning(refBlkVerm:refBlkVerm+mbSize-1, refBlkHorm:refBlkHorm+mbSize-1,refBlkDeptm);
            mbCount = mbCount + 1;
    end
    end
% imageComp=imresize(imageComp,[sr sc]);
end

%% Mask Cleaning/Smoothening/Dusting
function Mn=dustthemask(MTE)
se = strel('disk',4);
Mn=imopen(MTE,se);
Mn=imclose(Mn,se);
Mn=imfill(Mn,'holes');
end