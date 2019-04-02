clear all
close all
clc

AgeFile = importdata('AgeAndAncestryInformation.txt');
British = AgeFile(:,3);
indBritish = find(British==1);

AgeAtRecruit = AgeFile(:,5);
indAgeRecruit = find(AgeAtRecruit>0);

clear AgeFile

indBritishRecruit = intersect(indBritish,indAgeRecruit);
EnterEachAgeIndex = {}; %can be as old as something, but was not that old at recruitment

%The youngest British die at 40.8797, and the oldest British die at 78.5752, so we
%could calculate the death rate within interval 41 - 78. Maybe should use
%age starting 45-77 where there are more than 35 death in British.

for i = 40:1:69;
    indEnter = find(AgeAtRecruit >= i);
    indEnter = indEnter(find(AgeAtRecruit(indEnter) <i+5));
    EnterEachAgeIndex = [EnterEachAgeIndex,{intersect(indBritishRecruit,indEnter)}];
end



fid = fopen('CCRgenesVCF.txt');
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

rowThis = 12; %3,5,12,18
for i =1:rowThis;
    tline = fgetl(fid);
end

outGeno = strsplit(tline);
fclose(fid);
outGeno(1:10);
outGeno = outGeno(10:end);
indHetero = find(strcmp(outGeno,'0/1')==1);
indaa = find(strcmp(outGeno,'1/1')==1);
indAA = find(strcmp(outGeno,'0/0')==1);
Genotype = [];
for i = 1:30;
    Hetero = length(intersect(indHetero,EnterEachAgeIndex{i}));
    rateAA = length(intersect(indAA,EnterEachAgeIndex{i}));
    rateaa = length(intersect(indaa,EnterEachAgeIndex{i}));
    Genotype =[Genotype;rateaa,Hetero,rateAA];
end



Prob = Genotype./sum(Genotype,2)
p = Prob(:,1)+ Prob(:,2)/2;

bootHWE = [];
for j = 1:30;
    genotype = [zeros(1,Genotype(j,1)),ones(1,Genotype(j,2)),2*ones(1,Genotype(j,3))];
    N = sum(Genotype(j,:));
    HWE = [];
    for i = 1:1000;
        indRand = randsample(N,N,'true');
        bootGeno = genotype(indRand);
        ind0 = find(bootGeno == 0);
        ind1 = find(bootGeno == 1);
        pBoot = (length(ind0)+length(ind1)/2)/N;
        hwe = length(ind0)/(N*pBoot*pBoot);
        HWE = [HWE,hwe];
    end
    HWE = sort(HWE);
    bootHWE= [bootHWE;HWE(26),HWE(975)];
end


subplot('position',[0.15 0.19 0.32 0.73])
hold on
for i = 3:30;
    plot([39+i,39+i],bootHWE(i,:),'LineWidth',1,'color',[0.8 0.8 0.8]);
end
plot([42:69],Prob(3:end,1)./p(3:end)./p(3:end),'.-k','MarkerSize',20) ;

ylim([0.75 0.95])
set(gca,'YTick',[0.75:0.05:0.95]);
set(gca,'YTickLabel',[{'0.75'},{'0.80'},{'0.85'},{'0.90'},{'0.95'}]);
text(30,0.95,'a','FontSize',50);
box off
xlabel('Age');
ylabel('Observed deviation');
set(gca,'FontSize',20);


subplot('position',[0.6 0.19 0.32 0.73])

total = length(indaa)+length(indAA)+length(indHetero);
Frac = [length(indaa)/total,length(indHetero)/total,length(indAA)/total];
HWE = [];
hUKB = importdata('HazardDelta32Observed.txt'); %rateaa,rateHomo,rateHetero,rateAA
hUKBAll = importdata('HazardBritishnonBritishData.txt');
hUKBAll = hUKBAll(:,1);%The british hazard rates
hUK = importdata('deathRateUK.txt');%39-78
hUK = hUK(3:end); 
c = hUKBAll./hUK;
hUKB = hUKB(:,[1,3,4])./c;
for i = 1:length(c);
    aa_Surv(i) = prod(1-hUKB(1:i,1));
    Aa_Surv(i) = prod(1-hUKB(1:i,2));
    AA_Surv(i) = prod(1-hUKB(1:i,3));
end

fid = fopen('Delta32SurvivalProbCorrected.txt','w');
fprintf(fid,'age\t');
fprintf(fid,'%d\t',[41:77]);
fprintf(fid,'\naa_survival\t');
fprintf(fid,'%d\t',aa_Surv);
fprintf(fid,'\nAa_survival\t');
fprintf(fid,'%d\t',Aa_Surv);
fprintf(fid,'\nAA_survival\t');
fprintf(fid,'%d\t',AA_Surv);
fprintf(fid,'\n');
fclose(fid);

for i = 3:30;
    Fracnew = Frac .* [aa_Surv(i),Aa_Surv(i),AA_Surv(i)]/sum(Frac .* [aa_Surv(i),Aa_Surv(i),AA_Surv(i)]);
    pNew = Fracnew(1)+Fracnew(2)/2;
    HWE = [HWE; Fracnew./[pNew^2,pNew*(1-pNew)*2,(1-pNew)^2]];
end

hold on
plot([43:70],HWE(:,1),'.-k','MarkerSize',20);
ylim([0.84 0.87])
set(gca,'YTick',[0.84:0.01:0.87]);
set(gca,'YTickLabel',[{'0.84'},{'0.85'},{'0.86'},{'0.87'}]);
xlabel('Age');
ylabel('Predicted deviation');
text(30,0.87,'b','FontSize',50);
box off
set(gca,'FontSize',20);
set(gcf,'PaperPosition',[0 0 11 5])
saveas(1,'CCR5FigS1.png');
[r1 p1] = corr(Prob(3:end,1)./p(3:end)./p(3:end),HWE(:,1),'type','Spearman')
