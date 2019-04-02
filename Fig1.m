clear all
close all
clc


del32Candidate = importdata('HazardDelta32Bootstrap.txt');
%first convert the hazard rate per year into survival probability till each
%age
%the result from each bootsrap takes 37 rows (each row is one year interval). 
%so there are 1000*37 rows.
%The four columns, representing the genotypes rateHetero (Aa),rateHomo (aa and AA),rateaa,rateAA
Survive = [];
for i = 1:1000;
   rate = del32Candidate(1+37*(i-1):37*i,:);
   survive = 1-rate(1,:);
   for j = 2:37;
       survive = [survive;survive(end,:).*(1-rate(j,:))];
   end
   Survive = [Survive;survive];
end
subplot('position',[0.15 0.19 0.32 0.73])
hold on
%order rateHetero (Aa),rateHomo (aa and AA),rateaa,rateAA
AA_survive = reshape(Survive(:,4),37,1000)';
Aa_survive = reshape(Survive(:,1),37,1000)';
aa_survive = reshape(Survive(:,3),37,1000)';
aaCIsurvive = sort(aa_survive(:,:));
aaCIhighsurvive = aaCIsurvive(975,:);
aaCIlowsurvive = aaCIsurvive(26,:);
% Correct the CI in the first two years because there is the no death in aa
% in the first several years. This step is not necessary. 
aaCIhighsurvive(1:2) = aaCIhighsurvive(3:4);
aaCIlowsurvive(1:2) = aaCIlowsurvive(3:4);
aaCIint = aaCIhighsurvive-aaCIlowsurvive;


text(44,0.89,'log-rank {\itP} = 0.0089','Fontsize',15);

rate = importdata('HazardDelta32Observed.txt');
%The order here is rateaa,rateHomo,rateHetero,rateAA, slightly different
%from the bootstrap file. 
survive = 1-rate(1,:);
for j = 2:37;
   survive = [survive;survive(end,:).*(1-rate(j,:))];
end
AA_curve = survive(:,4);
Aa_curve = survive(:,3);
aa_curve = survive(:,1);

plot([41:77],AA_curve,'.--r','MarkerSize',10);
plot([41:77],Aa_curve,'.--b','MarkerSize',10);
plot([41:77],aa_curve,'.--k','MarkerSize',10);
l = legend('+/+','\Delta32/+','\Delta32/\Delta32')
legend boxoff
set(l,'FontSize',12);
pValueAa = length(find(aa_survive(:,37)>=Aa_survive(:,37)))/1000
pValueAA = length(find(aa_survive(:,37)>=AA_survive(:,37)))/1000
xlim([41 77])
ylim([0.8 1])
set(gca,'YTick',[0.8:0.05:1]);
set(gca,'YTickLabel',[{'0.80'},{'0.85'},{'0.90'},{'0.95'},{'1.00'}]);
text(28.5,1,'a','FontSize',50);
box on
xlabel('Age');
ylabel('Survival probability')
set(gca,'FontSize',20);


subplot('position',[0.6 0.19 0.32 0.73])
A = importdata('MatchMAFSNPsForHWE.txt');
candidateSNP = 'rs62625034';% this is the rs number for delta32
indCand = find(strcmp(A.textdata(:,1),candidateSNP)==1);
Data = A.data;
Pminor = (Data(:,1)+Data(:,2)/2)./sum(Data,2);
aa_devi = Data(:,1)./sum(Data,2)./Pminor.^2;
h=histogram(aa_devi,30);
set(h,'FaceColor',[0.5 0.5 0.5]);
set(h,'EdgeColor',[0.5 0.5 0.5]);


x = [0.72 0.72];
y = [0.27 0.19];
annotation('textarrow',x,y,'String',[{'\Delta32/\Delta32           '},{'{\itP} = 0.0034           '}],'FontSize',15)
text(-0.3,0.9*5933,'b','FontSize',50);
ylim([0 0.9*5933]);
ylabel('Frequency');
set(gca,'XTickLabel',[{'-0.5'},{'0.0'},{'0.5'}]);
xlabel('\itF');
%set(gca,'XTick',[0.8132:0.5:1.5])
set(gca,'YTick',[0:0.3:0.9]*5933);
set(gca,'YTickLabel',[{'0.0'},{'0.3'},{'0.6'},{'0.9'}]);
set(gca,'FontSize',20);
box on 
set(gcf,'PaperPosition',[0 0 11 5])
saveas(1,'CCR5Fig1.png');

% This commented part could be used to estimate the confidence interval for
% surviavl probability. 
% 
% for i = 20:37
%     if (length(find(AA_survive(:,i) >= aa_survive(:,i)))) >= 950 %&& (length(find(Aa_survive(:,i) >= aa_survive(:,i)))>=950)
%         i + 41,length(find(AA_survive(:,i) >= aa_survive(:,i)))%,length(find(Aa_survive(:,i) >= aa_survive(:,i)))
%     end
% end
% 
% for i = 20:37
%     if (length(find(Aa_survive(:,i) >= aa_survive(:,i)))>=950)
%         i + 41,length(find(Aa_survive(:,i) >= aa_survive(:,i)))
%     end
% end
% 
% for i = 1:37
%     if (length(find(Aa_survive(:,i) >= AA_survive(:,i)))) >= 950 || (length(find(Aa_survive(:,i) >= AA_survive(:,i)))>=950)
%         i + 41,length(find(AA_survive(:,i) >= aa_survive(:,i)))
%     end
% end
% 
% Sur1 = sort((1-aa_survive(:,20:37))./(1- AA_survive(:,20:37)) -1);
% Sur2 = sort((1-aa_survive(:,20:37))./(1- Aa_survive(:,20:37)) -1);
% Sur3 = sort((1-Aa_survive(:,20:37))./(1- AA_survive(:,20:37)) -1);
% aa_sort = sort(aa_survive(:,20:37));
% Aa_sort = sort(Aa_survive(:,20:37));
% AA_sort = sort(AA_survive(:,20:37));
