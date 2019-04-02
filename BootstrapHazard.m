clear all
close all
clc



AgeFile = importdata('AgeAndAncestryInformation_version2.txt');
British = AgeFile(:,3);
indBritish = find(British==1);
daysAtRecruit = AgeFile(:,4);
indBirth = find(daysAtRecruit>0);
AgeAtRecruit = AgeFile(:,5);

%-----deal with the approximation of using 15th of month as bd------------
AgeFile0 = importdata('AgeAndAncestryInformation.txt');
AgeRecruitOriginal = AgeFile0(:,5);
ind = find(AgeAtRecruit < AgeRecruitOriginal);
AgeAtRecruit(ind) = AgeRecruitOriginal(ind);

indProbmatic = find(AgeAtRecruit > AgeRecruitOriginal+2); %17 are problematic
%These are the ones should be removed perhaps, because age exit is also
%inaccurate
%among them one dead, and two are British
%the two British did not die

ind = find(AgeAtRecruit > AgeRecruitOriginal+1);
AgeAtRecruit(ind) = AgeRecruitOriginal(ind)+1;
%----deal with the approximation of using 15th of month as bd------------



indAgeRecruit = find(AgeAtRecruit>0);

AgeExit = (datenum(2016,02,16) - daysAtRecruit)/365.25 + AgeAtRecruit;
AgeAtDeath = AgeFile(:,6);
indDeath = find(AgeAtDeath>0);
clear AgeFile

indBritishRecruit = intersect(indBritish,indBirth);
indBritishRecruit = setdiff(indBritishRecruit,intersect(indProbmatic,indBritishRecruit));

l_boot = length(indBritishRecruit);

fidOut = fopen('HazardDelta32Bootstrap.txt','w');
rowThis = 12; %3,5,12,18
fid = fopen('CCRgenesVCF.txt');
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

for rowpass = 1:rowThis;
    tline = fgetl(fid);
end
outGeno = strsplit(tline);
outGeno(1:10)
outGeno = outGeno(10:end);
indHetero = find(strcmp(outGeno,'0/1')==1);
indaa = find(strcmp(outGeno,'1/1')==1);
indAA = find(strcmp(outGeno,'0/0')==1);
indNotMissing = find(strcmp(outGeno,'./.')==0);
indBritishRecruit = intersect(indBritishRecruit,indNotMissing);
indaa_recruit = intersect(indaa,indBritishRecruit);
indAA_recruit = intersect(indAA,indBritishRecruit);
indHetero_recruit = intersect(indHetero,indBritishRecruit);
l_aa = length(indaa_recruit);
l_Hetero = length(indHetero_recruit);
l_AA = length(indAA_recruit);

fclose(fid); 
SurviveCCR = [];
for rep = 1:1000
    indBritishRecruit_boot = [randsample(indaa_recruit,l_aa,true);...
        randsample(indAA_recruit,l_AA,true);randsample(indHetero_recruit,l_Hetero,true)];
    DeathIndexEachAge ={}; %death during age i to i+1
    EnterEachAgeIndex = {}; %can be as old as something, but was not that old at recruitment
    
    %The youngest British die at 40.8797, and the oldest British die at 78.5752, so we
    %could calculate the death rate within interval 41 - 78. Maybe should use
    %age starting 45-77 where there are more than 35 death in British.

    for i = 41:77;
        indEnter = find(AgeExit >= i+1);
        indEnter = indEnter(find(AgeAtRecruit(indEnter) <i));
        indEnter = intersect(indEnter,[find(AgeAtDeath >= i); find(AgeAtDeath==0)]);
        EnterEachAgeIndex = [EnterEachAgeIndex,{indBritishRecruit_boot(ismember(indBritishRecruit_boot,indEnter))}];
        indDeathi = find(AgeAtDeath >= i);
        indDeathi = indDeathi(find(AgeAtDeath(indDeathi)<i+1));
        DeathIndexEachAge = [DeathIndexEachAge,{indBritishRecruit_boot(ismember(indBritishRecruit_boot,intersect(indDeathi,indEnter)))}];
    end
    DeathRate = [];
    for i = 1:35;
        rateHetero = sum(ismember(DeathIndexEachAge{i},indHetero))/sum(ismember(EnterEachAgeIndex{i},indHetero));
        rateaa = sum(ismember(DeathIndexEachAge{i},indaa))/sum(ismember(EnterEachAgeIndex{i},indaa));
        rateAA = sum(ismember(DeathIndexEachAge{i},indAA))/sum(ismember(EnterEachAgeIndex{i},indAA));
        rateHomo = (sum(ismember(DeathIndexEachAge{i},indaa))+sum(ismember(DeathIndexEachAge{i},indAA)))/...
                (sum(ismember(EnterEachAgeIndex{i},indaa))+sum(ismember(EnterEachAgeIndex{i},indAA)));
        DeathRate =[DeathRate;rateHetero,rateHomo,rateaa,rateAA];
    end
    i = 36;
    rateHetero = (sum(ismember(DeathIndexEachAge{i},indHetero))+sum(ismember(DeathIndexEachAge{i+1},indHetero)))/...
        (sum(ismember(EnterEachAgeIndex{i},indHetero))+sum(ismember(EnterEachAgeIndex{i+1},indHetero)));
    rateaa = (sum(ismember(DeathIndexEachAge{i},indaa))+sum(ismember(DeathIndexEachAge{i+1},indaa)))/...
        (sum(ismember(EnterEachAgeIndex{i},indaa))+sum(ismember(EnterEachAgeIndex{i+1},indaa)));
    rateAA = (sum(ismember(DeathIndexEachAge{i},indAA))+sum(ismember(DeathIndexEachAge{i+1},indAA)))/...
        (sum(ismember(EnterEachAgeIndex{i},indAA))+sum(ismember(EnterEachAgeIndex{i+1},indAA)));
    rateHomo = (sum(ismember(DeathIndexEachAge{i},indaa))+sum(ismember(DeathIndexEachAge{i},indAA))+...
        sum(ismember(DeathIndexEachAge{i+1},indaa))+sum(ismember(DeathIndexEachAge{i+1},indAA)))/...
            (sum(ismember(EnterEachAgeIndex{i},indaa))+sum(ismember(EnterEachAgeIndex{i},indAA))+...
            sum(ismember(EnterEachAgeIndex{i+1},indaa))+sum(ismember(EnterEachAgeIndex{i+1},indAA)));
    DeathRate =[DeathRate;rateHetero,rateHomo,rateaa,rateAA;rateHetero,rateHomo,rateaa,rateAA];
    SurviveCCR = [SurviveCCR;prod(1-DeathRate)];
    for i = 1:37;
       % SurviveOut = [SurviveOut;prod(1-Collect((i-1)*37+1:i*37,:))];
        fprintf(fidOut,'%d\t',DeathRate(i,:));
        fprintf(fidOut,'\n');
    end
end

fclose(fidOut);
