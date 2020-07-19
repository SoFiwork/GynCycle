function [FSHVec, StartVec] = CreateFollicles(para,paraPoi,tb,te)

%create normal distributed fsh sensitivities for each foll
fileID2 = fopen('FSH.txt','w+');
fprintf(fileID2,'Number    FSH\n');
FSHdistri = makedist('Normal','mu',para(7),'sigma',para(8));
for i=1:1000
    fsh = random(FSHdistri);
    fprintf(fileID2,'%f %f\n',i,fsh);
end
fclose(fileID2);


%create poisson distributed starting times for each foll
fileID = fopen('StartTimesPoiss.txt','w+');
fprintf(fileID,'Start time\n');
TotIntervall = te-tb;
timevec=poissonproc(paraPoi(1),[tb,te]); 
arraysize=size(timevec);
for i=1:arraysize
    fprintf(fileID,'%f \n',timevec(i));
end
fclose(fileID);


%load StartNumbers and FSH Sensitivities from File
file = 'StartTimesPoiss.txt';
file2 = 'FSH.txt';
delimiterIn=' ';
headerlinesIn=1;
data=importdata(file,delimiterIn,headerlinesIn);
data2=importdata(file2,delimiterIn,headerlinesIn);
NumValStart = size(data.data(1:end,1));
NumValStart = max(NumValStart);
StartVec = zeros(NumValStart,1);
for i = 1:NumValStart
    StartVec(i) = data.data(i,1);
end
NumValFSH = size(data2.data(1:end,1));
NumValFSH = max(NumValFSH);
FSHVec = zeros(NumValFSH,1);
for i = 1:NumValFSH
    FSHVec(i) = data2.data(i,2);
end

end
