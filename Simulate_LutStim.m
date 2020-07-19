% Calculation and solving of the ODE
% needs (Parameter(programming, Poisson distribution, ODE calculations),
% integration starttime, endtime, initial values, ShowTime,
% poisson distributed starttimes of the follicles, normal distributed FSH
% sensitivities of the follicles, ShowStuff, SavePics, PicFile, tshow)
function [res, New, NoNew] = NoCohortSimH_LutStim(para,paraPoi,paraOde,Par,tb,te,StartValues,StartTimes,FSHVec,ShowStuff,SaveStuff,DirStuff)

%integration period
tspan=[tb,te];
%variable for the current time
t=tb;
%Timepoint of last ovulation
Tovu = 14;
%
%-----------------------------------------------------------------------
%
%initial values
y0Foll = StartValues(1);    %startsize of the follicle
y0E = StartValues(3);       %startvalue estradiol concentration in the ovaries 
y0P4 = StartValues(5);      %startvalue progesterone
y0LH = StartValues(7);      %startvalue LH
y0FSH = StartValues(9);     %startvalue FSH

y0 = StartValues';
%
%-----------------------------------------------------------------------
%
FollCounter = 1;                                            %index of the follicle we currently working on

Follicles = FollicleClass(y0Foll,FSHVec(FollCounter),t);    %class to save follicles and their properties
FollCounter = FollCounter + 1;
    
TimeCounter=1;
NextStart=StartTimes(TimeCounter);

%arrays to save times when new follicles emerge or when they can't emerge
NewFollicle = [];
NoNewFollicle = [];

%parameter for the equations-function
LastYValues = [];

%vector for max sizes of the foll
FollMax =[];
FollMaxTime=[];

%vector for #foll with a certain size
Foll10mm =[];
Foll15mm =[];
count15mm=0;
%
%-----------------------------------------------------------------------
%
%variables for drug Administration
dd1 = 0;
dosing_events1 = [];
%
%-----------------------------------------------------------------------
%
%E2 Concentration
E2.Time = t;
E2.Y = y0E;

%P4 Concentration
P4.Time = t;
P4.Y = y0(end-50);

%FSH Concentration
FSH.Time = t;
FSH.Y = y0(end-46);

%LH Concentration
LH.Time = t;
LH.Y = y0(end-48);

%activ GnRH complex
GnRHcomp.Time = t;
GnRHcomp.Y = y0(end-1);

%LH release
LHRez.Time = t;
LHRez.Y = y0(end-17);

%FSHRez
FSHRez.Time = t;
FSHRez.Y = y0(end-13);
FSHfRez.Y = y0(end-14);
FSHdRez.Y = y0(end-12);

%GnRHRez
GnRHRez.Time = t;
GnRHRez.Y = y0(end-3);

%FSH analogon im central compartment
FSHAnaC.Time = t;
FSHAnaC.Y = y0(end-5);

%FSH analogon im dosing compartment
LHAnaC.Time = t;
LHAnaC.Y = y0(end-20);

%Yall
solutions.Time = t;
solutions.Y=y0(2:end)';
%
%-----------------------------------------------------------------------
%
while (t<te) 
    
    t
    
    fshAll = y0(end-46)+y0(end-5);
    fshimp = fshAll^3/(fshAll^3 + 12.35^3);
    timevec=poissonproc(paraPoi(1)+4*paraPoi(1)*fshimp,[t,te]);
    
    %set integration period for the current follicle
     if (~isempty(timevec))
         NextStart=timevec(1);
         tspan=[t,NextStart]; 
     else
        tspan=[t,te];
     end

    %determine number of follicles
    NumFollicles=size(y0,1)-para(2); 
    %set mass matrix for DAE system
    n=length(y0);
    M=eye(n,n);
    M(NumFollicles+2,NumFollicles+2)=0; %alg. eq. for E2
    M(NumFollicles+4,NumFollicles+4)=0; %alg. eq. for P4
    M(NumFollicles+49,NumFollicles+49)=0; %alg. eq. for FSH Med
    M(NumFollicles+34,NumFollicles+34)=0; %alg. eq. for LH Med
    
    %event function stops the integration, when ever an ovulation takes
    %place within the intervall tspan
    options = odeset('Mass',M,'events',@(t,y)evaluate_follicle(t,y,paraOde(7),para,count15mm,Par,LH));
    
    %search for dosing events in [tspan(1),tspan(2)]:
    dosing_timeIdx=[];
    if ~isempty(dosing_events1)
        dosing_timeIdx=intersect(find(dosing_events1(1,:)>tspan(1)),find(dosing_events1(1,:)<=tspan(2)));
    end
    
    %solve differential equations
    para(1) = 0; %call testfunction to test
    Y=[];
    T=[];
    if ~isempty(dosing_timeIdx)
        tstart=tspan(1);
        tend=tspan(2);
        yInitial=y0;
        for i=1:length(dosing_timeIdx)
            tspan=[tstart;dosing_events1(1,dosing_timeIdx(i))];
            if tspan(1)~=tspan(2)
                [ ti, yi ] = ode15s(@(t,y)testfun_LutStim(t,y,Tovu,Follicles,para,paraOde,Par,dd1), tspan, yInitial, options);
                T=[T;ti(2:end)];
                Y=[Y;yi(2:end,:)];
                tstart=T(end);
                yInitial=Y(end,:);
            end
            dd1 = dosing_events1(2,dosing_timeIdx);
        end
        tspan=[T(end);tend];
         [ ti, yi ] = ode15s(@(t,y)testfun_LutStim(t,y,Tovu,Follicles,para,paraOde,Par,dd1), tspan, yInitial, options);
         T=[T;ti(2:end)];
         Y=[Y;yi(2:end,:)];
    else
        [T,Y] = ode15s(@(t,y)testfun_LutStim(t,y,Tovu,Follicles,para,paraOde,Par,dd1),tspan,y0,options);
    end
    
    %find the max. foll. diameter in each time step and the time

    [m,n]=size(Y(:,1:(end-para(2))));
    if n<2
        M=Y(:,1:(end-para(2)));
    else
        M=max(Y(:,1:(end-para(2)))');
        M=M';
    end
    
    FollMaxTime = [M T];
    
    count10mm=0;
    for i= 1:n 
        if ((Y(1,i)>= 10) && (Y(1,i) < 14))
            count10mm = count10mm+1;
        end
    end
    helpFoll10m=[count10mm t];
    Foll10mm=[Foll10mm; helpFoll10m];
    
    count15mm=0;
    for i= 1:n 
        if Y(1,i)>= 14
            count15mm = count15mm+1;
        end
    end
    helpFoll15m=[count15mm t];
    Foll15mm=[Foll15mm; helpFoll15m];
    
    %read values of follicles that were active during last run
    %(but don't save initial value->already saved previously)
    SumYFoll = zeros(max(size(T(2:end))),1);
    
    for i = 1:Follicles.NumActive
        %saves all times of the foll that was active during last run
        Follicles.Follicle{Follicles.Active(i)}.Time = [Follicles.Follicle{Follicles.Active(i)}.Time; T(2:end)];
        %saves all sizes of the foll that was active during last run
        Follicles.Follicle{Follicles.Active(i)}.Y = [Follicles.Follicle{Follicles.Active(i)}.Y; Y(2:end,i)];
        %dummy for area calculation
        SumYFoll = SumYFoll + Y(2:end,i); 
    end
    
    %Werte für die Medikamentengabe setzen
    if Par(599)== Tovu && Par(596) == 0 
        for i = 1:Follicles.NumActive   
            if Follicles.Follicle{Follicles.Active(i)}.Destiny == -1
            %matrix of dosing times and drugs added: row1: times, row2: drug, row 3, dose
                if t-Par(599) >= 1 && t-Par(599) < 4
                    Par(599) = ceil(t);   
                    Par(598)=Par(599)+11;  
                %Menopur
                    numDoses=Par(598)-Par(599)+1;
                    dosing_events1=[[Par(599):Par(598)];[1:numDoses]];
                    Par(596) = 1;
                end
            end
        end
    end
    
    %saves last Y values of the foll. from ode45
    LastYValues = Y(end,1:end)';
    
    %save values for E2 
    E2.Time = [E2.Time; T(2:end)];
    E2.Y = [E2.Y; Y(2:end,end-52)];
    %save values for P4
    P4.Time = [P4.Time; T(2:end)];
    P4.Y = [P4.Y; Y(2:end,end-50)];
    %save values for LH
    LH.Time = [LH.Time; T(2:end)];
    LH.Y = [LH.Y; Y(2:end,end-48)];
    %save values for FSH 
    FSH.Time = [FSH.Time; T(2:end)];
    FSH.Y = [FSH.Y; Y(2:end,end-46)];
    %save values for active GnRH complex
    GnRHcomp.Time = [GnRHcomp.Time; T(2:end)];
    GnRHcomp.Y = [GnRHcomp.Y; Y(2:end,end-1)];
    %save values for LH release
    LHRez.Time = [LHRez.Time; T(2:end)];
    LHRez.Y = [LHRez.Y; Y(2:end,end-17)];
    %FSHRez
    FSHRez.Time = [FSHRez.Time; T(2:end)];
    FSHRez.Y = [FSHRez.Y; Y(2:end,end-13)];
    FSHfRez.Y =[FSHfRez.Y; Y(2:end,end-14)];
    FSHdRez.Y =[FSHdRez.Y; Y(2:end,end-12)];
    %FSH analogon im central compartment
    FSHAnaC.Time = [FSHAnaC.Time; T(2:end)];
    FSHAnaC.Y = [FSHAnaC.Y; Y(2:end,end-5)];
    %save solutions
    solutions.Time = [solutions.Time; T(2:end)];
    solutions.Y = [solutions.Y; Y(2:end,NumFollicles+1:end)];
    
    %no ovulation (no event) occured
    if T(end)==tspan(2)
        %initialize new follicle
        %Set initial values for new follicle
        Follicle1.Y = y0Foll;
        if( ~isempty(FSHVec) )
            Follicle1.FSHSensitivity = FSHVec(FollCounter);
            FollCounter = FollCounter + 1;
        end
        FSHSSave = Follicles.ActiveFSHS;
        Follicles.ActiveFSHS = [Follicles.ActiveFSHS Follicle1.FSHSensitivity];
        %Test if Follicle(s) could survive
        %(slope of growth-function positive or negative)
        testyvalues = LastYValues(1:(end-para(2)));
        testyvalues = [testyvalues; Follicle1.Y; LastYValues(end+1-para(2):end)];
        para(1) = 1;
        testyslope = testfun_LutStim(T(end),testyvalues,Tovu,Follicles,para,paraOde,Par,dd1);
        %if follicle got chance to survive->initiate new follicle and update
        %follicles-vector
        if( testyslope(end-para(2)) > 0 )
            Follicle1.Time = [ T(end) ];
            Follicle1.TimeDecrease = 0;
            Follicle1.Destiny = -1;
            Follicles.Number = Follicles.Number + 1;
            Follicle1.Number = Follicles.Number;
            Follicles.NumActive = Follicles.NumActive + 1;
            Follicles.Active = [Follicles.Active Follicle1.Number];
            Follicles.Follicle = {Follicles.Follicle{1:end} Follicle1};
            NewFollicle = [NewFollicle T(end)];
            LastYValues = testyvalues;
        else
            %no chance to survive save in NoNewFollicle for statistic
            %and set back old FSH-Sensitivity vector
            Follicles.ActiveFSHS = FSHSSave;
            NoNewFollicle = [NoNewFollicle T(end)];
        end
        t=T(end);
        TimeCounter=TimeCounter+1;
        if TimeCounter<=length(StartTimes)
            NextStart=StartTimes(TimeCounter);
        else
            NextStart=te;
        end
     
    else %ovulation occured
        t=T(end);
    end
    
    %check on every stop of the integrator if status of follicles changed
    
    %helping variables
    ActiveHelp = [];
    %determine actual slope of growth of follicles
    para(1) = 0;
    res = testfun_LutStim(T(end),LastYValues,Tovu,Follicles,para,paraOde,Par,dd1);
    %reset vector of active FSH sensitivities
    Follicles.ActiveFSHS = [];  
 
    %loop over all active follicles to set new destiny (1 ovulation,
    %-1 unknown/growth, -2 decreasing but not yet dead, 
    %3 its big but there is not enough LH for ovulation)  
    for i = 1:Follicles.NumActive
        %Save y-values of i-th (current) follicle
        yCurFoll = LastYValues(i);
        %slope is negative so the follicle is dereasing in size
        if( res(i) <= 0 ) 
             Follicles.Follicle{Follicles.Active(i)}.Destiny = -2;
        end 
        
        %follicle is big, but doesn't ovulate yet because there is not enough LH
        if(yCurFoll >= (paraOde(7))) && (Y(end,end-48) < para(9) && ...
                Follicles.Follicle{Follicles.Active(i)}.Destiny == -1)
               Follicles.Follicle{Follicles.Active(i)}.Destiny = 3;
               Follicles.Follicle{Follicles.Active(i)}.TimeDecrease=t;
        end

        if Follicles.Follicle{Follicles.Active(i)}.Destiny == 3 && ...
                (t-Follicles.Follicle{Follicles.Active(i)}.TimeDecrease)>=para(10)
            Follicles.Follicle{Follicles.Active(i)}.Destiny =-2;
        end
        
        %LH hoch genug ist 
        if Y(end,end-48) >= para(9)
            if (yCurFoll >= paraOde(7)) && (Follicles.Follicle{Follicles.Active(i)}.Destiny==-1) ||...
               (yCurFoll >= paraOde(7)) && (Follicles.Follicle{Follicles.Active(i)}.Destiny==3)
                th = t-0.5;
                [val, idx] = min(abs(LH.Time-th)); 
               if (LH.Y(idx)) >= 20 
                    Follicles.Follicle{Follicles.Active(i)}.Destiny = 4;
                    Follicles.Follicle{Follicles.Active(i)}.TimeDecrease=t;
               end 
            end
        end
            
        %Follicle ovulates
        if (Follicles.Follicle{Follicles.Active(i)}.Destiny == 4) && ...
                Follicles.Follicle{Follicles.Active(i)}.TimeDecrease + 0.5 <= t
                Follicles.Follicle{Follicles.Active(i)}.Destiny = 1;
                Tovu=T(end);
                if Tovu > 80 && Par(596) == 0 
                    Par(599) = Tovu;
                end
        end
        
        if(Follicles.Follicle{Follicles.Active(i)}.Destiny ~= 1)
            %put the follicle back to the list of actives and its FSH
            ActiveHelp = [ActiveHelp Follicles.Active(i)];
            %sensitivity back in the FSH vector...
            Follicles.ActiveFSHS = [Follicles.ActiveFSHS Follicles.Follicle{Follicles.Active(i)}.FSHSensitivity];
        end
        
    end
    %Update list of active follicles
    Follicles.Active = ActiveHelp;
    %find out how many follicles are active...
    Follicles.NumActive = size(ActiveHelp,2);
    
    %determine new initial values for all differential equations to
    %solve
    y0old = [];
    for i = 1:(Follicles.NumActive)
        y0old = [y0old Follicles.Follicle{Follicles.Active(i)}.Y(end)];
    end
    y0old = y0old';
    y0 = [y0old;LastYValues(end+1-para(2):end)];
    %integration end reached -> end integration for real---
    t = T(end);
    if( te - t < 0.001 )
        t = te;
    end
    
    LastYValues = Y(end,1:end)';
    count = 0;
    if Par(596) == 1 
        for i = 1:Follicles.NumActive
            %Save y-values of i-th (current) follicle
            yCurFoll = LastYValues(i);
            if yCurFoll >= 18 
                count = count + 1; 
            end
        end
    end
    
    if Par(596) == 1 && t > Par(598)+1 || Par(596) == 1 && count >= 3
        break
    end
    
end

%vector to save informations about the ovulating follicle
%[number of follicle; starttime; time of ovulation; lifetime]
FollOvulInfo=[];
%vectors for plotting the FSH and Progestoron curves
Time = [0:0.1:te];
widthofline = 2;

%plotting
if( ShowStuff )
    hf=figure(1);
    clf;
    hold on;
end

% save t_start t_end destiny of all follicles
FollInfo = [];
%growth of the follicles
FollInfo2 = [];

for i = 1:Follicles.Number
    %fill follicle information variable...
    help = [Follicles.Follicle{i}.Time(1); Follicles.Follicle{i}.Time(end); Follicles.Follicle{i}.Destiny;i];
    FollInfo = [FollInfo help];
    
    FollInfo2 = [Follicles.Follicle{i}.Time Follicles.Follicle{i}.Y];
    
    if (SaveStuff)
        FileName = sprintf('Follicle%d.csv',i);
        fullFileName = fullfile(DirStuff, FileName);
        csvwrite(fullFileName,FollInfo2)
    end
    
    if(Follicles.Follicle{i}.Destiny==1)
        helpFOT=[i;Follicles.Follicle{i}.Time(1);Follicles.Follicle{i}.Time(end);...
            Follicles.Follicle{i}.Time(end)-Follicles.Follicle{i}.Time(1)];
        FollOvulInfo=[FollOvulInfo helpFOT];
    end    
    
    if(ShowStuff)
        h = plot(Follicles.Follicle{i}.Time,Follicles.Follicle{i}.Y,'Color',[0 0 0],...
                 'DisplayName','x1','LineWidth', widthofline);
    end
end

if(ShowStuff)
   %fsh
    hfsh = plot(FSH.Time,FSH.Y,'Color',[1/2 1 1/2],...
         'DisplayName','x1','LineWidth', widthofline);
     
   %LH  
   hLH = plot(LH.Time,LH.Y,'Color',[1 1/4 1/2],...
         'DisplayName','x1','LineWidth', widthofline);
  
   %P4
    hp4 = plot(P4.Time,P4.Y,'Color',[1 0 1],...
             'DisplayName','x1','LineWidth', widthofline);
    
    %threshold when you can measure the follicle size
    hTwo=plot(xlim,[4 4],'Color','r');
    
    %plot for the follicle size, ammount of FSH and ammount of P4 
    h=[h hfsh hTwo hp4 hLH]; %hOvul];
    xlabel('time in d','fontsize',15);
    ylabel('follicle diameter in mm','fontsize',15);
    %fontsize of plot ticks
    ax = gca;
    set(ax, 'Box', 'off' );
    ax.FontSize = 15; 
    set(gca,'linewidth',1.5);
    legend(h,{'follicle growth','FSH','measurable','P4', 'LH'},'fontsize',15,...
        'Location','NorthEast');%,'ovulation');

%     if SaveStuff
%         FileName = sprintf('FolliclePlot');
%         fullFileName = fullfile(DirStuff, FileName);
%         saveas(gca, fullFileName, 'png');
%         saveas(gca, fullFileName, 'fig');
%     end
    
    
    %plot for the numberof follicles >=4mm an the max. size at the timestep
    figure(2);
    clf;
    %yyaxis left
    plot(Foll10mm(:,2),Foll10mm(:,1),Foll15mm(:,2),Foll15mm(:,1))
    title('Foll. size','fontsize',24)
    ylabel('No. of follicles','fontsize',24)    
    xlabel('time in d','fontsize',24)
    legend({'10mm <= #follicle < 14mm', '#follicle >= 14mm'},'fontsize',24,'Location','NorthEast')
        
%     if SaveStuff
%         FileName = sprintf('FollSizePlot');
%         fullFileName = fullfile(DirStuff, FileName);
%         saveas(gca, fullFileName, 'png');
%         saveas(gca, fullFileName, 'fig');
%     end
    
    figure(3);
    plot(P4.Time,P4.Y, FSH.Time,FSH.Y,'LineWidth',2);
    set(gca,'fontsize',24);
    legend({'P4','FSH'},'fontsize',24,...
        'Location','NorthEastOutside');%,'ovulation');
    
%     if SaveStuff
%         FileName = sprintf('FSH_P4');
%         fullFileName = fullfile(DirStuff, FileName);
%         saveas(gca, fullFileName, 'png');
%         saveas(gca, fullFileName, 'fig');
%     end
    
    figure(4);
    plot(E2.Time,E2.Y, LH.Time, LH.Y, 'LineWidth',2);
    set(gca,'fontsize',24);
    legend({'E2', 'LH'},'fontsize',24,...
        'Location','NorthEastOutside');%,'ovulation');
    
%     if SaveStuff
%         FileName = sprintf('E2_LH');
%         fullFileName = fullfile(DirStuff, FileName);
%         saveas(gca, fullFileName, 'png');
%         saveas(gca, fullFileName, 'fig');
%     end
  
    sumFSH = FSH.Y + FSHAnaC.Y;
    figure(7);
    plot(FSH.Time,FSH.Y,FSHAnaC.Time,FSHAnaC.Y,FSH.Time,sumFSH, 'LineWidth',2);
    set(gca,'fontsize',24);
    legend({'FSH','FSHc','Sum FSH'},'fontsize',24,...
        'Location','NorthEastOutside');

%     figure(6);
%     plot(test.Time,test.Y, 'LineWidth',2);
%     set(gca,'fontsize',24);
%     legend({'LHpit'},'fontsize',24,...
%         'Location','NorthEastOutside');
    

    NewFollicle;
    NoNewFollicle;
    FollOvulInfo
    
    %indexes of solutions are number from Model28_ODE + 1
    solutions = [solutions.Time solutions.Y];
    
    if (SaveStuff)
        FileName = sprintf('NumericalSolutions_allComponants.csv');
        fullFileName = fullfile(DirStuff, FileName);
        csvwrite(fullFileName,solutions)
    end
    
    %Cycle length
    OvuT = FollOvulInfo(3,:);
    Cyclelength = diff(OvuT)
    Cyclelengthmean = mean(Cyclelength)
    NumCycles = length(Cyclelength)

    FollperCycle =[];
    for i = 1:NumCycles
        t1 = OvuT(i);
        t2 = OvuT(i+1);
        count = 0;
        tp = length(FollInfo(1,:));
        for j = 1:tp
            if FollInfo(1,j) > t1 && FollInfo(1,j) < t2
                count = count +1;
            end
        end
        FollperCycle =[FollperCycle count];
    end
    
    FollperCyclemean = mean(FollperCycle)
    
    a = sum(FollperCycle);
    rest = n - a;
    
    CycleInfo = [[0 Cyclelength]; [rest FollperCycle]; OvuT];
    
%     if (SaveStuff)
%         FileName = sprintf('CycleInfo.csv');
%         fullFileName = fullfile(DirStuff, FileName);
%         csvwrite(fullFileName,CycleInfo)
%     end
        
    [val, idx] = min(abs(E2.Time-(Par(599)-1)));
    E2dm1 = E2.Y(idx);
    [val, idx] = min(abs(E2.Time-(Par(599)+1)));
    E2d1 = E2.Y(idx);
    [val, idx] = min(abs(E2.Time-(Par(599)+5)));
    E2d6 = E2.Y(idx);
    E2dend= E2.Y(end);

    [val, idx] = min(abs(P4.Time-(Par(599)-1)));
    P4dm1 = P4.Y(idx);
    [val, idx] = min(abs(P4.Time-(Par(599)+1)));
    P4d1 = P4.Y(idx);
    [val, idx] = min(abs(P4.Time-(Par(599)+5)));
    P4d6 = P4.Y(idx);
    P4dend= P4.Y(end);
    
    [val, idx] = min(abs(LH.Time-(Par(599)-1)));
    LHdm1 = LH.Y(idx);
    [val, idx] = min(abs(LH.Time-(Par(599)+1)));
    LHd1 = LH.Y(idx);
    [val, idx] = min(abs(LH.Time-(Par(599)+5)));
    LHd6 = LH.Y(idx);
    LHdend= LH.Y(end);
    
    [val, idx] = min(abs(FSH.Time-(Par(599)-1)));
    FSHdm1 = sumFSH(idx);
    [val, idx] = min(abs(FSH.Time-(Par(599)+1)));
    FSHd1 = sumFSH(idx);
    [val, idx] = min(abs(FSH.Time-(Par(599)+5)));
    FSHd6 = sumFSH(idx);
    FSHdend= sumFSH(end);

    MedInfo = [count10mm count15mm E2dm1 E2d1 E2d6 E2dend P4dm1 P4d1 P4d6 P4dend ...
                LHdm1 LHd1 LHd6 LHdend FSHdm1 FSHd1 FSHd6 FSHdend dosing_events1(1,1) t];


    if Par(596) == 1 && (SaveStuff)
        FileName = sprintf('MedInfo.csv');
        fullFileName = fullfile(DirStuff, FileName);
        M = load(fullFileName);
        M = [M; MedInfo];
        csvwrite(fullFileName,M);
    end

end

end