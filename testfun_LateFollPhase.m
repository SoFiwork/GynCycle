%ODE function to solve
function f=testfun_RandStim(t,y,Tovu,Follicles,para,paraOde,Par,dd1)%[f NumFoll MaxSize]

%determine number of active follicles
NumFollicles=size(y,1)-para(2);   

if( NumFollicles > 0 )
    x= y(1:NumFollicles);    
else
    x=0;
end

if(NumFollicles > 0  && para(1)==0)
    for i = 1:(NumFollicles) 
        if( (Follicles.Follicle{Follicles.Active(i)}.Destiny == -2 ))
            x(i)=0;
        end
    end
end
%
%-----------------------------------------------------------------------
%
dy = Model28_LutStim(t, y, Par);
f=dy;
%
%-----------------------------------------------------------------------
%
[r,c] = size(y);
%fshAll = y(r-46)+y(r-5); 
fshrezcomp = y(r-13);
p4all = y(r-50);
SumV = sum(x.^paraOde(1));

 for i = 1:(NumFollicles)
    
    %FSH sensitivity of the follicles
    fFSH=Follicles.ActiveFSHS(i);
    fsize = y(i);
    
    %growth rate
    gammaP4 = paraOde(2)*((1/(1+(p4all/3)^5))+(fshrezcomp^10/(0.59^10+fshrezcomp^10)));
    
    %negative Hill function for FSH with kappa(proportion of self harm)
    kappaFSH=paraOde(5)*(0.85^25/(0.85^25+fshrezcomp^25));

    xiP4=paraOde(3);
    ffsh = (fshrezcomp)^1/(fshrezcomp^1+(fFSH)^1);

    X =ffsh*(xiP4-y(i))*y(i)*(gammaP4-(kappaFSH*(SumV-(paraOde(4)*(y(i)^paraOde(1))))));
     if (para(1)== 1)
         if(X<=0)
             NoFoll=X;
         end
     end
             
     if( para(1) == 0 )                                                  
        %if the size of the foll. is decreasing (or constant), 
        %the size increases very slow and the follicle is 2 or more days
        %alive
        %or the foll. is big but alive for two or more days and has not
        %ovulated, then set destiny to decrease (-2)
        %and make it decreasing in size faster                                       
        if((X <= 0.01)|| ( Follicles.Follicle{Follicles.Active(i)}.Destiny == -2 )||...
                ((X<=0.1) && (( t - Follicles.Follicle{Follicles.Active(i)}.Time(1) ) >= 2 )) && ( Follicles.Follicle{Follicles.Active(i)}.Destiny ~= 4 ) ||... 
                ( Follicles.Follicle{Follicles.Active(i)}.Destiny == 3 && ((t- Follicles.Follicle{Follicles.Active(i)}.TimeDecrease)>=para(10))))
            if( Follicles.Follicle{Follicles.Active(i)}.Destiny ~= -2)
                %set time the follicle starts to decrease & set destiny to decrease
                  Follicles.Follicle{Follicles.Active(i)}.Destiny = -2;
                  Follicles.Follicle{Follicles.Active(i)}.TimeDecrease = t;
            end
            %to decrease the size of the follicle faster
            f(i)= -0.05*y(i)*(t-Follicles.Follicle{Follicles.Active(i)}.TimeDecrease);
          else
            %if not dying use normal equation
            f(i)=X;
        end
     else
        %if called to test use normal equation
        f(i)=X;
     end
     
 end
%
%-----------------------------------------------------------------------
%
if(NumFollicles > 0  && para(1)==0)
    for i = 1:(NumFollicles) 
        if( (Follicles.Follicle{Follicles.Active(i)}.Destiny == 4 ))
            x(i)=0;
        end
    end
end
SF = pi*sum((x.^5)./(x.^5+15^5).*(x.^2));
f(NumFollicles+2)=y(NumFollicles+2)-13*(1+0.02*SF)-150*exp(-0.06*(t-(Tovu+7))^2);

%Hemmung der Estradioproduktion -> keine Ahnung, wie...
%%Vielleicht concentration von Med berechnen und dann neg. Hill-Term
% if dd2 > 0 
%     f(NumFollicles+2)=(y(NumFollicles+2)-13*(1+0.02*SF)-150*exp(-0.06*(t-(Tovu+7))^2))*0.5;
% end
%
%-----------------------------------------------------------------------
%
%Calculation of P4 values
f(NumFollicles+4)=y(NumFollicles+4)-15*exp(-0.06*(t-(Tovu+7))^2);
%
%-----------------------------------------------------------------------
%
%Calculation of FSHAnaC
%
%-----------------------------------------------------------------------
%
if t >= Par(599) 

    n = dd1;
    H = 0;
    J = 0;
     
     if n < 6
         for i = 1:n
             dt = Par(599) + i - 1 ;
             h = ((Par(579)*(Par(578)^2))/((Par(578)-Par(577))^2)) ...
                        * [exp(-Par(578)*(t-dt)) * (Par(577)*(t-dt) ...
                        -Par(578)*(t-dt)-1)+exp(-Par(577)*(t-dt))];   
             H = H + h;
             
              dt = Par(599) + i - 1 ;
              j = ((Par(576)*(Par(575)^2))/((Par(575)-Par(574))^2)) ...
                    * [exp(-Par(575)*(t-dt)) * (Par(574)*(t-dt) ...
                    -Par(575)*(t-dt)-1)+exp(-Par(574)*(t-dt))];    
              J = J + j; 
         end
     else
         for i = 1:n
             dt = Par(599) + i - 1 ;
             h = ((Par(579)*(3/2)*(Par(578)^2))/((Par(578)-Par(577))^2)) ...
                        * [exp(-Par(578)*(t-dt)) * (Par(577)*(t-dt) ...
                        -Par(578)*(t-dt)-1)+exp(-Par(577)*(t-dt))];   
             H = H + h;
             
             dt = Par(599) + i - 1 ;
             j = ((Par(576)*(3/2)*(Par(575)^2))/((Par(575)-Par(574))^2)) ...
                    * [exp(-Par(575)*(t-dt)) * (Par(574)*(t-dt) ...
                    -Par(575)*(t-dt)-1)+exp(-Par(574)*(t-dt))];    
             J = J + j;
         end
     end
     f(NumFollicles+49)=y(NumFollicles+49)-H;
     f(NumFollicles+34)=y(NumFollicles+34)-J;
else 
    f(NumFollicles+49)=y(NumFollicles+49)-0;
    f(NumFollicles+34)=y(NumFollicles+34)-0;
 
end
%
%-----------------------------------------------------------------------
%
end