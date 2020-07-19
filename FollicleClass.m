classdef FollicleClass < handle

    properties
        
        Number; % number of active foll
        Follicle; % number of created foll
        %%(Number, Time, Size, Density (1- not clear jet, 0 - foll died, 
        %%1 - foll ovulated, 2 - foll is decreasing in size), FSHSensetivity, TimeDecrease)
        Active; % number of currently active foll
        NumActive; % all foll structures 
        ActiveFSHS; %all foll FSH sensitivities
        
    end

    methods
        function obj = FollicleClass(y0Follicle,FSHSensitivity,t)
     
            Foll.Number = 1;
            Foll.Time = [t];
            Foll.Y = y0Follicle;
            %the foll destinys: -1: not clear yet, 0:died, 1:ovulated, 
            %-2:decreasing in size but did not yet died 
            Foll.Destiny = -1;
            Foll.FSHSensitivity = FSHSensitivity;
            %time the size of the foll started to decrease
            Foll.TimeDecrease = 0;
            
            obj.Active = [1];
            obj.Number = 1;
            obj.NumActive = 1;
            obj.Follicle = {Foll};
            obj.ActiveFSHS = [Foll.FSHSensitivity];
        end
        
    end
end