function dy = ODE_Model_NormalCycle(t,y,Par)
%
%-----------------------------------------------------------------------
%
[r,c] = size(y);
dy=zeros(r,c);
%
%-----------------------------------------------------------------------
%
i_RP_LH    =  r-49;
i_LH       =  r-48;
i_RP_FSH   =  r-47;
i_FSH      =  r-46;
i_P4       =  r-50;

i_FSHfoll  =  r-45;  

i_RLH      =  r-18;

i_RFSH     =  r-14;
i_RFSH_des =  r-12;
i_FSHR     =  r-13;

i_GnRH     =  r-4;
i_RecGa    =  r-3;
i_RecGi    =  r-2;
i_GReca    =  r-1;
i_GReci    =  r;  

y_e2 = y(r-52);
y_p4 = y(r-50);
y_lh  = y(r-48);
y_fsh = y(r-46);
%
%-----------------------------------------------------------------------
%
%GnRH frequency and mass 
%
  yGfreq = 16.d0 / ( 1.0d0 + ( y_p4 / Par(203) ) ^ Par(204) ) ...
              * ( 1.0d0 + y_e2 ^ Par(206) ...
              / ( Par(205) ^ Par(206) + y_e2 ^ Par(206) ) );
%
  yGmass = y_e2 ^ Par(209) ...
              / ( Par(208) ^ Par(209) + y_e2 ^ Par(209) ) ...
              + Par(210) ^ Par(211) ...
              / ( Par(210) ^ Par(211) + y_e2 ^ Par(211) );
%
%-----------------------------------------------------------------------
%
% LH in pituitary
%
  hp_e2 = ( y_e2 / Par(3) ) ^ Par(6) ... 
          / ( 1.0d0 + ( y_e2 / Par(3) ) ^ Par(6) );
            
  f_LH_prod1 = Par(1) + Par(2) * hp_e2; 
  f_LH_prod2 = 1.0d0 + ( y_p4 / Par(4) ) ^ Par(7); 
  f_LH_prod  = f_LH_prod1 / f_LH_prod2; 
%
  f_LH_rel = ( Par(16) + Par(5) * ... 
               ( ( y(i_GReca)) / Par(8) ) ^ Par(9) ...
             / ( 1.0d0 + ( ( y(i_GReca) ) ...
               / Par(8) ) ^ Par(9) ) ); %...
%
  dy(i_RP_LH) = f_LH_prod - f_LH_rel * y(i_RP_LH);
%
%-----------------------------------------------------------------------
%
% LH in the blood
%     
  dy(i_LH) =   ( 1.0d0 / Par(12)) * f_LH_rel * y(i_RP_LH) ...
             - Par(230) * y_lh * y(i_RLH) ...
             - Par(231) * y_lh;
%
%-----------------------------------------------------------------------
%
% FSH in pituitary
%
  f_FSH_prod1 = Par(21);
%
  f_FSH_prod2 = 1.d0 + (y(i_P4)/5) ^ 2;
%
  hm_freq = 1.0d0 / ( 1.0d0 + ( yGfreq / Par(11) ) ^ Par(13) );
%
  f_FSH_prod = f_FSH_prod1 / f_FSH_prod2 * hm_freq;
%
  f_FSH_rel =   Par(17) + Par(28) ...
              * ( ( y(i_GReca) ) ... 
                / Par(18) ) ^ Par(20) ...
              / ( 1.0d0 + ( ( y(i_GReca) ) ... 
                / Par(18) ) ^ Par(20) ); 
%     
   dy(i_RP_FSH) = f_FSH_prod - f_FSH_rel * y(i_RP_FSH);
%
%-----------------------------------------------------------------------
%
%  FSH in the blood
%      
   dy(i_FSH)  =   1.0d0 / Par(12) * f_FSH_rel * y(i_RP_FSH) ...  
               - Par(243) * y_fsh ...
               - Par(241) * y_fsh ;
           
%
%-----------------------------------------------------------------------
%
%  FSH in ovaries
%
    dy(i_FSHfoll)  =  Par(243) * (y_fsh) * Par(246)...
               - Par(240) * y(i_FSHfoll) * y(i_RFSH) ...
               - Par(245) * (i_FSHfoll);
%
%-----------------------------------------------------------------------
%
% FSH receptor mechanisms
%
% FSH free receptors  
%
   dy(i_RFSH) =   Par(242) * y(i_RFSH_des) ...
                - Par(240) * y(i_FSHfoll) * y(i_RFSH);
%
% bound FSH receptors
%
  dy(i_FSHR) = Par(240) * y(i_FSHfoll) * y(i_RFSH) - Par(244) * y(i_FSHR);
%
% desensitized FSH receptors
%
  dy(i_RFSH_des) = Par(244) * y(i_FSHR) - Par(242) * y(i_RFSH_des); 
%
%-----------------------------------------------------------------------
%
% GnRH in pituitary 
%
  dy(i_GnRH) =   Par(301) * yGmass * yGfreq ... 
               - Par(302) * y(i_GnRH) * y(i_RecGa) ...
               + Par(303) * y(i_GReca) ...
               - Par(300) * y(i_GnRH);
%
%-----------------------------------------------------------------------
%
% GnRH receptor mechanisms
%
% active GnRH receptor
%
  dy(i_RecGa) =   Par(303) * y(i_GReca) ...
                - Par(302) * y(i_GnRH) * y(i_RecGa) ...  
                - Par(306) * y(i_RecGa) ...
                + Par(307) * y(i_RecGi);
%
% inactive GnRH receptor	 
%
  dy(i_RecGi) =   Par(311) ...
                + Par(306) * y(i_RecGa) ...
                - Par(307) * y(i_RecGi) ...
                + Par(305) * y(i_GReci) ...
                - Par(308) * y(i_RecGi);
%
% active GnRH-receptor complex	 
%
  dy(i_GReca) =   Par(302) * y(i_GnRH) * y(i_RecGa) ...
                - Par(303) * y(i_GReca) ...
                - Par(309) * y(i_GReca) ...
                + Par(310) * y(i_GReci);
%
% inactive GnRH-receptor complex
%
  dy(i_GReci) =   Par(309) * y(i_GReca) ... 
                - Par(310) * y(i_GReci) ...
                - Par(305) * y(i_GReci) ...
                - Par(304) * y(i_GReci);
%
%-----------------------------------------------------------------------
%
end