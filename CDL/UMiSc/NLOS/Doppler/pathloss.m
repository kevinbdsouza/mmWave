function [PL] = pathloss(d2D,Fc)

    % Set Parameters for path loss
    scenario = 'UMiSc'; 
    cVel = 3e8;
            
        
    %Apply scenario psecific pathloss and shadowing 
    if strcmp(scenario,'RMa')
        hBS = 35;
        hUT = 1.5;
        W = 20;
        h = 5;
        dBP = 2*pi*hBS*hUT*Fc/cVel;
        d3D = sqrt((hBS-hUT)^2 + d2D^2);
        
        PL1 = 20*log10(40*pi*d3D*Fc/3) + min(0.03* (h^1.72), 10)* log10(d3D) - min(0.044*h^1.72,14.77) + 0.002*log10(h)*d3D;
        PL2 = 20*log10(40*pi*dBP*Fc/3) + min(0.03* (h^1.72), 10)* log10(dBP) - min(0.044*h^1.72,14.77) + 0.002*log10(h)*dBP + 40*log10(d3D/dBP);
           
           if (d2D >= 10 && d2D <= dBP)
           PLLOS =  PL1;
           SFLOS = 4;  %#ok<*NASGU>
           elseif (d2D >= dBP && d2D <= 10000)
           PLLOS = PL2;  
           SFLOS = 6;
           end  
        
           PLNLOS = 161.04 - 7.1*log10(W) + 7.5*log10(h) - (24.37 - 3.7*(h/hBS)^2)*log10(hBS) + (43.42 - 3.1*log10(hBS))*(log10(d3D) - 3) ...
           + 20*log10(Fc) - (3.2*(log10(11.75*hUT))^2 - 4.97);
            
           
           if (d2D >= 10 && d2D <= 5000)
              PLNLOS = max(PLLOS,PLNLOS);
              SFNLOS = 8;
           end 
            
        
        
    elseif strcmp(scenario,'UMiSc') 
        hBS = 10;
        hUT = 1.5;
        hE = 1;
        dBPh = 4*(hBS - hE)*(hUT - hE)*Fc/cVel;
        d3D = sqrt((hBS-hUT)^2 + d2D^2);
        
        PL1 = 32.4 + 21*log10(d3D) + 20*log10(Fc);
        PL2 = 32.4 + 40*log10(d3D) + 20*log10(Fc) - 9.5*log10((dBPh^2) + (hBS - hUT)^2);
        
                   
           if (d2D >= 10 && d2D <= dBPh)
           PLLOS =  PL1;
           SFLOS = 4;
           elseif (d2D >= dBPh && d2D <= 5000)
           PLLOS = PL2;  
           SFLOS = 4;
           end 
            
        
           PLNLOS = 35.3*log10(d3D) + 22.4 + 21.3*log10(Fc) - 0.3*(hUT - 1.5);
                      
                     
            if (d2D >= 10 && d2D <= 5000)
                PLNLOS = max(PLLOS,PLNLOS);
                SFNLOS = 7.82;
            end 
            
            
        
    
    elseif strcmp(scenario,'UMa')
        hBS = 25;
        hUT = 1.5;
        if d2D <= 18
            g2D = 0;
        elseif d2D > 18
            g2D = (5/4)*(d2D/100)^3 * exp(-d2D/150);
        end 
        if hUT < 13
            Cd2DhUT = 0;
        elseif hUT >= 13 && hUT <= 23
            Cd2DhUT = ((hUT - 13)/10)^1.5 * g2D;
        end 
        
        if (1/(1+Cd2DhUT)) >= 0.5
            hE = 1;
        else 
            uni = [12,15,18,21];
            pos = randi(length(uni));
            hE = uni(pos);
        end 
        
        dBPh = 4*(hBS - hE)*(hUT - hE)*Fc/cVel;
        d3D = sqrt((hBS-hUT)^2 + d2D^2);
        PL1 = 32.4 + 20*log10(d3D) + 20*log10(Fc);
        PL2 = 32.4 + 40*log10(d3D) + 20*log10(Fc);
        
        
            if (d2D >= 10 && d2D <= dBPh)
            PLLOS =  PL1;
            SFLOS = 4;
            elseif (d2D >= dBPh && d2D <= 5000)
            PLLOS = PL2; 
            SFLOS = 4;
            end 
            
        
            PLNLOS = 13.54 + 39.08*log10(d3D) + 20*log10(Fc) - 0.6*(hUT - 1.5);
            
                       
            if (d2D >= 10 && d2D <= 5000)
                PLNLOS = max(PLLOS,PLNLOS);
                SFNLOS = 6;
            end 
            
         
        
        
    end 
    
    PL = PLNLOS + normrnd(0,SFNLOS);     % apply pathloss and shadowing
        
end 


        