%%%Task 1%%%
clc;
close all;
%calculate Resistance
Re=input('Enter conductor resistivity :\n ');        %calc resistivity
len=input('Enter conductor length (Km) :\n ');       %Entering Length in Km

   % Making sure the length is between 0 and 250 km
    while len==0 || len > 250
        
        fprintf('\n[Invalid input]: Please renter your choice\n');
        len=input('Enter conductor length (Km) :\n '); 
        
    end
    
Dcm=input('Enter conductor diameter (cm) :\n ');     %calc diameter in cm
Dm=Dcm/100;
A=pi*((Dm/2).^(2));                            %cross sectionalarea of theconductor in m2
R=(Re*len)/A;         %resistance dc
R2=R*1.1;            %resistance ac

while 1
    system=input('\n Select the desired system :\n\n -> 1  symmetrical \n -> 2  unsymmetrical\n'); 
 
    switch system
        
      case 1 %symmetrical
 
        while 1
            bundle=input('\n select number bundle per phase :\n ->1 one \n->2 two \n->3 Three \n->4 Four\n');
          
            switch bundle
              
            case 1
               %one bundle
              D1=input('\nEnter D of three phase line :\n ');
              r=Dm/2;
              GMR1=r*exp(-0.25);
              Lphase=(2*10^(-7))*log(D1/GMR1);
              Cphase=(2*pi*8.85*10^(-12))/(log(D1/r));
              break;
              
            case 2
                %Two bundle
               D2=input('\nEnter D of three phase line :\n');
               d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR2=sqrt(d*re);
               Lphase=(2*10^(-7))*log(D2/GMR2);
               Cphase=(2*pi*8.85*10^(-12))/(log(D2/sqrt(d*r)));
               break;
               
             case 3
               %Three bundle
               
               D3=input('\nEnter D of three phase line :\n');
               d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR3=nthroot(d*d*re,3);
               Lphase=(2*10^(-7))*log(D3/GMR3);
               Cphase=(2*pi*8.85*10^(-12))/(log(D3/nthroot(d*d*r,3)));
               break;
               
               case 4
                   %Four bundle
                   
               D4=input('\nEnter D of three phase line :\n');
               d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR4=nthroot(d*d*(d*sqrt(2))*re,4);
               Lphase=(2*10^(-7))*log(D4/GMR4);
               Cphase=(2*pi*8.85*10^(-12))/(log(D4/nthroot(d*d*(d*sqrt(2))*r,4)));
               break;
               
               otherwise
                fprintf('\n[Invalid input]: Please renter your choice\n');
                continue;  
            end
        end
        break;
        
    case 2 %not symmetrical
        
         while 1
            bundle=input('\n select number bundle per phase :\n ->1 one \n->2 two \n->3 Three \n->4 Four\n');
          
            switch bundle
                
              case 1                 
               %one bundle
              D1=input('\nEnter D1 of three phase line (m) :\n ');
              D2=input('\nEnter D2 of three phase line (m) :\n ');
              D3=input('\nEnter D3 of three phase line (m) :\n ');
              GMD1=nthroot(D1*D2*D3,3);
              r=Dm/2;
              GMR1=r*exp(-0.25);
              Lphase=(2*10^(-7))*log(GMD1/GMR1);
              Cphase=(2*pi*8.85*10^(-12))/(log(GMD1/r));
              break;
              
              case 2
                %Two bundle
              D1=input('\nEnter D1 of three phase line (m) :\n ');
              D2=input('\nEnter D2 of three phase line (m) :\n ');
              D3=input('\nEnter D3 of three phase line (m) :\n ');
              GMD2=nthroot(D1*D2*D3,3);
              d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR2=sqrt(d*re);
               Lphase=(2*10^(-7))*log(GMD2/GMR2);
               Cphase=(2*pi*8.85*10^(-12))/(log(GMD2/sqrt(d*r)));
               break;
               
               case 3
               %Three bundle 
              D1=input('\nEnter D1 of three phase line (m) :\n ');
              D2=input('\nEnter D2 of three phase line (m) :\n ');
              D3=input('\nEnter D3 of three phase line (m) :\n ');
              GMD3=nthroot(D1*D2*D3,3);
               d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR3=nthroot(d*d*re,3);
               Lphase=(2*10^(-7))*log(GMD3/GMR3);
               Cphase=(2*pi*8.85*10^(-12))/(log(GMD3/nthroot(d*d*r,3)));
               break;
               
               case 4
                   %Four bundle
                   
              D1=input('\nEnter D1 of three phase line (m) :\n ');
              D2=input('\nEnter D2 of three phase line (m) :\n ');
              D3=input('\nEnter D3 of three phase line (m) :\n ');
              GMD4=nthroot(D1*D2*D3,3);
               d=input('\nEnter d between bundles :\n');
               r=Dm/2;
               re=r*exp(-.25);
               GMR4=nthroot(d*d*(d*sqrt(2))*re,4);
               Lphase=(2*10^(-7))*log(GMD4/GMR4);
               Cphase=(2*pi*8.85*10^(-12))/(log(GMD4/nthroot(d*d*(d*sqrt(2))*r,4)));
               break;
               
                 otherwise
                fprintf('\n[Invalid input]: Please renter your choice\n');
                continue;
            end
         end
         break;
         
         otherwise
                fprintf('\n[Invalid input]: Please renter your choice\n');
                continue;
    end
end



f=input('Enter frequancy (Hz) :\n ');
L_tot = Lphase * len ; %calculating the total inductance
C_tot = Cphase * len ; %calculating the total capacitance
Z = R* len + 1i*2*pi*f*L_tot;
Y = 1i*2*pi*f*C_tot;


if len < 80
    fprintf('\n The Transmission line is Short Model\n');
    A= 1 ;
    B= Z ;
    C= 0 ;
    D= 1;
    
elseif 80 <= len && len<= 250 
    
   while 1 
     model=input('choose your model:\n (1)Pi model\n (2)T model\n ');
    switch model 
        case 1 
            fprintf('\n The Transmission line is pi Model\n');
            A= 1+ Z*Y./2;
            B= Z;
            C= Y.*(1+Z*Y./4);
            D= 1+ Z*Y./2;
           break;
            
        case 2
            fprintf('\n The Transmission line is T Model\n');
            A= 1+ Z*Y./2;
            B = Z.*(1+Z*Y./4);
            C = Y;
            D= 1+ Z*Y./2;
            break;
            
    end
   end
  
end
 
fprintf('\n Parameter A = %.4f + %.4fi',real(A),imag(A));
fprintf('\n Parameter B = %.4f + %.4fi',real(B),imag(B));
fprintf('\n Parameter C = %.4f + %.4fi',real(C),imag(C));
fprintf('\n Parameter D = %.4f + %.4fi',real(D),imag(D)); 

%Task3
% Taking Recieving line voltage
Vr=input('\nEnter recieving line voltage(KV): \n');Vr=Vr*1000;
% Asking the user to select either Case I or Case II
selectedCase = input('choose your case:\n (1)Case I\n (2)Case II:\n ');
% Making sure the user selected a valid case
while selectedCase ~= 1 && selectedCase ~= 2
    selectedCase = input('[INVALID INPUT],choose your case:\n (1)Case I\n (2)Case II:\n  ');
end

switch selectedCase
    case 1
           %%caseI
            % Calculating Phase recieving voltage
            Vr_phase= Vr/(3^(1/2));
            % Power factor
            pfr=0.8;
            % Initializing an array that contains the active recieving
            % power values from 0 to 100000 
            powerR_lag=0:100000;
            % Getting the magnitude of the recieving current
            Ir_mag= powerR_lag/(3*pfr*Vr_phase);
            % Calculating the recieving current magnitude and phase
            Ir=Ir_mag*(cosd(-36.87)+1i*sind(-36.87));
            % calculating the sending end voltage
            Vs= A* Vr_phase+B* Ir;
            % calculating the sending end current
            Is= C* Vr_phase+D* Ir;
            % calculating the no load recieving voltage
            Vr_noload= Vs/ A;
            % calculating the voltage regulation value
            VoltReg=(abs(Vr_noload)-abs(Vr_phase))/abs(Vr_phase);
            % Getting the angle by which the sending current lags or leads the
            % sending voltage
            theta=angle(Vs)-angle(Is);
            % Caculating the sending power
            PowerS =3.*abs(Vs).*abs(Is).*cos(theta);
            % Getting the efficiency
            eff=(powerR_lag./PowerS)*100;
            Active_power=real(powerR_lag);

            cftool
figure;
plot(Active_power,eff);
xlabel('Active power (Watt)')
ylabel('Efficiency')
title('The Relation between the Active Power and the Efficiency')
figure;
plot(Active_power,VR);
xlabel('Active power (Watt)')
ylabel('Voltage Regulation')
title('The Relation between the Active Power and the Voltage Regulation')
    case 2
         %case 2
            % Active power is constant
            Pr=100e3;
            % Power factor values in an array from 0.3 to 1
            pf=0.3:0.01:1;
            % Getting the phase recieving voltage
            Vr_phase= Vr/(3^(1/2));
            % getting the magnitude of the recieving end current
            Ir_mag = Pr./(3*pf*Vr_phase);
            %for lagging pf 
            lag_VR = 0;
            lag_eff = 0;
            % solving for both lead and lag
            for k=1:2
                % Getting the recieving end current
                Ir=Ir_mag.*(cos((-1)^(k)*acos(pf))+1i*sin((-1)^(k)*acos(pf)));
                % Getting the sending end voltage
                Vs= A* Vr_phase+B.* Ir;
                % Getting the sending end current
                Is= C* Vr_phase+D.* Ir;
                % Getting the no load recieving end voltage
                Vr_noload= Vs./ A;
                % getting the voltage regulation value
                VoltReg=((abs(Vr_noload)-abs(Vr_phase))/abs(Vr_phase))*100;
                % Getting the angle by which the sending current lags or leads the
                % sending voltage
                theta=angle(Vs)-angle(Is);
                % Caculating the sending power
                PowerS =3.*abs(Vs).*abs(Is).*cos(theta);
                % Getting the efficiency
                eff=(Pr./PowerS).*100;
                % Plotting
                cftool
                figure;
                subplot(2,1,1);
                plot(pf,eff);
                xlabel('Effeciency versus Power Factor');
                subplot(2,1,2);
                plot(pf,VoltReg);
                xlabel('Voltage Regulation versus Power Factor');
                if k==1
                    suptitle('Lagging power factor');
                    lag_VR = VoltReg;
                    lag_eff = eff;
                else
                    suptitle('Leading power factor');
                end
            end
end
        







  
       
