function Complete_I_V_ramprate()
format shortG

%%%%%%%---Define parameters---%%%%%%%%
rel_permittivity =25*8.85*10^(-12);%[Nm^2/C^2]
boltzmann_constant = 1.381e-23;%[J/K]
%%%______Circuitary_______%%%
ramp_rate =[1e2 1e4 1e6];
%define source voltage
V_pos_amp = 1.25;%[V]
V_neg_amp = -1.75;%[V]
step_volt = 0.05;%[V]
%%%_____Electric Conductivity_____%%%
alpha_fil=-0.05;EC_HfO2x=5e3;%[kS/m]
alpha_gap=0.05;EC_gap=3e3;%[kS/m]
atomic_vibration=1e-13;%[s]
%%%_____Thermal Conductivity_____%%%
k_eff=10;
%%%______Chemical Energy______%%%
surface_tension = 0.01;%[J/m^2]
delta_mu_SET_J = 100e8;%[GJ/m^3]
delta_mu_RESET_J = 65e8;%[GJ/m^3]
beta_1 = 3.5e8;%[GJ/m^3]
beta_2 = 5e8;%[GJ/m^3]
delta_W_uc = 1;%eV
delta_W_i = 0.1;%eV
delta_W_mc = 0.3;%eV
%%%______filament Nucleation_____%%%
height_filament = 5e-9;%[m]
nucleation_barrier = 2.5*1.6e-19;%[J]
lambda = 6.6487; %following e-field of 1V/5nm
critical_nucleation_radius = 2.9e-9;%[m]
min_filament_radius = 0.6e-9;%[m]
alpha = min_filament_radius/critical_nucleation_radius;
%%%____create tables for database___%%%
DATA_OFF = table;
DATA_ON = table;
DATA_OFF_neg = table;
DATA_ON_neg = table;
I_V = table;
DATA_SET = table;
DATA_RESET = table;
DATA_min_SET = table;
DATA_min_RESET = table;

%define # of cycles
cycle = 2;

%set ramp rate (pulse length)
for ramp_rate_loop = 1 : length(ramp_rate)
    %loop for cycle
    %s_g and s_r are the manually picked stable gap length and filament radius for
    %the inital OFF state and doesnot follow thermodyanmics
    %trick is to run the program twice so that in the second cycle the
    %program starts with thermodynamically choosen gap and filament size
    s_g=2.7e-9;%[m]
    s_r=2.5e-9;%[m]
    for cycle_loop = 1 : cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%Positive Voltage%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%OFF_MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model_OFF= RRAM_OFF;
        model_OFF.hist.disable;
        %model_OFF.param.set('delta_E',random_E);
        model_OFF.param.set('V_amp', V_pos_amp);
        model_OFF.param.set('h_g',s_g);
        model_OFF.param.set('r_f',s_r);
        model_OFF.param.set('alpha_fil',alpha_fil);
        model_OFF.param.set('alpha_gap',alpha_gap);
        model_OFF.param.set('k_eff',k_eff);
        model_OFF.param.set('EC_HfO2x', EC_HfO2x);
        model_OFF.param.set('EC_gap', EC_gap);
        model_OFF.param.set('R',ramp_rate(1,ramp_rate_loop));
        model_OFF.param.set('t_rise',V_pos_amp/ramp_rate(1,ramp_rate_loop));
        model_OFF.study('std1').feature('time').set('tunit', 's');
        model_OFF.study('std1').feature('time').set('tlist', 'range(0,t_rise/30, 3/2*t_rise)');
        model_OFF.study('std1').run;
        time = 0 : (V_pos_amp/ramp_rate(1,ramp_rate_loop))/30 :3/2*V_pos_amp/ramp_rate(1,ramp_rate_loop);
        time_s = time(:);
        avg_gap_Temperature = mphmean(model_OFF, 'T', 2, 'selection',[3]);
        avg_gap_Temperature_K = avg_gap_Temperature(:);
        device_VoltageOFF_V = mphglobal(model_OFF, 'cir.IvsU1_v');
        device_CurrentOFF_uA = 1e6*abs(mphglobal(model_OFF, 'cir.IvsU1_i'));
        source_VoltageOFF_V = mphglobal(model_OFF, 'cir.V1_v');
        load_VoltageOFF_V = mphglobal(model_OFF, 'cir.R1_v');
        Threshold_Voltage_V = zeros(length(time),1);
        Voltage_wiggle_V = zeros(length(time),1);
        %find pulse length and temperaute dependent threshold voltage
        for j = 1 : length(time)
            Voltage_wiggle_V(j,1) = (height_filament*nucleation_barrier)/(boltzmann_constant*avg_gap_Temperature_K(j,1))*sqrt((3*pi^3*alpha^3*lambda*nucleation_barrier)/(32*rel_permittivity*critical_nucleation_radius^3));
            Threshold_Voltage_V(j,1) = Voltage_wiggle_V(j,1)/ log(Voltage_wiggle_V(j,1)/(ramp_rate(1,ramp_rate_loop)*atomic_vibration));
        end
        tbl_OFF = table(time_s,avg_gap_Temperature_K,Threshold_Voltage_V,device_VoltageOFF_V,device_CurrentOFF_uA,source_VoltageOFF_V,load_VoltageOFF_V);
        % when the device voltage exceeds the threshold voltage
        % the device switches
        count_c = 0;
        count_b = 0;
        Threshold_coloumn = zeros(length(time),1);
        for loop = 1: length(time)
            count_c = count_c+1;
            if Threshold_Voltage_V(loop)<device_VoltageOFF_V(loop)
                count_b = count_b+1;
                Threshold_coloumn(count_b) = count_c;
            end
        end
        % store device I-V upto threshold voltage then switch
        device_Voltage_V = device_VoltageOFF_V(1:Threshold_coloumn(1),1);
        device_current_uA = device_CurrentOFF_uA(1:Threshold_coloumn(1),1);
        I_V_OFF = table(device_Voltage_V,device_current_uA);
        DATA_OFF = vertcat(DATA_OFF,tbl_OFF);
        I_V = vertcat(I_V,I_V_OFF);
        
        %%%%%%%%%%%%%%%%%%%%%%%%Postitive Voltage%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%SET_MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        model_SET = RRAM_SET;
        model_SET.hist.disable;
        format shortG
        %Assign all the constant value
        height_filament = 5e-9;
        
        %Assign the range of source voltage with step increament strting from
        %threshold voltage
        I_V_SET1 = table;
        min_Current = 100e-6;
        max_Current = 230e-6;
        step_increament_Current = 10e-6;
        
        %Sweeps source voltage from minimum to maximum voltages in step increament
        for Current_Source_loop = min_Current:step_increament_Current:max_Current
            min_fil_radius =0.5e-9;max_fil_radius =9e-9; step_increament_radius = 2.5e-9;
            row =  min_fil_radius:step_increament_radius:max_fil_radius;
            FE_diff_SET = zeros(length(row),1);
            DATA3 = table;
            stable_radius = [];
            %%%%________Locate Minimium__________%%%%%%%%
            for radius_filament_loop = min_fil_radius:step_increament_radius:max_fil_radius
                r_f=radius_filament_loop; I_S= Current_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                    max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                    Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_SET(r_f,I_S,R);
                Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                radius_filament_m = r_f;Source_Current_A=Current_Source_loop;
                tbl=table(ramp_rate(1,ramp_rate_loop),Source_Current_A,source_VoltageSET_V,device_current_uA1,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,radius_filament_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                        max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Free_Energy_eV,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                        Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                DATA_SET = vertcat(DATA_SET,tbl);
                DATA3 = vertcat(DATA3,tbl);
                [num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(source_VoltageSET_V),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA1),'|',num2str(device_Voltage_V1),'|',num2str(radius_filament_m),'|',num2str(max_filament_Temperature_K),'|',num2str(delta_mu)]
                
            end
            for count_e = 1 : length(row)-1
                FE_diff_SET(count_e,1) = (DATA3.Free_Energy_eV(count_e+1)-DATA3.Free_Energy_eV(count_e));
            end
            %locates the minimum in free energy using the difference and records the row corresponding to it
            for count_f = 1 : length(row)-2
                if (FE_diff_SET(count_f)<0) && (FE_diff_SET(count_f+1)>0)
                    stable_radius = DATA3.radius_filament_m(count_f+1);
                end
            end
            %if minimum exist then uses Brents minimization for speedy
            %convergence
            if (isempty(stable_radius)==0)
                min_fil_radius =stable_radius - step_increament_radius;
                max_fil_radius =stable_radius + step_increament_radius;
                ITMAX=20; tol=1e-2; ZEPS=1e-11; CGOLD=0.3819660;
                ax=min_fil_radius; cx=max_fil_radius; bx=0.5*(min_fil_radius+max_fil_radius)+1e-9;
                a=ax; b=cx; v=bx;
                w=v; x=v; e=0;
                r_f=x; I_S=Current_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                    max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                    Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] =FE_SET(r_f,I_S,R);
                Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                fx=Free_Energy_eV; fv=fx; fw=fx;
                radius_filament_m = r_f ;
                Source_Current_A=Current_Source_loop;
                for iter = 1:ITMAX
                    [num2str(cycle_loop),'|',num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(Source_Current_A),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA1),'|',num2str(device_Voltage_V1),'|',num2str(radius_filament_m),'|',num2str(max_filament_Temperature_K),'|',num2str(avg_dielectric_Temperature_K),'|',num2str(delta_mu)]
                    tbl=table(ramp_rate(1,ramp_rate_loop),Source_Current_A,source_VoltageSET_V,device_current_uA1,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,radius_filament_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                        max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Free_Energy_eV,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                        Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                    DATA_SET = vertcat(DATA_SET,tbl);
                    xm = 0.5*(a+b); tol1=tol*abs(x)+ZEPS; tol2=2*tol1;
                    if (abs(x-xm)<=(tol2-0.5*(b-a))) %CONDITION 1
                        r_f=x; I_S=Current_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                        [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                            max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] =FE_SET(r_f,I_S,R);
                        Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                        radius_filament_m = r_f ;
                        tbl=table(ramp_rate(1,ramp_rate_loop),Source_Current_A,source_VoltageSET_V,device_current_uA1,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,radius_filament_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                            max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Free_Energy_eV,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                        DATA_SET = vertcat(DATA_SET,tbl);
                        [num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(Source_Current_A),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA1),'|',num2str(device_Voltage_V1),'|',num2str(radius_filament_m),'|',num2str(max_filament_Temperature_K),'|',num2str(avg_dielectric_Temperature_K),'|',num2str(delta_mu)]
                        break
                    end
                    if (abs(e)>tol1) %CONDITION 2
                        r=(x-w)*(fx-fw); q=(x-v)*(fx-fw); p=(x-v)*q-(x-w)*r; q = 2*(q-r);
                        if (q > 0)
                            p=-p;
                        end
                        q = abs(q); etemp=e; e=d;
                        if (abs(p)>=abs(0.5*q*etemp))||(p<=q*(a-x))||(p>=q*(b-x)) %CONDITON 2.1
                            if (x>=xm)
                                e=a-x;
                            else
                                e=b-x;
                            end
                            d=CGOLD*e;
                            if(abs(d)>=tol1)
                                u=x+d;
                            else
                                u=x+abs(tol1)*sign(d);
                            end
                            r_f=u; I_S=Current_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                            [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                                max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                                Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_SET(r_f,I_S,R);
                            Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                            fu=Free_Energy_eV; radius_filament_m=r_f;
                            if (fu<fx)
                                if(u>=x)
                                    a=x;
                                else
                                    b=x;
                                end
                                v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                            else
                                if (u<x)
                                    a=u;
                                else
                                    b=u;
                                end
                                if(fu<=fw)||(w==x)
                                    v=w; fv=fw; w=u; fw=fu;
                                elseif (fu<=fv)||(v==x)||(v==w)
                                    v=u; fv=fu;
                                end
                                continue
                            end
                        else
                            d=p/q; u=x+d;
                            if (u-a<tol2)||(b-u<tol2)
                                d=abs(tol1)*sign(xm-x);
                            end
                            if(abs(d)>=tol1)
                                u=x+d;
                            else
                                u=x+abs(tol1)*sign(d);
                            end
                            r_f=u; I_S=Current_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                            [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                                max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                                Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_SET(r_f,I_S,R);
                            Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                            fu=Free_Energy_eV; radius_filament_m =r_f;
                            if (fu<fx)
                                if(u>=x)
                                    a=x;
                                else
                                    b=x;
                                end
                                v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                            else
                                if (u<x)
                                    a=u;
                                else
                                    b=u;
                                end
                                if(fu<=fw)||(w==x)
                                    v=w; fv=fw; w=u; fw=fu;
                                elseif (fu<=fv)||(v==x)||(v==w)
                                    v=u; fv=fu;
                                end
                                continue
                            end
                        end
                    end
                    if (x>=xm) %CONDITION 3
                        e=a-x;
                    else
                        e=b-x;
                    end
                    d = CGOLD*e;
                    if(abs(d)>=tol1) %CONDITON 4
                        u=x+d;
                    else
                        u=x+abs(tol1)*sign(d);
                    end
                    r_f=u;I_S =Current_Source_loop;R= ramp_rate(1,ramp_rate_loop);
                    [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                        max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                        Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_SET(r_f,I_S,R);
                    Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                    fu=Free_Energy_eV; radius_filament_m =r_f;
                    if (fu<fx)
                        if(u>=x)
                            a=x;
                        else
                            b=x;
                        end
                        v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                    else
                        if (u<x)
                            a=u;
                        else
                            b=u;
                        end
                        if(fu<=fw)||(w==x)
                            v=w; fv=fw; w=u; fw=fu;
                        elseif (fu<=fv)||(v==x)||(v==w)
                            v=u; fv=fu;
                        end
                    end
                end
            end
            I_V_tbl = tbl;
            %ends the current loop if the total voltage exceed the input
            %voltage
            if source_VoltageSET_V<=V_pos_amp
                I_V_SET1 = vertcat(I_V_SET1,I_V_tbl);
            else
                break
            end
        end
        fil_temp = I_V_SET1.max_filament_Temperature_K;
        SET_Temp = fil_temp(end);
        device_Voltage_V = I_V_SET1.device_Voltage_V1;
        device_current_uA = I_V_SET1.device_current_uA1;
        I_V_SET=table(device_Voltage_V,device_current_uA);
        DATA_min_SET=vertcat(DATA_min_SET,I_V_SET1);
        I_V = vertcat(I_V,I_V_SET);
        s_r=I_V_SET1.radius_filament_m(end);
        
        %%%%%%%%%%%%%%%%%%%%%%%%Postitive Voltage%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%ON MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        model_ON = RRAM_ON;
        model_ON.hist.disable;
        format shortG
        model_ON.param.set('r_f',s_r);
        model_ON.param.set('alpha_fil',alpha_fil);
        model_ON.param.set('EC_HfO2x', EC_HfO2x);
        model_ON.param.set('R',ramp_rate(1,ramp_rate_loop));
        model_ON.param.set('V_amp', V_pos_amp);
        model_ON.param.set('t_rise',V_pos_amp/ramp_rate(1,ramp_rate_loop));
        model_ON.geom('geom1').run;
        model_ON.study('std1').feature('time').set('tunit', 's');
        model_ON.study('std1').feature('time').set('tlist', 'range(0,t_rise/30, 3/2*t_rise)');
        model_ON.study('std1').run;
        time = 0:(V_pos_amp/ramp_rate(1,ramp_rate_loop))/30:3/2*V_pos_amp/ramp_rate(1,ramp_rate_loop);
        time_row = time(:);time_s=flip(time_row);
        V_D_ON = mphglobal(model_ON, 'cir.IvsU1_v'); device_VoltageON_V = flip(V_D_ON);
        I_D_ON= 1e6*abs(mphglobal(model_ON, 'cir.IvsU1_i'));device_CurrentON_uA = flip(I_D_ON);
        V_S_ON = mphglobal(model_ON, 'cir.V1_v'); source_VoltageON_V = flip(V_S_ON);
        V_L_ON= mphglobal(model_ON, 'cir.R1_v');load_VoltageON_V = flip(V_L_ON);
        tbl_ON=table(time_s,device_VoltageON_V,device_CurrentON_uA,source_VoltageON_V,load_VoltageON_V);
        device_Voltage_V=device_VoltageON_V; device_current_uA=device_CurrentON_uA;
        I_V_ON = table(device_Voltage_V,device_current_uA);
        DATA_ON = vertcat(DATA_ON,tbl_ON);
        I_V = vertcat(I_V,I_V_ON);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%Negative voltage%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ON MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        model_ON_neg = RRAM_ON;
        model_ON_neg.hist.disable;
        format shortG
        model_ON_neg.param.set('r_f',s_r);
        model_ON_neg.param.set('alpha_fil',alpha_fil);
        model_ON_neg.param.set('EC_HfO2x', EC_HfO2x);
        model_ON_neg.param.set('R',ramp_rate(1,ramp_rate_loop));
        model_ON_neg.param.set('V_amp', V_neg_amp)
        model_ON_neg.param.set('t_rise',abs(V_neg_amp/ramp_rate(1,ramp_rate_loop)));
        model_ON_neg.geom('geom1').run;
        model_ON_neg.study('std1').feature('time').set('tunit', 's');
        model_ON_neg.study('std1').feature('time').set('tlist', 'range(0,t_rise/30, 3/2*t_rise)');
        model_ON_neg.study('std1').run;
        device_VoltageON_neg_V = mphglobal(model_ON_neg, 'cir.IvsU1_v');
        device_CurrentON_neg_uA = -1e6*abs(mphglobal(model_ON_neg, 'cir.IvsU1_i'));
        source_VoltageON_neg_V = mphglobal(model_ON_neg, 'cir.V1_v');
        load_VoltageON_neg_V = mphglobal(model_ON_neg, 'cir.R1_v');
        max_fil_TempON_neg_K = mphmax(model_ON_neg, 'T', 2, 'selection',3);
        time_neg = 0 : abs(V_neg_amp/ramp_rate(1,ramp_rate_loop))/30 :3/2*abs(V_neg_amp/ramp_rate(1,ramp_rate_loop));
        time_neg_s = time_neg(:);
        tbl_ON_neg=table(time_neg_s,device_VoltageON_neg_V,device_CurrentON_neg_uA,source_VoltageON_neg_V,load_VoltageON_neg_V);
        ON_negcount_coloumn = zeros(length(time_neg),1);
        count_i = 0;
        count_j2 = 0;
        %include voltages that corresponds to filament temperature smaller
        %than temperature corresponding to I_SET
        for loop = 1: length(time_neg)
            count_i = count_i+1;
            if max_fil_TempON_neg_K(loop)< SET_Temp
                count_j2 = count_j2+1;
                ON_negcount_coloumn(count_j2) = count_i;
            else
                break
            end
        end
        if count_j2 > 0
            ON_negcount_nonzero = nonzeros(ON_negcount_coloumn);
            start_neg = ON_negcount_nonzero(length(ON_negcount_nonzero),1);
        else
            [~, start_neg] = min(device_VoltageON_neg_V);
        end
        RESET_Voltage_V = device_VoltageON_neg_V(start_neg,1);
        device_Voltage_V = device_VoltageON_neg_V(1:start_neg,1);
        device_current_uA = device_CurrentON_neg_uA(1:start_neg,1);
        I_V_ON_neg = table(device_Voltage_V,device_current_uA);
        DATA_ON_neg = vertcat(DATA_ON_neg,tbl_ON_neg);
        I_V = vertcat(I_V,I_V_ON_neg);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESET MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        model_RESET = RRAM_RESET;
        model_RESET.hist.disable;
        format shortG
        %Assign the range of source voltage with step increament
        min_Voltage = RESET_Voltage_V;
        max_Voltage = V_neg_amp;
        step_increament_Voltage = -step_volt;
        %create empty table for database
        I_V_RESET1 = table;
        %Sweeps source voltage from minimum to maximum voltages in step increament
        for Voltage_Source_loop = min_Voltage:step_increament_Voltage:max_Voltage
            
            min_gap_length = 1e-9; max_gap_length = 4e-9; step_increament_gap = 1e-9;
            row =  min_gap_length:step_increament_gap:max_gap_length;
            FE_diff_RESET = zeros(length(row),1);
            DATA2 = table;
            stable_length = [];
            for length_gap_loop = min_gap_length:step_increament_gap:max_gap_length
                h_g=length_gap_loop; V_S= Voltage_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                    max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                    Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] =FE_RESET(h_g,s_r,V_S,R);
                Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                length_gap_m = h_g; Source_Voltage_V=Voltage_Source_loop;
                tbl=table(ramp_rate(1,ramp_rate_loop),Source_Voltage_V,device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,length_gap_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                    max_gap_Temperature_K,avg_gap_Temperature_K,Free_Energy_eV,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                    Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                DATA_RESET = vertcat(DATA_RESET,tbl);
                DATA2 = vertcat(DATA2,tbl);
                [num2str(cycle_loop),'|',num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(Source_Voltage_V),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA2),'|',num2str(device_Voltage_V2),'|',num2str(length_gap_m),'|',num2str(max_filament_Temperature_K),'|',num2str(avg_gap_Temperature_K),'|',num2str(delta_mu)]
                
            end
            for count_e = 1 : length(row)-1
                FE_diff_RESET(count_e,1) = (DATA2.Free_Energy_eV(count_e+1)-DATA2.Free_Energy_eV(count_e));
            end
            %locates the minimum in free energy using the difference and records the row corresponding to it
            for count_f = 1 : length(row)-2
                if (FE_diff_RESET(count_f)<0) && (FE_diff_RESET(count_f+1)>0)
                    stable_length = DATA2.length_gap_m(count_f+1);
                end
            end
            %if minimum exist then uses Brents minimization for speedy
            %convergence
            if (isempty(stable_length)==0)
                min_gap_length = stable_length - step_increament_gap;
                max_gap_length = stable_length + step_increament_gap;
                ITMAX=20; tol=1e-2; ZEPS=1e-11; CGOLD=0.3819660;
                ax=min_gap_length; cx=max_gap_length; bx=0.5*(min_gap_length+max_gap_length)+0.2e-9;
                a=ax; b=cx; v=bx;
                w=v; x=v; e=0;
                h_g=x; V_S=Voltage_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                    max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                    Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] =FE_RESET(h_g,s_r,V_S,R);
                Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                fx=Free_Energy_eV; fv=fx; fw=fx;
                length_gap_m = h_g ;
                Source_Voltage_V=Voltage_Source_loop;
                for iter = 1:ITMAX
                    [num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(Source_Voltage_V),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA2),'|',num2str(device_Voltage_V2),'|',num2str(length_gap_m),'|',num2str(avg_filament_Temperature_K),'|',num2str(avg_gap_Temperature_K),'|',num2str(delta_mu)]
                    tbl=table(ramp_rate(1,ramp_rate_loop),Source_Voltage_V,device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,length_gap_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                        max_gap_Temperature_K,avg_gap_Temperature_K,Free_Energy_eV,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                        Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                    DATA_RESET = vertcat(DATA_RESET,tbl);
                    xm = 0.5*(a+b); tol1=tol*abs(x)+ZEPS; tol2=2*tol1;
                    if (abs(x-xm)<=(tol2-0.5*(b-a))) %CONDITION 1
                        h_g=x; V_S=Voltage_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                        [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                            max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] =FE_RESET(h_g,s_r,V_S,R);
                        Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                        length_gap_m = h_g ;
                        tbl=table(ramp_rate(1,ramp_rate_loop),Source_Voltage_V,device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,length_gap_m,max_filament_Temperature_K,avg_filament_Temperature_K,...
                            max_gap_Temperature_K,avg_gap_Temperature_K,Free_Energy_eV,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu);
                        DATA_RESET = vertcat(DATA_RESET,tbl);
                        [num2str(ramp_rate(1,ramp_rate_loop)),'|',num2str(Source_Voltage_V),'|',num2str(device_Resistance_Ohm),'|', num2str(device_current_uA2),'|',num2str(device_Voltage_V2),'|',num2str(length_gap_m),'|',num2str(avg_filament_Temperature_K),'|',num2str(avg_gap_Temperature_K),'|',num2str(delta_mu)]
                        break
                    end
                    if (abs(e)>tol1) %CONDITION 2
                        r=(x-w)*(fx-fw); q=(x-v)*(fx-fw); p=(x-v)*q-(x-w)*r; q = 2*(q-r);
                        if (q > 0)
                            p=-p;
                        end
                        q = abs(q); etemp=e; e=d;
                        if (abs(p)>=abs(0.5*q*etemp))||(p<=q*(a-x))||(p>=q*(b-x)) %CONDITON 2.1
                            if (x>=xm)
                                e=a-x;
                            else
                                e=b-x;
                            end
                            d=CGOLD*e;
                            if(abs(d)>=tol1)
                                u=x+d;
                            else
                                u=x+abs(tol1)*sign(d);
                            end
                            h_g=u; V_S=Voltage_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                            [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                                max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                                Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_RESET(h_g,s_r,V_S,R);
                            Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                            fu=Free_Energy_eV; length_gap_m=h_g;
                            if (fu<fx)
                                if(u>=x)
                                    a=x;
                                else
                                    b=x;
                                end
                                v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                            else
                                if (u<x)
                                    a=u;
                                else
                                    b=u;
                                end
                                if(fu<=fw)||(w==x)
                                    v=w; fv=fw; w=u; fw=fu;
                                elseif (fu<=fv)||(v==x)||(v==w)
                                    v=u; fv=fu;
                                end
                                continue
                            end
                        else
                            d=p/q; u=x+d;
                            if (u-a<tol2)||(b-u<tol2)
                                d=abs(tol1)*sign(xm-x);
                            end
                            if(abs(d)>=tol1)
                                u=x+d;
                            else
                                u=x+abs(tol1)*sign(d);
                            end
                            h_g=u; V_S=Voltage_Source_loop; R=ramp_rate(1,ramp_rate_loop);
                            [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                                max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                                Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_RESET(h_g,s_r,V_S,R);
                            Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                            fu=Free_Energy_eV; length_gap_m =h_g;
                            if (fu<fx)
                                if(u>=x)
                                    a=x;
                                else
                                    b=x;
                                end
                                v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                            else
                                if (u<x)
                                    a=u;
                                else
                                    b=u;
                                end
                                if(fu<=fw)||(w==x)
                                    v=w; fv=fw; w=u; fw=fu;
                                elseif (fu<=fv)||(v==x)||(v==w)
                                    v=u; fv=fu;
                                end
                                continue
                            end
                        end
                    end
                    if (x>=xm) %CONDITION 3
                        e=a-x;
                    else
                        e=b-x;
                    end
                    d = CGOLD*e;
                    if(abs(d)>=tol1) %CONDITON 4
                        u=x+d;
                    else
                        u=x+abs(tol1)*sign(d);
                    end
                    h_g=u;V_S =Voltage_Source_loop;R= ramp_rate(1,ramp_rate_loop);
                    [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
                        max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
                        Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_RESET(h_g,s_r,V_S,R);
                    Free_Energy_eV=Thermal_energy_filament_eV+Thermal_energy_dielectric_eV+Electro_Static_Energy_eV+Surface_Energy_eV+Volume_Energy_eV;
                    fu=Free_Energy_eV; length_gap_m =h_g;
                    if (fu<fx)
                        if(u>=x)
                            a=x;
                        else
                            b=x;
                        end
                        v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
                    else
                        if (u<x)
                            a=u;
                        else
                            b=u;
                        end
                        if(fu<=fw)||(w==x)
                            v=w; fv=fw; w=u; fw=fu;
                        elseif (fu<=fv)||(v==x)||(v==w)
                            v=u; fv=fu;
                        end
                    end
                end
                I_V_tbl = tbl;
                I_V_RESET1 = vertcat(I_V_RESET1,I_V_tbl);
            end
        end
        
        if (isempty(I_V_RESET1)==0)
            device_Voltage_V = I_V_RESET1.device_Voltage_V2;
            device_current_uA = I_V_RESET1.device_current_uA2;
            I_V_RESET=table(device_Voltage_V,device_current_uA);
            DATA_min_RESET=vertcat(DATA_min_RESET,I_V_RESET1);
            I_V = vertcat(I_V,I_V_RESET);
            s_g = I_V_RESET1.length_gap_m(end);
        else
            s_g = 1e-9;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%Negative voltage%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%OFF_MODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        model_OFF_neg= RRAM_OFF;
        model_OFF_neg.hist.disable;
        model_OFF_neg.param.set('V_amp', V_neg_amp);
        model_OFF_neg.param.set('r_f',s_r);
        model_OFF_neg.param.set('h_g',s_g);
        model_OFF_neg.param.set('EC_HfO2x', EC_HfO2x);
        model_OFF_neg.param.set('EC_gap', EC_gap);
        model_OFF_neg.param.set('alpha_fil',alpha_fil);
        model_OFF_neg.param.set('alpha_gap',alpha_gap);
        model_OFF_neg.param.set('k_eff',k_eff);
        model_OFF_neg.param.set('R',ramp_rate(1,ramp_rate_loop));
        model_OFF_neg.param.set('t_rise',abs(V_neg_amp/ramp_rate(1,ramp_rate_loop)));
        model_OFF_neg.study('std1').feature('time').set('tunit', 's');
        model_OFF_neg.study('std1').feature('time').set('tlist', 'range(0,t_rise/30, 3/2*t_rise)');
        model_OFF_neg.study('std1').run;
        %       max_dielec_Temperature_neg= mphmax(model_OFF_neg, 'T', 2, 'selection',[4 9]);
        %       max_dielec_Temperature_neg_K = max_dielec_Temperature_neg(:);
        V_D_OFFneg = mphglobal(model_OFF_neg, 'cir.IvsU1_v'); device_VoltageOFF_neg_V = flip(V_D_OFFneg);
        I_D_OFFneg = -1e6*abs(mphglobal(model_OFF_neg, 'cir.IvsU1_i'));device_CurrentOFF_neg_uA = flip(I_D_OFFneg);
        V_S_OFFneg = mphglobal(model_OFF_neg, 'cir.V1_v');source_VoltageOFF_neg_V = flip(V_S_OFFneg);
        V_L_OFFneg = mphglobal(model_OFF_neg, 'cir.R1_v');load_VoltageOFF_neg_V =flip(V_L_OFFneg);
        time = 0:(V_pos_amp/ramp_rate(1,ramp_rate_loop))/30:3/2*V_pos_amp/ramp_rate(1,ramp_rate_loop);
        time_row = time(:);time_s=flip(time_row);
        tbl_OFF = table(time_s,device_VoltageOFF_neg_V,device_CurrentOFF_neg_uA,source_VoltageOFF_neg_V,load_VoltageOFF_neg_V);
        device_Voltage_V = device_VoltageOFF_neg_V;
        device_current_uA = device_CurrentOFF_neg_uA;
        I_V_OFF_neg = table(device_Voltage_V,device_current_uA);
        DATA_OFF_neg = vertcat(DATA_OFF_neg,tbl_OFF);
        I_V = vertcat(I_V,I_V_OFF_neg);
        
        
    end
end

save('Complete_I_V_ramprate_DATA.mat','I_V','DATA_OFF','DATA_min_SET','DATA_SET','DATA_ON','DATA_ON_neg','DATA_min_RESET','DATA_RESET','DATA_OFF_neg')

%%%%%%%%%%%%%%%%SET Control%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [device_current_uA1,source_VoltageSET_V,Load_VoltageSET_V,device_Voltage_V1,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
            max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_SET(r_f,I_S,R)
        %sets the Source voltage and filament radius in COMSOL file and solves the problem
        model_SET.param.set('I_S',I_S);
        model_SET.param.set('r_f',r_f);
        A=exp(-alpha_fil*log((V_pos_amp)/(R*atomic_vibration)));
        model_SET.param.set('A',A);
        model_SET.param.set('EC_HfO2x', EC_HfO2x);
        model_SET.geom('geom1').run;
        model_SET.study('std1').run;
        %asks COMSOL to calculate mean filament temperature, temperature distribution in the device, potential distribution, current density, and electric field in cylindrical coordinate
        avg_filament_Temperature_K = mphmean(model_SET, 'T' , 2,'selection',3);
        max_filament_Temperature_K = mphmax(model_SET, 'T', 2, 'selection',3);
        avg_dielectric_Temperature_K = mphmean(model_SET, 'T' , 2,'selection',7);
        max_dielectric_Temperature_K = mphmax(model_SET, 'T', 2, 'selection',7);
        device_Voltage_V1 = mphglobal(model_SET, 'cir.IvsU1_v');
        device_current_uA1 = 1e6*(mphglobal(model_SET, 'cir.IvsU1_i'));
        Load_VoltageSET_V = mphglobal(model_SET, 'cir.R1_v');
        device_Resistance_Ohm = 1e6*device_Voltage_V1/device_current_uA1;
        source_VoltageSET_V = device_Voltage_V1+Load_VoltageSET_V;
        delta_mu = delta_mu_SET_J+beta_1*(1/delta_W_uc-1/delta_W_i)*(1/1.6e-19)*boltzmann_constant*avg_dielectric_Temperature_K*log(V_pos_amp/(R*atomic_vibration));
        if delta_mu < 1e7
            delta_mu = 1e7;
        end
        Electro_Static_Energy_eV = 6.2415e+18*mphint2(model_SET,'((ec.Er)^2+(ec.Ez)^2+(ec.Ephi)^2)*8.85*10^(-12)*25*2*pi*r*(1/2)',2,'selection',7);
        Thermal_energy_electrode_eV = 6.2415e+18*mphint2(model_SET,'2*pi*5.22e3*545.33*r*(T-293.15)',2,'selection',[2,5]);
        Thermal_energy_OEL_eV= 6.2415e+18*mphint2(model_SET,'2*pi*13.31e3*144*r*(T-293.15)',2,'selection',4);
        Thermal_energy_filament_eV = 6.2415e+18*mphint2(model_SET,'2*pi*12e3*140*r*(T-293.15)',2,'selection',3);
        Thermal_energy_dielectric_eV = 6.2415e+18*mphint2(model_SET,'2*pi*10e3*120*r*(T-293.15)',2,'selection',7);
        Surface_Energy_eV = 6.2415e+18* 2*pi*r_f*height_filament*surface_tension;
        Volume_Energy_eV= 6.2415e+18*pi*r_f^2*height_filament*delta_mu;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESET control%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [device_current_uA2,Load_VoltageSET_V,device_Voltage_V2,device_Resistance_Ohm,max_filament_Temperature_K,avg_filament_Temperature_K,...
            max_gap_Temperature_K,avg_gap_Temperature_K,max_dielectric_Temperature_K,avg_dielectric_Temperature_K,Thermal_energy_electrode_eV,Thermal_energy_filament_eV,Thermal_energy_dielectric_eV,...
            Thermal_energy_OEL_eV,Electro_Static_Energy_eV,Surface_Energy_eV,Volume_Energy_eV,delta_mu] = FE_RESET(h_g,s_r,V_S,R)
        %sets the Source voltage and filament radius in COMSOL file and solves the problem
        model_RESET.param.set('V_S',V_S);
        model_RESET.param.set('h_g',h_g);
        model_RESET.param.set('r_f',s_r);
        A=exp(-alpha_fil*log(-V_neg_amp/(R*atomic_vibration)));
        B=exp(-alpha_gap*log(-V_neg_amp/(R*atomic_vibration)));
        model_RESET.param.set('A',A);
        model_RESET.param.set('B',B);
        model_RESET.param.set('k_eff',k_eff);
        model_RESET.param.set('EC_HfO2x', EC_HfO2x);
        model_RESET.param.set('EC_gap', EC_gap);
        model_RESET.geom('geom1').run;
        model_RESET.study('std1').run;
        %asks COMSOL to calculate mean filament temperature, temperature distribution in the device, potential distribution, current density, and electric field in cylindrical coordinate
        %multiple values of source voltage (-0.1,V_S/2,V_S) is required to solve thus taking the
        %last value which corresponds to the voltage of interest
        avgfilT = mphmean(model_RESET, 'T' , 2,'selection',4); avg_filament_Temperature_K = avgfilT(end);
        maxfilT= mphmax(model_RESET, 'T', 2, 'selection',4); max_filament_Temperature_K =maxfilT(end);
        avggapT= mphmean(model_RESET, 'T' , 2,'selection',3); avg_gap_Temperature_K = avggapT(end);
        maxgapT= mphmax(model_RESET, 'T', 2, 'selection',3);max_gap_Temperature_K = maxgapT(end);
        avgdielT= mphmean(model_RESET, 'T' , 2,'selection',8);avg_dielectric_Temperature_K=avgdielT(end);
        max_dielT= mphmax(model_RESET, 'T', 2, 'selection',8);max_dielectric_Temperature_K=max_dielT(end);
        
        devV= mphglobal(model_RESET, 'cir.IvsU1_v');device_Voltage_V2 = devV(end);
        devI= 1e6*(mphglobal(model_RESET, 'cir.IvsU1_i'));device_current_uA2 =devI(end);
        LoadV= mphglobal(model_RESET, 'cir.R1_v');Load_VoltageSET_V =LoadV(end);
        device_Resistance_Ohm = 1e6*device_Voltage_V2/device_current_uA2;
        
        delta_mu = delta_mu_RESET_J + beta_2*(1/delta_W_uc-1/delta_W_mc)*(1/1.6e-19)*boltzmann_constant*avg_filament_Temperature_K*log(abs(V_neg_amp)/(R*atomic_vibration));
        if delta_mu < 1e7
            delta_mu = 1e7;
        end
        El_Stat= 6.2415e+18*mphint2(model_RESET,'((ec.Er)^2+(ec.Ez)^2+(ec.Ephi)^2)*8.85*10^(-12)*25*2*pi*r*(1/2)',2,'selection',8);
        Electro_Static_Energy_eV = El_Stat(end);
        Th_El= 6.2415e+18*mphint2(model_RESET,'2*pi*5.22e3*545.33*r*(T-293.15)',2,'selection',[2,6]);
        Thermal_energy_electrode_eV = Th_El(end);
        Th_OEL= 6.2415e+18*mphint2(model_RESET,'2*pi*13.31e3*144*r*(T-293.15)',2,'selection',5);
        Thermal_energy_OEL_eV = Th_OEL(end);
        Th_fil= 6.2415e+18*mphint2(model_RESET,'2*pi*12e3*140*r*(T-293.15)',2,'selection',4);
        Thermal_energy_filament_eV = Th_fil(end);
        Th_di= 6.2415e+18*mphint2(model_RESET,'2*pi*10e3*120*r*(T-293.15)',2,'selection',[3,8]);
        Thermal_energy_dielectric_eV=Th_di(end);
        Surface_Energy_eV = 6.2415e+18* 2*pi*s_r*h_g*surface_tension;
        Volume_Energy_eV= 6.2415e+18*pi*s_r^2*h_g*delta_mu;
    end
end
