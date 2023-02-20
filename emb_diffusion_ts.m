%% 1D diffusion of H2O through a rhyolitic melt to one bubble during two-stage decompression (optimized)

% 1D diffusion of H2O and CO2 through a rhyolitic melt to one bubble during two-stage decompression 
% Optimized for first and second stages of decompression (dPdt_stage1 and dPdt_stage1) and transition pressure (P_int)

% Script created by Hosseini et al. (in prep)

clear all % Clear any existing variables in the workspace
close all % Close all figures
clc % Clear command window

%% ASSUMPTIONS
% All diffusion occuring in one dimension from interior to outlet of embayment
% No volume change of bubble due to mass influx (i.e., melt in embayment is incompressible)
% Melt and bubble are initially in equilibrium 
% Isothermal
% Cartesian coordinates

%% ENTER VARIABLE INPUT PARAMETERS
% The following are all of the input parameters that can be adjusted for each sample.
              
P_i = 150; % Initial pressure [MPa]
P_f = 30; % Final (quench/fragmentation) pressure [MPa]
lengthi = 150; % Distance from interior of embayment to glass-bubble interface [um]
radius = 0.1; % Bubble radius (1/2 mm) [um]
H2O_i = 3.3; % Initial H2O concentration from melt inclusion or embayment interior plateau [wt. %]
rho = 2350; % Density of rhyolite [kg/m^3]-estimated density of obsidian from White and Harris (1997) 
TC = 780; % Temperature [deg C]
errorH2O = 0.2; % Analytical uncertainty in H2O measurement [wt %]
errorCO2 = 20; % Analytical uncertainty in CO2 measurement [ppm]

%% ENTER MEASURED H2O AND CO2 CONCENTRATIONS
% The following arrays store the measured H2O and CO2 concentrations, as well as the distance along the 
% embayment at which measurements were made. Note: distances and measurements go from outlet/bubble to interior

% SAMPLE NAME: E3-2
MeasH2OCon = [2.32 2.65 2.86 2.99 3.04 3.10 3.14 3.16 2.88]; % Measured H2O concentration profile
MeasCO2Con = [180.80 233.70 239.82 251.69 257.63 263.56 263.63 263.66 280.80]; % Measured CO2 concentration profile
MeasDist = [20 30 40 50 60 70 80 90 100]; % Distance along profile, from outlet to interior

%% ENTER RANGE OF DPDT FOR THE FIRST AND SECOND STAGES OF DECOMPRESSION AND THE RANGE OF INTERMEDIATE PRESSURES TO ITERATE THROUGH. MIDDLE NUMBER IS STEP SIZE.

dPdt_stage1 = (0.005:0.001:0.5); % dPdt (stage 1) to iterate through to find best fit to data [MPa/s]
dPdt_stage2 = (0.005:0.001:0.5); % dPdt (stage 2) to iterate through to find best fit to data [MPa/s]
P_int = (50:10:150); % P_int to iterate through to find best fit to data [MPa]

%% CALCULATE DPDT_STAGE1, DPDT_STAGE2, AND P_INT THAT GENERATE THE LEAST MISFIT WITH THE MEASURED DATA.

misfit = zeros(length(dPdt_stage1),length(dPdt_stage2),length(P_int));

% Iterate through array of dPdt_stage1, dPdt_stage2, and P_int
for idpdt_1 = 1:length(dPdt_stage1)
    disp(['Running First Stage Decompression ',num2str(idpdt_1),' of ',num2str(length(dPdt_stage1))])
    
    for idpdt_2 = 1:length(dPdt_stage2)
    disp(['Running Second Stage Decompression ',num2str(idpdt_2),' of ',num2str(length(dPdt_stage2))])
    
        for int_p = 1:length(P_int)
        disp(['Running Intermediate Pressure ',num2str(int_p),' of ',num2str(length(P_int))])

        % Call the 1D diffusion function (diffusion_function_ts.m)
        [nodes, H2O_array, CO2_array] = diffusion_function_ts(P_i, P_int(int_p), P_f, lengthi,...
            radius, H2O_i, rho, TC, dPdt_stage1(idpdt_1), dPdt_stage2(idpdt_2));

        % Analyze the output
        misfitModel = 0;
        for im = 1:length(MeasDist)

            % Iterate through measured H2O and CO2 concentrations and distances
            d = MeasDist(im); 
            H2O = MeasH2OCon(im);
            CO2 = MeasCO2Con(im);

            % Find the node that is closest to our measured distance
            % Tilde indicates that minimum value is discarded, and only the index of the minimum value is used
            [~,iN] = min(abs(nodes - d));
            H2O_model = H2O_array(iN);
            CO2_model = CO2_array(iN);

            misfitModel = misfitModel + ((((H2O_model - H2O) / errorH2O).^2) + ((CO2_model - CO2) / errorCO2).^2) / length(MeasCO2Con);
        end
        misfit(idpdt_1,idpdt_2,int_p) = misfitModel;
        end
    end
end

[~,loc] = min(misfit(:));
[IDPDT_1,IDPDT_2,INT_P] = ind2sub(size(misfit),loc);

[nodes, H2O_array, CO2_array] = diffusion_function_ts(P_i, P_int(INT_P), P_f, lengthi, radius, H2O_i, rho, TC, dPdt_stage1(IDPDT_1), dPdt_stage2(IDPDT_2));

disp(['Minimum misfit = ', num2str(min(misfit(:)))])
disp(['First stage dP/dt = ',num2str(dPdt_stage1(IDPDT_1))])
disp(['Second stage dP/dt  = ',num2str(dPdt_stage2(IDPDT_2))])
disp(['Intermediate pressure = ',num2str(P_int(INT_P))])

if min(misfit(:)) < 1
    disp(['Good fit achieved with parameters.'])
else
    disp(['Good fit not achieved with parameters.']) 
end

%% PLOT DIFFUSION PROFILES

figure(1)
clf
hold on

h = subplot(2,1,1);
meH2O = plot(MeasDist,MeasH2OCon,'s'); % Plot the measured H2O profile
error = ones(1,length(MeasH2OCon))*errorH2O; % Initiate an array containing the H2O analytical uncertainty
e = errorbar(MeasDist,MeasH2OCon,error,'s','MarkerEdgeColor','k','MarkerFaceColor','#EDB120');
e.Color = 'black';
e.MarkerSize = 10;
hold on
plot(nodes,H2O_array,'k-','LineWidth',3) % Plot the modeled H2O profile
hold on
xlim([0 lengthi])
ylim([1 max(MeasH2OCon)+1])
uistack(e,'top');
xlabel('Distance (\mum)','fontweight','bold','FontSize',12)
ylabel('H_{2}O (wt. %)','fontweight','bold','FontSize',12)

i = subplot(2,1,2);
meCO2 = plot(MeasDist,MeasCO2Con,'s'); % Plot the measured CO2 profile
error = ones(1,length(MeasCO2Con))*errorCO2; % Initiate an array containing the CO2 analytical uncertainty
e = errorbar(MeasDist,MeasCO2Con,error,'s','MarkerEdgeColor','k','MarkerFaceColor','#EDB120');
e.Color = 'black';
e.MarkerSize = 10;
hold on
plot(nodes,CO2_array,'k-','LineWidth',3) % Plot the modeled CO2 profile
hold on
xlim([0 lengthi])
ylim([0 max(MeasCO2Con)+200])
uistack(e,'top');
xlabel('Distance (\mum)','fontweight','bold','FontSize',12)
ylabel('CO_{2} (ppm)','fontweight','bold','FontSize',12)