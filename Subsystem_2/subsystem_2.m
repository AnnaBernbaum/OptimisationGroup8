close all
clear all

tic

%% Material 1 - 1060 Alloy

syms yield_all

% Loading the datapoints file and reading the values for the first material
M = csvread('datapoints.csv');
X_all = M(1:10);
Y_all = transpose(M(1:10,2));

% Initial model for fitting the datapoints
modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

% Initial Guesses
beta0_all = [3;4;4;13];

% Non-linear regression
beta_all = fitnlm(X_all,Y_all,modelfun,beta0_all);

% Extracts the coefficients
result_all = beta_all.Coefficients.Estimate;

% Plot the datapoints
set(figure,'color','w');
subplot(2,2,1);
scatter(X_all,Y_all)
title('1060 Aluminium Alloy','fontsize',18);
xlabel("Thickness (mm)",'fontsize',16);
ylabel("Max Von Mises Stress (MPa)",'fontsize',16); 

% Find the minimum thickness for where the material fails
y_all = 27.5;
reversemodelfun_all = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
thickness_all = reversemodelfun_all(result_all, y_all)

% Plotting the regresstion and line for the material yield strength
xinterval = [0.25 1.375];
hold on
plot(X_all, modelfun(result_all, X_all))
fplot(subs(yield_all,27.5), xinterval)
legend({"Data points","Regression model","Material yield strength = 27.5 MPa"},'fontsize',16)

%% Material 2 - ABS

syms yield_abs

% Reading the values for the second material
X_abs = M(11:20);
Y_abs = transpose(M(11:20,2));

% Initial model for fitting the datapoints
modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

% Initial Guesses
beta0_abs = [3;4;4;13];

% Non-linear regression
beta_abs = fitnlm(X_abs,Y_abs,modelfun,beta0_abs);

% Extracts the coefficients
result_abs = beta_abs.Coefficients.Estimate;

% Plot the datapoints
subplot(2,2,2);
scatter(X_abs,Y_abs)
title('ABS','fontsize',18);
xlabel("Thickness (mm)",'fontsize',16);
ylabel("Max Von Mises Stress (MPa)",'fontsize',16);

% Find the minimum thickness for where the material fails
y_abs = 35;
reversemodelfun_abs = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
thickness_abs = reversemodelfun_abs(result_abs, y_abs)

% Plotting the regresstion and line for the material yield strength
xinterval = [0.25 1.375];
hold on
plot(X_abs, modelfun(result_abs, X_abs))
fplot(subs(yield_abs,35), xinterval)
legend({"Data points","Regression model","Material yield strength = 35 MPa"},'fontsize',18)

%% Material 3 - Copper Alloy

syms yield_cop

% Reading the values for the third material
X_cop = M(21:30);
Y_cop = transpose(M(21:30,2));

% Initial model for fitting the datapoints
modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

% Initial guesses
beta0_cop = [3;4;4;13];

% Non-linear regression
beta_cop = fitnlm(X_cop,Y_cop,modelfun,beta0_cop);

% Extracts the coefficients
result_cop = beta_cop.Coefficients.Estimate;

% Find the minimum thickness for where the material fails
y_cop = 34;
reversemodelfun_cop = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
thickness_cop = reversemodelfun_cop(result_cop, y_cop)

% Plot the datapoints
subplot(2,2,3);
scatter(X_cop,Y_cop)
title('Copper Alloy','fontsize',18);
xlabel("Thickness (mm)",'fontsize',16);
ylabel("Max Von Mises Stress (MPa)",'fontsize',16);

% Plotting the regresstion and line for the material yield strength
xinterval = [0.25 1.375];
hold on
plot(X_cop, modelfun(result_cop, X_cop))
fplot(subs(yield_cop,34), xinterval)
legend({"Data points","Regression model","Material yield strength = 34 MPa"},'fontsize',18)

%% Zinc Alloy

syms yield_zin

% Reading the values for the fourth material
X_zin = M(31:40);
Y_zin = transpose(M(31:40,2));

% Initial model for fitting the datapoints
modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

% Initial guesses
beta0_zin = [3;4;4;13];

% Non-linear regression
beta_zin = fitnlm(X_zin,Y_zin,modelfun,beta0_zin);

% Extracts the coefficients
result_zin = beta_zin.Coefficients.Estimate;

% Find the minimum thickness for where the material fails
y_zin = 35.9;
reversemodelfun_zin = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
thickness_zin = reversemodelfun_zin(result_zin, y_zin)

% Plot the datapoints
subplot(2,2,4);
scatter(X_zin,Y_zin)
title('Zinc-Aluminium Alloy','fontsize',18);
xlabel("Thickness (mm)",'fontsize',16);
ylabel("Max Von Mises Stress (MPa)",'fontsize',16);

% Plotting the regresstion and line for the material yield strength
xinterval = [0.25 1.375];
hold on
plot(X_zin, modelfun(result_zin, X_zin))
fplot(subs(yield_zin,35.9), xinterval)
legend({"Data points","Regression model","Material yield strength = 35.9"},'fontsize',18)

%% Multi-Objective Optimisation

% Initiate table for results
data = table();

% Arrays for looping through materials/radii and their values
materials = ["Aluminium","ABS","Copper","Zinc-Aluminium Alloy"];
thickness = [0:0.1:1.5];
price = [1.62 2.25 4.48 2.42];
density = [2900 1213 8943 7000];
yield_max = [27.5 35 34  35.9];

% Fetch the beta values from each regression model
betas = [result_all result_abs result_cop result_zin];

% Looping through each combination of materials and radii 
for x = 1:length(materials)
    for i = 1:length(thickness)
        
        % Finding the current values for each iteration
        current_material = materials(x);
        thickness_current = thickness(i);
        current_beta = betas(:,x);
        current_stress = modelfun(current_beta,thickness_current);
        
        % Calculating the volume for current iteration
        area = (pi()*(0.28+(thickness_current/1000))^2) - (pi()*0.28^2);
        cylinder = area*0.4;
        sides = (pi()*(0.28+(thickness_current/1000))^2)*(thickness_current/1000);
        volume = cylinder + sides;
        
        % Calculating the cost for current iteration
        price_current = price(x);
        density_current = density(x);
        mass = volume * density_current;
        cost = price_current * mass;
        
        % Calculating omega_d
        r_p = 0.0716;
        omega_p = 250;
        r_d = 0.28+(i/1000);
        omega_d = (r_p*omega_p)/r_d;
   
        % Add the results from current iteration to data array 
        new_row = {current_material, thickness_current, cost, omega_d};
        
        % But only the solutions that are successful 
        if current_stress < yield_max(x)
            data = [data;new_row];
            
        end
        
    end
end

toc 
t = toc
