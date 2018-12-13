close all
clear all


%% Material 1, 1060 Alloy

syms yield_all
syms x_all

M_all = csvread('alloy.csv');
X_all = M_all(1:10);
Y_all = transpose(M_all(1:10,2));

modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));
beta0_all = [3;4;4;13];
beta_all = fitnlm(X_all,Y_all,modelfun,beta0_all);

result_all = beta_all.Coefficients.Estimate;

figure(1)
set(figure,'color','w');
subplot(2,2,1);
scatter(X_all,Y_all)
title('1060 Alloy'),xlabel("Thickness (mm)"),ylabel("Max Von Mises Stress (MPa)"); 

y_all = 27.5;
reversemodelfun_all = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
radius_all = reversemodelfun_all(result_all, y_all)

xinterval = [0.25 1.375];
hold on
plot(X_all, modelfun(result_all, X_all))
fplot(subs(yield_all,27.5), xinterval)
legend("Data points","Regression model","Material yield strength = 27.5 MPa")

%% 

% Material 2 - ABS

syms yield_abs

M_abs = csvread('alloy.csv');


X_abs = M_abs(11:20);
Y_abs = transpose(M_abs(11:20,2));

modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

beta0_abs = [3;4;4;13];

beta_abs = fitnlm(X_abs,Y_abs,modelfun,beta0_abs);

result_abs = beta_abs.Coefficients.Estimate;

subplot(2,2,2);
scatter(X_abs,Y_abs)
title('ABS'),xlabel("Thickness (mm)"),ylabel("Max Von Mises Stress (MPa)");

y_abs = 35;
reversemodelfun_abs = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
radius_abs = reversemodelfun_abs(result_abs, y_abs)

xinterval = [0.25 1.375];
hold on
plot(X_abs, modelfun(result_abs, X_abs))
fplot(subs(yield_abs,35), xinterval)
legend("Data points","Regression model","Material yield strength = 35 MPa")

%% Copper

syms yield_cop

M_cop = csvread('alloy.csv');

% Material 1060 Alloy
X_cop = M_cop(21:30);
Y_cop = transpose(M_cop(21:30,2));

modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

beta0_cop = [3;4;4;13];

beta_cop = fitnlm(X_cop,Y_cop,modelfun,beta0_cop);

result_cop = beta_cop.Coefficients.Estimate;

y_cop = 34;
reversemodelfun_cop = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
radius_cop = reversemodelfun_cop(result_cop, y_cop)

subplot(2,2,3);
scatter(X_cop,Y_cop)
title('Copper Alloy'),xlabel("Thickness (mm)"),ylabel("Max Von Mises Stress (MPa)");

xinterval = [0.25 1.375];
hold on
plot(X_cop, modelfun(result_cop, X_cop))
fplot(subs(yield_cop,34), xinterval)
legend("Data points","Regression model","Material yield strength = 34 MPa")

%% Zinc

syms yield_zin

M_zin = csvread('alloy.csv');

% Material 1060 Alloy
X_zin = M_zin(31:40);
Y_zin = transpose(M_zin(31:40,2));

modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

beta0_zin = [3;4;4;13];

beta_zin = fitnlm(X_zin,Y_zin,modelfun,beta0_zin);

result_zin = beta_zin.Coefficients.Estimate;

y_zin = 35.9;
reversemodelfun_zin = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
radius_zin = reversemodelfun_zin(result_zin, y_zin)

subplot(2,2,4);
scatter(X_zin,Y_zin)
title('Zinc-Aluminium Alloy'),xlabel("Thickness (mm)"),ylabel("Max Von Mises Stress (MPa)");

xinterval = [0.25 1.375];
hold on
plot(X_zin, modelfun(result_zin, X_zin))
fplot(subs(yield_zin,35.9), xinterval)
legend("Data points","Regression model","Material yield strength = 35.9")


%% Multi-Objective


% Aluminium Alloy
% outer loop for radius
% inner loop for material
data = table()



materials = ["Aluminium","ABS","Copper","Zinc-Aluminium Alloy"]
radii = [0:0.1:1.5]

price = [1.62 2.25 4.48 2.42]
density = [2900 1200 8943 7000]
yield_max = [27.5 35 34  35.9]

betas = [result_all result_abs result_cop result_zin]

for x = 1:length(materials)
    for i = 1:length(radii)
        current_material = materials(x)
        radius_current = radii(i);
        current_beta = betas(:,x);
        current_stress = modelfun(current_beta,radius_current)
        
        area = (pi()*(0.28+(radius_current/1000))^2) - (pi()*0.28^2)
        cylinder = area*0.4
        sides = (pi()*(0.28+(radius_current/1000))^2)*(radius_current/1000)
        volume = cylinder + sides
        
        price_current = price(x)
        density_current = density(x)
        mass = volume * density_current
        cost = price_current * mass
        
        r_p = 0.0716
        omega_p = 250
        r_d = 0.28+(i/1000)
        omega_d = (r_p*omega_p)/r_d
   
        % Add data to 
        new_row = {current_material, radius_current, cost, omega_d}
        
        
        if current_stress < yield_max(x)
            data = [data;new_row]
            
        end
        
    end
end





