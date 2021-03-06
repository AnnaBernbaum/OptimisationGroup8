function betas = drum_stress_regression(filename)
    syms yield_all
   
    % Material 1 - 1060 Alloy
    M  = csvread(filename);
    X_all = M(1:10);
    Y_all = transpose(M(1:10,2));

    % Initial model for fitting the datapoints
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));
    beta0_all = [3;4;4;13];
    
    % Non-linear regression
    beta_all = fitnlm(X_all,Y_all,modelfun,beta0_all);
    
    % Extracts the coefficients
    result_all = beta_all.Coefficients.Estimate;

    % Find the minimum thickness for where the material fails
    y_all = 27.5;
    reversemodelfun_all = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
    thickness_all = reversemodelfun_all(result_all, y_all);

    % Material 2 - ABS
    syms yield_abs

    X_abs = M(11:20);
    Y_abs = transpose(M(11:20,2));

    % Initial model for fitting the datapoints
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));
    beta0_abs = [3;4;4;13];
    
    % Non-linear regression
    beta_abs = fitnlm(X_abs,Y_abs,modelfun,beta0_abs);
    
    % Extracts the coefficients
    result_abs = beta_abs.Coefficients.Estimate;

    % Find the minimum thickness for where the material fails
    y_abs = 35;
    reversemodelfun_abs = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
    thickness_abs = reversemodelfun_abs(result_abs, y_abs);


    % Material 3 - Copper Alloy
    syms yield_cop

    X_cop = M(21:30);
    Y_cop = transpose(M(21:30,2));

    % Initial model for fitting the datapoints
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));
    beta0_cop = [3;4;4;13];
    
    % Non-linear regression
    beta_cop = fitnlm(X_cop,Y_cop,modelfun,beta0_cop);

    % Extracts the coefficients
    result_cop = beta_cop.Coefficients.Estimate;

    % Find the minimum thickness for where the material fails
    y_cop = 34;
    reversemodelfun_cop = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
    thickness_cop = reversemodelfun_cop(result_cop, y_cop);

    % Material 4 - Zinc Alloy
    syms yield_zin
    
    X_zin = M(31:40);
    Y_zin = transpose(M(31:40,2));

    % Initial model for fitting the datapoints
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));
    beta0_zin = [3;4;4;13];

    % Non-linear regression
    beta_zin = fitnlm(X_zin,Y_zin,modelfun,beta0_zin);
    
    % Extracts the coefficients
    result_zin = beta_zin.Coefficients.Estimate;

    % Find the minimum thickness for where the material fails
    y_zin = 35.9;
    reversemodelfun_zin = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));
    thickness_zin = reversemodelfun_zin(result_zin, y_zin);
    
    betas = [result_all result_abs result_cop result_zin];
end
