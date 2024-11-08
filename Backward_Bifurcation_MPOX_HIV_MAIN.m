clear;  % Clear all variables from the workspace to start fresh.

% Define the parameters used in the model
Lambda = 1520000.*0.003; mu = 0.003; alpha = 1/8.5; gamma = 1/14; nuepsilon = 0;
p_m = 0.19; p_h = 0.012; eta = 0.15; rho = 1.2;
c = 0.85; Omega = 0.75; beta_11 = (p_m.*c.*mu)/Lambda;

sigma = 1.5; 
sigmaT = 1.5;


rend = 100000;  % Define the number of points for calculations
r1st = linspace(0,1.1,rend);  % Create a linearly spaced vector from 0 to 1.1 with 'rend' number of points
mylist1 = zeros(1, length(r1st));  % Preallocate a list of zeros for storing results of sol(1)
mylist2 = zeros(1, length(r1st));  % Preallocate a list of zeros for storing results of sol(2)
mylist3 = zeros(1, length(r1st));  % Preallocate a list of zeros for storing results of sol(3)

% Define components of variables (K1 to K8) appearing in polynomial coefficients
K1 = alpha + mu; K2 = gamma + mu; K3 = (Omega.*mu./(1 - Omega)) + mu; K4 = mu; K5 =  (Omega.*mu./(1 - Omega)) + alpha + mu; 
K6 = gamma + rho.*(Omega.*mu./(1 - Omega)) + mu; K7 = alpha + mu; K8 = gamma + mu;  
R_h = p_h.*c.*(eta.*(Omega.*mu./(1 - Omega)) + K4)/(K3.*K4) %HIV Control Reproduction number
lambda_h = mu.* (R_h - 1); %Force of infection for HIV at the HIV endemic regime

R_mh = beta_11.*((alpha.*Lambda)./((lambda_h + mu).*K1.*K2) + (sigma.*alpha.*Lambda.*lambda_h)./((lambda_h + mu).*K3.*K5.*K6) ... 
       + (sigma.*alpha.*(Omega.*mu./(1 - Omega)).*Lambda.*lambda_h)./((lambda_h + mu).*K3.*K5.*K7.*K8) ...
       + (sigma.*alpha.*rho.*(Omega.*mu./(1 - Omega)).*Lambda.*lambda_h)/((lambda_h + mu).*K3.*K5.*K6.*K8)...
       + (sigmaT.*alpha.*(Omega.*mu./(1 - Omega)).*Lambda.*lambda_h)./((lambda_h+ mu).*K3.*K4.*K7.*K8)); %InvR0 for fixed values of the paramters

% Loop over each value in r1st to calculate polynomial coefficients and solve
for i = 1:rend
    R0 = r1st(i);  % Here, we create a loop to vary Invasion R0 within the range: linspace(0,1.1,rend);

    % The coefficients of the cubic polynomial (D1, D2, D3, D4)
    D1 = K1.*K2.*K5.*K6.*K7.*K8.*sigma.*sigmaT;
    D2 = K1.*K2.*K5.*K6.*K7.*K8.*(sigma.*sigmaT.*(lambda_h + mu) + (K4.*sigma + K3.*sigmaT)) - alpha.*beta_11.*Lambda.*K5.*K6.*K7.*K8.*sigma.*sigmaT;
    D3 = K1.*K2.*K3.*K4.*K5.*K6.*K7.*K8 + K1.*K2.*K5.*K6.*K7.*K8.*mu.*R_h.*(K4.*sigma + K3.*sigmaT) ...
        - (alpha.*beta_11.*Lambda.*(K4.*sigma + K3.*sigmaT).*K5.*K6.*K7.*K8 ...
        + alpha.*beta_11.*Lambda.*mu.*(R_h - 1).*sigma.*sigmaT.*K1.*K2.*(K7.*rho.*(Omega.*mu./(1 - Omega)) + K6.*(Omega.*mu./(1 - Omega)) + K7.*K8));
    D4 = K1.*K2.*K3.*K4.*K5.*K6.*K7.*K8.*mu.*R_h.*(1 - R0);
    
    mypoly = [D1, D2, D3, D4];  % Form the polynomial from the coefficients
    sol = roots(mypoly);  % Calculate the roots of the polynomial
    
    mylist1(i) = sol(1);  % Store the first root in mylist1
    mylist2(i) = sol(2);  % Store the second root in mylist2
    mylist3(i) = sol(3);  % Store the third root in mylist3
    
    % Ensure that the second and third roots are real and non-negative
    if imag(mylist2(i)) ~= 0
        mylist2(i) = 0;
    else
        mylist2(i) = max(0, mylist2(i));
    end
    if imag(mylist3(i)) ~= 0
        mylist3(i) = 0;
    else
        mylist3(i) = max(0, mylist3(i));
    end
end

% Find the indices where the values in mylist2 and mylist3 are non-zero
v1stpos = find(mylist2 ~= 0);
index1 = min(v1stpos);  % Index of the first non-zero value in mylist2
index2 = max(v1stpos);  % Index of the last non-zero value in mylist2
m1stpos = find(mylist3 ~= 0);
index3 = min(m1stpos);  % Index of the first non-zero value in mylist3
index4 = max(m1stpos);  % Index of the last non-zero value in mylist3

% Define lines for plotting purposes
s = linspace(0, 1, rend);  % Linearly spaced vector for plotting a horizontal line at y=0
z = linspace(1, 1.1, rend);  % Linearly spaced vector for plotting another line at y=0
y1 = 0.*s;  % Define a horizontal line at y=0 for the first range
y2 = 0.*z;  % Define a horizontal line at y=0 for the second range

% Plot the results
hold on
plot(r1st(index1:index2), mylist2(index1:index2), 'b--', r1st(index3:index4), mylist3(index3:index4), 'r--', 'LineWidth', 2);
plot(s, y1, 'b-', 'LineWidth', 2);
plot(z, y2, 'r-', 'LineWidth', 2);

% Label the axes with LaTeX formatting
xlabel('$\bar{R}^{mh}_{c}$', 'Interpreter', 'latex');
ylabel('$\lambda^{e}_{m}$', 'Interpreter', 'latex');

% Add a box around the plot and hold for additional plots
box on;
hold on;

% Add a text box annotation to the plot with specific formatting
t = annotation('textbox', 'String', '\sigma = 1.5, \sigma^T = 1.5');
t.FontSize = 15;  % Set the font size to 15
t.FontWeight = 'bold';  % Set the font weight to bold
t.Rotation = 25;  % Rotate the textbox by 30 degrees
t.EdgeColor = 'none';  % Remove the box around the text
