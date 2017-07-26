
format long;
L = 0.003;
C = 1e-10;
R = 300;    % for underdamped
%R = 11000;  %for overdamped
%R = sqrt(4*L/C);    %for critically damped
xf = 1e-4;  %underdamped
%xf = 1e-5; % for overdamped and critically damped
x0 = 0;
n = 500;
f11 = @(x,y1,y2) y2;
f12 = @(x,y1,y2) -y1/(L*C) - (R*y2)/L;
h = (xf-x0)/n;
h24 = h / 24;

x(1) = x0;
y11(1)=3;
y12(1)= 0;

%for overdamp
%fdexact = @(x) 18*(-1.6667e+6)*exp((-1.6667e+6)*x) + -15*(-2e+6)*exp((-2e+6)*x);   
%fexact = @(x) 18*exp((-1.6667e+6)*x) + -15*exp((-2e+6)*x);
%for critically damped
%fexact = @(x) 3*exp((-1.8257e+6)*x) + (5.472e+6)*x*exp((-1.8257e+6)*x);
%fdexact = @(x) 3*(-1.8257e+6)*exp((-1.8257e+6)*x) + (5.472e+6)*(-1.8257e+6)*x*exp((-1.8257e+6)*x) + (5.472e+6)*exp((-1.8257e+6)*x);
%for under damped
fexact = @(x) 3*(exp((-5e+4)*x))*cos((1.825e+6)*x) + 0.082192*(exp((-5e+4)*x))*sin((1.825e+6)*x);
fdexact = @(x) 3*(-5e+4)*(exp((-5e+4)*x))*cos((1.825e+6)*x) + 0.082192*(-5e+4)*(exp((-5e+4)*x))*sin((1.825e+6)*x) + 3*(1.825e+6)*(exp((-5e+4)*x))*sin((1.825e+6)*x) + 0.082192*(1.825e+6)*(exp((-5e+4)*x))*cos((1.825e+6)*x);

for i=1:3 
    x(i+1) = x(i) + h;
    k_1 = f11(x(i),y11(i),y12(i));
    L_1 = f12(x(i),y11(i),y12(i));
    k_2 = f11(x(i)+0.5*h,y11(i)+0.5*h*k_1,y12(i)+0.5*h*L_1);
    L_2 = f12(x(i)+0.5*h,y11(i)+0.5*h*k_1,y12(i)+0.5*h*L_1);
    k_3 = f11((x(i)+0.5*h),(y11(i)+0.5*h*k_2),(y12(i)+0.5*h*L_2));
    L_3 = f12((x(i)+0.5*h),(y11(i)+0.5*h*k_2),(y12(i)+0.5*h*L_2));
    k_4 = f11((x(i)+h),(y11(i)+k_3*h),(y12(i)+L_3*h)); % Corrected        
    L_4 = f12((x(i)+h),(y11(i)+k_3*h),(y12(i)+L_3*h));

    y11(i+1) = y11(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    y12(i+1) = y12(i) + (1/6)*(L_1+2*L_2+2*L_3+L_4)*h;  % main equation
end

for i = 4 : n % main phase
    x(i+1) = x(i) + h;
    y11(i+1) = y11(i) + (55 * f11(x(i),y11(i),y12(i)) - 59 * f11(x(i-1),y11(i-1),y12(i-1)) + 37 * f11(x(i-2),y11(i-2),y12(i-2)) - 9 * f11(x(i-3),y11(i-3),y12(i-3))) * h24; % predictor
    y12(i+1) = y12(i) + (55 * f12(x(i),y11(i),y12(i)) - 59 * f12(x(i-1),y11(i-1),y12(i-1)) + 37 * f12(x(i-2),y11(i-2),y12(i-2)) - 9 * f12(x(i-3),y11(i-3),y12(i-3))) * h24;
    
    y11(i+1) = y11(i) + (9 * f11(x(i+1), y11(i+1),y12(i+1)) + 19 * f11(x(i),y11(i),y12(i)) - 5 * f11(x(i-1),y11(i-1),y12(i-1)) + f11(x(i-2),y11(i-2),y12(i-2))) * h24; % corrector
    y12(i+1) = y12(i) + (9 * f12(x(i+1), y11(i+1),y12(i+1)) + 19 * f12(x(i),y11(i),y12(i)) - 5 * f12(x(i-1),y11(i-1),y12(i-1)) + f12(x(i-2),y11(i-2),y12(i-2))) * h24;
end;
for i=1:n+1 %exact solution
    yexact(i) = feval(fexact,(x(i)));
    yd(i) = feval(fdexact,x(i));
end
% voltage Vs time
plot(x,y11,'-',x,yexact,'-o'); grid on; title('Using adam4 method');
xlabel('Voltage ---->');ylabel('Time ---->'); legend('y-adm4','y-exact');
% voltage gradient vs voltage
%plot(y11,y12,'-',yexact,yd,'-o'); grid on; title('Using adm4 method');
%xlabel('Voltage ---->');ylabel('dV/dt'); legend('y-adm4','y-exact');
