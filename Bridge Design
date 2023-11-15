P = 400;
x_train = 110;

x_wheels = [52 228 392 568 732 908] + x_train;

L = 1200;
x = 0:L;

P_wheel = -P/6

P_x = zeros(1, L+1);
P_x(end) = -P_wheel*sum(x_wheels)/1200;
P_x(1) = -(P_wheel*length(x_wheels) + P_x(end));
P_x(x_wheels + 1) = P_wheel

    
SFD = cumsum(P_x)
BMD = cumsum(SFD);
BMD = [0 BMD(1:end-1)]

disp('Sign convention used for SFD is up is positive.')
figure
hold on
plot(x, SFD, 'b-')
plot(x, zeros(1,1201), 'k-')
title('SFD (N)')
ylabel('- Force (N) +')
xlabel('Position of Bridge (mm)')
hold off

disp('Standard sign convention was used. For BMD is down is positive.')
figure
hold on
title('BMD (N*mm)')
plot(x,BMD,'b-')
plot(x, zeros(1,1201), 'k-')
xlabel('Position of Bridge (mm)')
ylabel('+ Moment (N*mm) -')
set(gca, 'YDir','reverse')
hold off
SFD and BMD Envelope
SFD_max = SFD
SFD_min = SFD
BMD_env = BMD
for i = 0:240
    x_wheels = [52 228 392 568 732 908] + i;
    
    P_x = zeros(1, L+1);
    P_x(end) = -P_wheel*sum(x_wheels)/1200;
    P_x(1) = -(P_wheel*length(x_wheels) + P_x(end));
    P_x(x_wheels + 1) = P_wheel;
    
    SFD_i = cumsum(P_x);
    BMD_i = cumsum(SFD_i);
    BMD_i = [0 BMD_i(1:end-1)];

    SFD_max = max(SFD_i,SFD_max);
    SFD_min = min(SFD_i,SFD_min);

    BMD_env = max(BMD_i, BMD_env);
end

SFD_env = [SFD_max(1:600) SFD_min(601:1201)]
BMD_env
BMD_max_loc = 0;
for i = 0:240
    x_wheels = [52 228 392 568 732 908] + i;
    
    P_x = zeros(1, L+1);
    P_x(end) = -P_wheel*sum(x_wheels)/1200;
    P_x(1) = -(P_wheel*length(x_wheels) + P_x(end));
    P_x(x_wheels + 1) = P_wheel;
    
    SFD = cumsum(P_x);
    BMD = cumsum(SFD);
    BMD = [0 BMD(1:end-1)];

    if max(BMD) == max(BMD_env)
        BMD_max_loc = i;
        break;
    end
end
BMD_max_loc
max_V = max(SFD_env)
max_M = max(BMD_env)

figure
hold on
plot(x,SFD_env,'b-')
plot(x, zeros(1,1201), 'k-')
hold off

figure
hold on
plot(x, BMD_env, 'r-')
plot(x, zeros(1,1201), 'k-')
set(gca, 'YDir','reverse')
hold off

Calculating I and y_bar
num_param = 3
dim_beam_wide = [75 1.27 6.27 100] %"unique dimentions"
dim_beam_h = [1.27 117.5 1.27 1.27*3]
A_total = 813*1016;
b_1 = dim_beam_wide(1)-1.27;
a = zeros(1,1201);
a(1:300) = 100;
a(301:450) = 150;
a(451:750) = 300;
a(751:900) = 150;
a(901:1201) = 100;
dim_A = dim_beam_h.*dim_beam_wide;
A_used = sum(dim_A)*1250/1.27+10*dim_beam_wide(1)*dim_beam_h(2);
d_A = A_total - A_used
for i = 2:3 %change this to change repeated dims
    dim_A(i) = dim_A(i)*2;
end
y_i = cumsum(dim_beam_h);
for i = 1:4
    y_i(i) = y_i(i) - dim_beam_h(i)/2;
end

y_bar = sum(y_i.*dim_A)/sum(dim_A)
y_bot = y_bar
y_top = sum(dim_beam_h)-y_bar

I = sum((dim_A.*(dim_beam_h.^2))./12 + ((y_i - y_bar).^2).*dim_A)
Q Calcs
Qcent = dim_A(1)*(y_bar-y_i(1)) + (y_bar-dim_beam_h(1))^2*(dim_beam_wide(2));
Qglue = dim_A(4)*(y_i(4)-y_bar);
Applied stress
b = dim_beam_wide(2)*2;
b_g = dim_beam_wide(3)*2;
S_top = BMD_env*(y_top)/I;
S_bot = BMD_env*(y_bot)/I;
T_cent = (SFD_env*Qcent)/(I*b);
T_glue = (SFD_env*Qglue)/(I*b_g);
Material and Thin Plate Buckling Capacities
E = 4000;
mu = 0.2;
S_tens  = 30;
S_comp  =  6;
T_max   =  4;
T_gmax  =  2;
S_buck1 =  StressBucking(dim_beam_h(4),b_1,mu,E,1);
S_buck2 =  StressBucking(dim_beam_h(4),(dim_beam_wide(4)-b_1)/2,mu,E,2);
S_buck3 =  StressBucking(dim_beam_wide(2),(y_i(3)-y_bar),mu,E,3);
T_buck  =  ShearBuckling(dim_beam_wide(2), a, y_i(3)-y_i(1), mu, E);
6. FOS
FOS_tens  =  abs(S_tens ./ S_bot);
FOS_comp  =  abs(S_comp ./ S_top);
FOS_shear =  abs(T_max ./ T_cent);
FOS_glue  =  abs(T_gmax ./ T_glue);
FOS_buck1 =  abs(S_buck1 ./ S_top);
FOS_buck2 =  abs(S_buck2 ./ S_top);
FOS_buck3 =  abs(S_buck3 ./ S_top);
FOS_buckV =  abs(T_buck ./ T_cent);

min_FOS_tens  = min(FOS_tens)
min_FOS_comp  = min(FOS_comp)
min_FOS_shear = min(FOS_shear)
min_FOS_glue  = min(FOS_glue)
min_FOS_buck1 = min(FOS_buck1)
min_FOS_buck2 = min(FOS_buck2)
min_FOS_buck3 = min(FOS_buck3)
min_FOS_buckV = min(FOS_buckV)
7. Min FOS and the failure load Pfail 
minFOS      =  min(abs([FOS_tens FOS_comp FOS_shear FOS_glue FOS_buck1 FOS_buck2 FOS_buck3 FOS_buckV]))
Pf          =  minFOS * P
8. Vfail and Mfail 
Mf_tens  =  FOS_tens .* BMD_env
Mf_comp  =  FOS_comp .* BMD_env
Vf_shear =  FOS_shear.* SFD_env
Vf_glue  =  FOS_glue .* SFD_env
Mf_buck1 =  FOS_buck1.* BMD_env
Mf_buck2 =  FOS_buck2.* BMD_env
Mf_buck3 =  FOS_buck3.* BMD_env
Vf_buckV =  FOS_buckV.* SFD_env
9. Output plots of Vfail and Mfail 
subplot(2,3,1) 
hold on; grid on; grid minor; 
plot(x, abs(Vf_shear), 'r')  
plot(x,-abs(Vf_shear),'r-')
plot(x, SFD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Matboard Shear Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)') 
hold off

subplot(2,3,2)
hold on; grid on; grid minor;
plot(x, abs(Vf_glue), 'r')  
plot(x,-abs(Vf_glue),'r-')
title('Shear Force Diagram vs Shear Force Capacities')
plot(x, SFD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Glue Shear Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)') 
hold off

subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, abs(Vf_buckV), 'r')  
plot(x,-abs(Vf_buckV),'r-')
plot(x, SFD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Matboard Shear Buckling Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Shear Force (N)') 
hold off

subplot(2,3,4)
hold on; grid on; grid minor;
plot(x, abs(Mf_tens), 'r')  
plot(x, abs(Mf_comp),'b-')
plot(x, BMD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Matboard Tension Failure','Matboard Compression Failure') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (N*mm)') 
set(gca, 'YDir','reverse')
hold off

subplot(2,3,5)
hold on; grid on; grid minor;
title('Bending Moment Diagram vs Bending Moment Capacities')
plot(x, abs(Mf_buck1), 'r')  
plot(x, abs(Mf_buck2),'b-')
plot(x, BMD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Matboard Buckling Failure, Top Flange-Mid','Matboard Compression Failure, Top Flange-Sides') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (N*mm)') 
set(gca, 'YDir','reverse')
hold off

subplot(2,3,6)
hold on; grid on; grid minor;
title('Bending Moment Diagram vs Bending Moment Capacities')
plot(x, abs(Mf_buck3), 'r')  
plot(x, BMD_env, 'k'); 
plot([0, L], [0, 0], 'k', 'LineWidth', 2) 
legend('Matboard Buckling Failure, Webs') 
xlabel('Distance along bridge (mm)') 
ylabel('Bending Moment (N*mm)') 
set(gca, 'YDir','reverse')
hold off

Functions for buckling fail
function [S_buck] = StressBucking(t,b,mu, E, type)
    if type == 1
        K = 4;
    elseif type == 2
        K = 0.425;
    else
        K = 6;
    end
    S_buck = (K*pi^2*E)/(12*(1-mu^2))*(t/b)^2;
end

function [T_buck] = ShearBuckling(t,a,b,mu, E)
    T_buck = (5*pi^2*E)/(12*(1-mu^2)).*((t./a).^2+(t/b)^2);
end



