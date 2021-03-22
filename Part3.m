% ELEC 4700
% Name: Oritseserundede Eda 
% Student Number: 100993421
% Assignment 3 
% Part 3
% In this part, we use the coupled simulations for the device and 
% trajectory investigations  

clc
clear
set(0,'DefaultFigureWindowStyle','docked')

global C
global Vtotal Vx Vy x y
global Ecount
Ecount =1000; 
C.mo = 9.10938215e-31;
C.k = 1.3806504e-23; 
qcharge = -1.60217662e-19;
temp =300;
eff_mass = 0.26*C.mo;
length = 200e-9;
width = 100e-9; 
therm_volt = sqrt((2*C.k*temp)/eff_mass);
t_step = 10e-15; 
frame = 100*t_step; 
x = zeros(Ecount, 2);  
y = zeros(Ecount, 2);  
temp = zeros(1,2); 
Time = 0;
visible_ecount = 50; 
tmn = 0.2e-12;
PScat = 1 - exp(-t_step/tmn);
V_Histogram = zeros(Ecount, 1);
b_x = [80e-9 80e-9 120e-9 120e-9 80e-9];

b_y1 = [100e-9 60e-9 60e-9 100e-9 100e-9];
b_y2 = [40e-9 0 0 40e-9 40e-9]; 
spec = true; 
inside_box = true; 
s_map = 10e-9;
d_map = zeros(width/s_map, length/s_map);
t_map = zeros(width/s_map, length/s_map);

wid_x = 30;
len_y = 20;
change_x = length/wid_x;
change_y = width/len_y;
outside_conduct = 1;
inside_conduct = 01e-2;
conductivity = zeros(wid_x,len_y);
G = sparse (wid_x*len_y, wid_x*len_y);
V = zeros(1, wid_x*len_y);
volt_x = 0.1;

for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
        if (i > (0.3*wid_x) || i < (0.6*wid_x)) && (j > (0.6*len_y) || j < (0.3*len_y))
            conductivity(i,j) = inside_conduct;
        else
            conductivity(i,j) = outside_conduct;
        end
        
    end
end
for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
        if (i == 1)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            V(n) = volt_x;
            G(n,n) = 1;
        elseif (i == wid_x)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            V(n) = 0;
            G(n,n) =1;
        elseif (j == 1)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j+1))/2);
            G(n, nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G(n, nyp) = ((conductivity(i,j) + conductivity(i,j+1))/2);
        elseif (j == len_y)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j-1))/2);
            G(n,nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G(n,nym) = ((conductivity(i,j) + conductivity(i,j-1))/2);
        else
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j-1))/2) - ((conductivity(i,j) + conductivity(i,j+1))/2);
            G(n,nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G(n,nym) = ((conductivity(i,j) + conductivity(i,j-1))/2);
            G(n,nyp) = ((conductivity(i,j) + conductivity(i,j+1))/2);
        end
    end
end

soln = G\V';
surf = zeros(wid_x,len_y);
for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1) * len_y;
        nxp = j + ((i+1) - 1) * len_y;
        nym = (j-1) + (i - 1) * len_y;
        nyp = (j+1) + (i - 1) * len_y;
        surf(i,j) = soln(n);
    end
end
[E_x, E_y] = gradient(-surf);
Force_x = qcharge*E_x;
Force_y = qcharge*E_y;
Acceleration_x = Force_x /eff_mass;
Acceleration_y = Force_y /eff_mass;

for i = 1:Ecount
    x(i,1) = rand()*200e-9;
    y(i,1) = rand()*100e-9;
    inside_box = true;
    while inside_box == true
        if (x(i) >= 40e-9 && x(i) <= 120e-9) && (y(i) >= 60e-9 ||...
                y(i) <= 40e-9)
            x(i,1) = rand * 200e-9;
            y(i,1) = rand * 100e-9;
        else
            inside_box = false;
        end
    end
    
end

for i = 1:Ecount
    
Vx(1:Ecount) = therm_volt * randn;
Vy(1:Ecount) = therm_volt * randn;
end

figure(8)
subplot(2,1,1);
plot(b_x, b_y1, b_x, b_y2)
axis([0 length 0 width]);
title('Question 3');
xlabel('x');
ylabel('y');
hold on;

while Time < frame
    subplot(2,1,1)
    for j = 1:Ecount
        leaking = true;
        if PScat> rand
                Vx(j) = therm_volt * randn;
                Vy(j) = therm_volt * randn;
        end
        x_index = round((x(j,2)/length) * 30);
        y_index = round((y(j,2)/width)*20);
        if x_index < 1
            x_index = 1;
        elseif x_index > 30 
                x_index = 30;
        end
        if y_index < 1
            y_index = 1;
        elseif y_index > 20
            y_index = 20;
        end
        
        Vx(j) =  Vx(j) + Acceleration_x(x_index,y_index)*t_step;
        Vy(j) =  Vy(j) + Acceleration_y(x_index,y_index)*t_step;
        x(j,2) = x(j,1);
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (t_step * Vx(j));
        y(j,1) = y(j,1) + (t_step * Vy(j));
        
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) >= 60e-9
            
                if y(j,2) < 60e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 60e-9;
                    y(j,2) = 60e-9;
                    
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                    
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
                
            if spec == true
                
                x(j,1) = x(j,2) + Vx(j)*t_step;
                y(j,1) = y(j,2) + Vy(j)*t_step;
            else
                
             Vx(j) = therm_volt * randn;
             Vy(j) = therm_volt * randn;
             
             while leaking == true
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) < 60e-9 && Vy(j) >= 0)      
                     Vx(j) = therm_volt * randn;
                     Vy(j) = therm_volt * randn;
                 else
                     leaking = false;
                 end
             end
             x(j,1) = x(j,2) + Vx(j)*t_step;
             y(j,1) = y(j,2) + Vy(j)*t_step;
            end
        end
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) <= 40e-9
                if y(j,2) > 40e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 40e-9;
                    y(j,2) = 40e-9; 
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
            if spec == true
                x(j,1) = x(j,2) + Vx(j)*t_step;
                y(j,1) = y(j,2) + Vy(j)*t_step;
            else
             Vx(j) = therm_volt * randn;
             Vy(j) = therm_volt * randn;
             while leaking == true
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) > 40e-9 && Vy(j) <= 0)
                     Vx(j) = therm_volt * randn;
                     Vy(j) = therm_volt * randn;
                 else
                     leaking = false;
                 end
             end
             x(j,1) = x(j,2) + Vx(j)*t_step;
             y(j,1) = y(j,2) + Vy(j)*t_step;
            end
        end
        if x(j,1) > length
            x(j,2) = 0;
            x(j,1) = t_step * Vx(j);
        end
        if x(j,1) < 0
            x(j,2) = length;
            x(j,1) = x(j,2) + (t_step * Vx(j));
        end
        if y(j,1) > width || y(j,1) < 0
            Vy(j) = -Vy(j);
        end
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
        if j < visible_ecount
        plot(XPlot,YPlot);
        end
        
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);    
    end

    average_temp = temp(1,2)/Ecount;
    temp_plot = [temp(1,1) average_temp];
    time_plot = [(Time - t_step) Time];
    subplot(2,1,2);
    plot(time_plot, temp_plot);
    temp(1,1) = average_temp;
    average_temp = 0;
    temp(1,2) = 0;
    pause(1e-19)
    Time = Time + t_step;
end 
for i = 1:(length/s_map)
    for j = 1:(width/s_map)
        for m = 1:Ecount
            if(x(m,1) > s_map*(i -1)) && ...
                    (x(m,1) < s_map*(i)) && ...
                    (y(m,1) > s_map*(j - 1)) && ...
                    (y(m,1) < s_map*(j))
             
                Vtotal(m) = sqrt(Vx(m)^2 + Vy(m)^2);
                
                d_map(j, i) = d_map(j, i) + 1;
                t_map(j, i) = t_map(j,i) + ...
                    (eff_mass*Vtotal(m)^2)/(2*C.k);
            end
            t_map(j,i) = t_map(j,i)/d_map(j,i);
        end
    end
end
figure(9)
imagesc(d_map)
title('The Density-Mapping of all electrons in the frame')
xlabel('x');
ylabel('y');
set(gca, 'Ydir', 'Normal')
title('The Electron Count')
