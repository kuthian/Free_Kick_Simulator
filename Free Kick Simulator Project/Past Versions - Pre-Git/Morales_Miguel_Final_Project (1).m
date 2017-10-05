clc

%%%%%%%%%%%%%%%%%%%
%Defining Variables


%user prompt for initial variables
prompt = {'Distance from the goal (m): ','Rise Angle (deg): ','Initial Velocity (m/s): ','Spin (rotation/s): ', 'Shot Angle (deg): '};
def={'16','45','14','20','8'};  
answer = inputdlg(prompt,'Input',1,def);
[a1, a2, a3, a4, a5] = answer{:};
net = str2num(a1); 
theta_deg = str2num(a2);
v_0 = str2num(a3);
rot = str2num(a4);   
phi_deg = str2num(a5);

clf
hold on

t_step= 0.01; %time step between each calculation
num_step = 5/t_step; %the number of steps
% theta_deg = 9.1; %angle of the kick in degrees
theta_rad = theta_deg*3.14159*(1/180); %angle of the kick in radians
% phi_deg = -7; %angle of the kick in degrees
phi_rad = phi_deg*3.14159*(1/180); %angle of the kick in radians
w = rot*pi*2; %angular velocity on ball
s0= 4e-4; %average drag coeff
a_y = -9.81; %initial acceleration in the y direction
a_x = 0; %initial acceleration in the x direction
a_z = 0; %initial acceleration in the z direction
% v_0 = 40; %initial magnitude of the velocity
v_0x = v_0*cos(theta_rad); %initial velocity in the x direction
v_0y = v_0*sin(theta_rad); %initial velocity in the y direction
v_0z = v_0x*tan(phi_rad); %initial velocity in the z direction
d_0y = 0; %initial displacement in the y direction
d_0x = 0; %initial displacement in the x direction
d_0z = 11; %initial displacement in the z direction
r = 0.109; %radius of the soccer ball
p =1.6; %average density of the soccer ball
C = 0.15; %drag coefficient
A = 3.14159*(r^2); %silhouette area of soccer ball
D = p*C*A*0.5; %constant
m = 0.5; %mass of soccer ball
react = 0.20; 
react_step= 0.20/t_step;
% net=30;

%%%%%%%%%%%%%%%%%%%%
%time array
t_vals = 0:num_step;

%%%%%%%%%%%%%%%%%%%%
%acceleration arrays
ay_vals = 0:num_step;
ay_vals(1) = a_y;

ax_vals = 0:num_step;
ax_vals(1) = a_x;

az_vals = 0:num_step;
az_vals(1) = a_z;

%%%%%%%%%%%%%%%%%%%%
%velocity arrays
v_vals = 0:num_step;

vy_vals = 0:num_step;
vx_vals = 0:num_step; 

vy_vals(1) = v_0y;
vx_vals(1) = v_0x;

vz_vals = 0:num_step; 
vz_vals(1) = v_0z;

%%%%%%%%%%%%%%%%%%%%
%position arrays
d_x = 0:num_step;
d_x(1) = d_0x;

d_y = 0:num_step;
d_y(1) = d_0y;

d_z = 0:num_step;
d_z(1) = d_0z;

%%%%%%%%%%%%%%%%%%%%
%objects

%Soccer Net Display
line([net net], [7.35 7.35] , [0 2.44] ,'LineWidth',4)
line([net net], [14.65 14.65] , [0 2.44] ,'LineWidth',4)
line([net net], [7.35 14.65] , [2.44 2.44] ,'LineWidth',4)
line([net net+2], [7.35 7.35] , [0 0] ,'LineWidth',4)
line([net net+2], [14.65 14.65] , [0 0] ,'LineWidth',4)
line([net+2 net], [7.35 7.35] , [0 2.44] ,'LineWidth',4)
line([net+2 net], [14.65 14.65] , [0 2.44] ,'LineWidth',4)
line([net+2 net+2], [7.35 14.65] , [0 0] ,'LineWidth',4)


%Grid vertical
for i=0:23;
line([net+2 net], [7.85+i*(0.3) 7.85+i*(0.3)] , [0 2.44])
end

%Grid horizontal
for i=0:7;
line([net+2-i*(0.3) net+2-i*(0.3)], [7.35 14.65] , [0+i*(0.3) 0+i*(0.3)])
%Grid side horizontal
line([net net+2-i*(0.3)], [7.35 7.35] , [0+i*(0.3) 0+i*(0.3)])
line([net net+2-i*(0.3)], [14.65 14.65] , [0+i*(0.3) 0+i*(0.3)])

line([net+i*(0.3) net+i*(0.3)], [7.35 7.35] , [0 2.44-i*(0.3)])
line([net+i*(0.3) net+i*(0.3)], [14.65 14.65] , [0 2.44-i*(0.3)])
end


%wall of players Display
if net >15
line([9.1 9.1], [d_0z-1 d_0z-1] , [0 1.74] ,'LineWidth',2)
line([9.1 9.1], [d_0z+1 d_0z+1] , [0 1.74] ,'LineWidth',2)
line([9.1 9.1], [d_0z-1 d_0z+1] , [1.74 1.74] ,'LineWidth',2)
end

%%initial goalkeeper%%
left_lim= d_0z - 0.25;
right_lim = d_0z + 0.25; 

a = line([net net], [left_lim left_lim] , [0 1.74],'Color','r' ,'LineWidth',2);
b = line([net net], [right_lim right_lim] , [0 1.74],'Color','r' ,'LineWidth',2);
c = line([net net], [left_lim right_lim] , [1.74 1.74],'Color','r' ,'LineWidth',2);
% bicep
% d = line([net net], [left_lim+0.5 right_lim+0.3] , [1.44 1.44],'Color','r' ,'LineWidth',2);
% tricep
% e = line([net net], [left_lim+0.5 right_lim+0.5] , [1.34 1.34],'Color','r' ,'LineWidth',2);
% 
% forearm
% f = line([net net], [left_lim+1 right_lim+0.5] , [1.34 1.74],'Color','r' ,'LineWidth',2);
% g = line([net net], [left_lim+0.8 right_lim+0.3] , [1.44 1.74],'Color','r' ,'LineWidth',2);
% 
% other side
% hand
% h = line([net net], [left_lim - 0.5 right_lim - 0.8] , [1.74 1.74],'Color','r' ,'LineWidth',2);
% 
% bicep
% i = line([net net], [left_lim-0.3 right_lim-0.5] , [1.44 1.44],'Color','r' ,'LineWidth',2);
% tricep
% j = line([net net], [left_lim-0.5 right_lim-0.5] , [1.34 1.34],'Color','r' ,'LineWidth',2);
% 
% forearm
% k = line([net net], [left_lim-0.5 right_lim-1] , [1.34 1.74],'Color','r' ,'LineWidth',2);
% l = line([net net], [left_lim-0.3 right_lim-0.8] , [1.44 1.74],'Color','r' ,'LineWidth',2);
% hand
% h = line([net net], [left_lim - 0.5 right_lim - 0.8] , [1.74 1.74],'Color','r' ,'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%
%Graph specifics
grid on
axis square
view(3);
title('Soccer Ball Trajectory')
xlabel('x-position (m)')
ylabel('y-position (m)')
zlabel('z-position (m)')
axislims = axis;
axis([0 31 axislims(3) axislims(4)]);
axis([-1 net 0 20 0 10 ])

check = 0; 
check_wall = 0;
check_keeper = 0;


for (i=1: num_step)    
    
    %array calculations
    t_vals(i+1)= t_vals(i)+t_step;
    v_vals(i) = sqrt((vx_vals(i)^2)+(vx_vals(i)^2));
    ax_vals(i) = -(D/m)*v_vals(i)*vx_vals(i);
    ay_vals(i) = a_y-(D/m)*v_vals(i)*vy_vals(i);
    az_vals(i) = a_z-(s0*w/m)*vx_vals(i); 
    vx_vals(i+1)= vx_vals(i)+ax_vals(i)*t_step;
    vy_vals(i+1)= vy_vals(i)+ay_vals(i)*t_step;
    vz_vals(i+1)= vz_vals(i)+az_vals(i)*t_step;
    d_x(i+1)= d_x(i) + vx_vals(i+1)*t_step + 0.5*ax_vals(i)*(t_step^2);
    d_y(i+1)= d_y(i) + vy_vals(i+1)*t_step + 0.5*ay_vals(i)*(t_step^2);
    d_z(i+1)= d_z(i) + vz_vals(i+1)*t_step + 0.5*az_vals(i)*(t_step^2);
   
    %%plot animation%%
    if i>1
        delete(h);
    end    
    figure(1);
    h = plot3(d_x(i+1),d_z(i+1),d_y(i+1),'ok' );
    
    %Goalkeeper AI%%   
    if (right_lim > 7.35 && left_lim < 14.65 && i > react_step )
    left_lim=d_0z+0.25-d_0z+d_z(i-react_step);
    right_lim =d_0z-0.25-d_0z+d_z(i-react_step); 
    
        if i>react_step
        delete(a,b,c);
        end  
        
    a = line([net net], [left_lim left_lim] , [0 1.74],'Color','r' ,'LineWidth',2);
    b = line([net net], [right_lim right_lim] , [0 1.74],'Color','r' ,'LineWidth',2);
    c = line([net net], [left_lim right_lim] , [1.74 1.74],'Color','r' ,'LineWidth',2);
    end

    pause(0.01)
    
    %check to see if the ball has hit the wall
    if net >15
    if (check_wall == 0 && d_x(i+1) > 9.1 && d_y(i+1) < 1.74 && d_z(i+1) < (d_0z+1) && d_z(i+1) > (d_0z-1))
        final =i+1;
        disp('The ball has hit the wall of players')
        break;
    else if check_wall ==0 && d_x(i+1) > 9.1
       check_wall = 1;
        end
       
    end
    end
    
    %check to see if the ball has been caught by goalkeeper
    
    if (check_keeper == 0 && d_x(i+1) > net && d_y(i+1) < 1.85 && d_z(i+1) < (left_lim+0.8) && d_z(i+1) > (right_lim-0.8))
        final =i+1;
        disp('The ball has been caught by goalkeeper.')
        plot3(d_x(final),d_z(final),d_y(final),'or','LineWidth',10)
        break;
    else if check_keeper ==0 && d_x(i+1) > net
       check_keeper = 1;
        end
       
    end
    
    %check to see if the ball has gone into the net
    if (check == 0 && d_x(i+1) > net && d_y(i+1) < 2.44 && d_z(i+1) < 14.65 && d_z(i+1) > 7.35)
        final =i+1;
        disp('GOAL!')
        plot3(d_x(final),d_z(final),d_y(final),'oc','LineWidth',10)
        break;
    else if check ==0 && d_x(i+1) > net
       check = 1;
        end       
        
    end
    
    %check to see if the ball has gone over the wall and net and hit the ground
    if (d_y(i+1) < 0)
        final = i+1;
        disp('The ball has gone over the net and hit the ground.')
        plot3(d_x(final),d_z(final),d_y(final),'xr','LineWidth',20)
        break;
    end
    
   
end

%calculaion of final velocity
v_final = sqrt((vx_vals(final)^2)+(vx_vals(final)^2));

%trajectory trail
plot3(d_x(1:final),d_z(1:final),d_y(1:final),'.k' );

%%%%%%%%%%
%display of final variables
disp(['The final velocity is ', num2str(v_final),'. The final position at time ',num2str(t_vals(final)), 's is [', num2str(d_x(final)),',',num2str(d_z(final)),',', num2str(d_y(final)),'].'])