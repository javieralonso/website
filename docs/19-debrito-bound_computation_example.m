clear all 
close all
clc 
 
%% ellipsoid and circle dimensions
a = 10; 
b = 2; 
r = 0.4;
 
M = 400; 
dtheta = 2*pi / M; 
theta_M = (0 : dtheta : 2*pi)'; 


minimal_positive_root = @(delta) (2*(delta + r)^2*(2*a*b + a*(delta + r) + b*(delta + r)))/((a + b)*(a + b + 2*delta + 2*r))-r^2; 
x0 = 0.01; % initial guess must be always positive
delta = fzero(minimal_positive_root,x0)

for i = 1 : M+1 
    theta = theta_M(i); 
    
    %% ellipse coordinates
    x_M(i) = a * cos(theta); 
    y_M(i) = b * sin(theta); 
     
    alpha = a+delta+r; 
    beta = b+delta+r; 

    %% bounding ellipse
    
    x_1_M(i) = alpha * cos(theta); 
    y_1_M(i) = beta * sin(theta); 
    
    %% Minkowsky sum of ellipse (a,b) and circle r
    x_2_M(i) = a*cos(theta) + r*cos(theta)/(sqrt((cos(theta))^2 + (a^2/b^2)*(sin(theta))^2)); 
    y_2_M(i) = b*sin(theta) + r*sin(theta)/(sqrt((b^2/a^2)*(cos(theta))^2 + (sin(theta))^2));    
     

    %% previously used bounding ellipse
    a_3 = a + r; 
    b_3 = b + r; 
    x_3_M(i) = a_3 * cos(theta); 
    y_3_M(i) = b_3 * sin(theta); 
     
    %% circle coordinates
    x_4_M(i) = r * cos(theta); 
    y_4_M(i) = r * sin(theta); 
end 

h=figure; 
hold all; 
box on; 
grid on; 
axis equal; 
plot(x_M, y_M, '-r') 

plot(x_1_M, y_1_M, '-k') 
plot(x_2_M, y_2_M, '-b') 
plot(x_3_M, y_3_M, '-r')
legend("Ellipse","Minimal Bounding ellipse","Minkowski Sum","Bound Ellipse + Circle","Circle")
% legend(h,'off')
plot(x_4_M, y_4_M, '-g') 
% circle(0,0,a+r)


K = randi([floor(M/4/4), ceil(M/4/2)]); 
theta = theta_M(K); 
x = a * cos(theta); 
y = b * sin(theta); 
plot([0,x], [0, y], '-k') 
normal = r*[2*cos(theta)/a; 2*sin(theta)/b] / norm([2*cos(theta)/a; 2*sin(theta)/b]); 
plot([0, normal(1)], [0, normal(2)], '-b') 
plot([x, normal(1)+x], [y, normal(2)+y], '-b')
circle(normal(1)+x, normal(2)+y, r)
hold all
 
% theta_T = theta_M(1 : ceil(M/2)+1); 
% theta_T = theta_M; 
% delta_x_T = cos(theta_T) .* (r ./ (sqrt((cos(theta_T)).^2 + (a^2/b^2)*(sin(theta_T)).^2)) - 1); 
% delta_y_T = sin(theta_T) .* (r ./ (sqrt((b^2/a^2)*(cos(theta_T)).^2 + (sin(theta_T)).^2)) - 1); 
% figure; 
% hold all; 
% grid on; 
% box on; 
% plot(delta_x_T, '-k'); 
% plot(delta_y_T, '-r'); 
%  
% delta_a_T = r ./ (sqrt((cos(theta_T)).^2 + (a^2/b^2)*(sin(theta_T)).^2)); 
% delta_b_T = r ./ (sqrt((b^2/a^2)*(cos(theta_T)).^2 + (sin(theta_T)).^2)); 
% figure; 
% hold all; 
% grid on; 
% box on; 
% plot(delta_a_T, '-k'); 
% plot(delta_b_T, '-r'); 

%%Curvature calculation
% k=a*b/(sqrt(a^2/2+b^2/2)^3)
% k_r=(a+r)*(b+r)/(sqrt((a+r)^2/2+(b+r)^2/2)^3)
% t=0:0.01:2*pi
% figure;
% plot(t,a.*b./(sqrt(a.^2.*cos(t).^2+b.^2.*sin(t).^2).^3))
% hold on;
% ar=a+r;
% br=b+r;
% plot(t,ar.*br./(sqrt(ar.^2.*cos(t).^2+br.^2.*sin(t).^2).^3),'b')
% figure;
% plot(t,a^2.*cos(t).^2+b^2.*sin(t).^2)

function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
%     hold off
end



