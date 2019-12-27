clc;
clear;
t=linspace(0,20,500);
a = pi/2; 
b = 0; 
theta0=[a b 0 0];I1=1;I2=4;m1=1;m2=2;L1=1;L2=2;Lc1=0.5;Lc2=1;g=9.8;damp=1.5;
[theta1, theta2,dtheta1,dtheta2]= ode45_ex(theta0,I1,I2,m1,m2,L1,L2,Lc1,Lc2,g,damp);
x1=-L1*cos(theta1);
y1=L1*sin(theta1);
x2=-L1*cos(theta1)-L2*cos(theta1+theta2);
y2=L1*sin(theta1)+L2*sin(theta1+theta2);

figure(1);
hold on;
plot(t,theta1,'r');
plot(t,theta2,'b');
legend('theta1','theta2');
 xlabel('t(s)');ylabel('theta(rad)');
 title('theta response');
 hold off;
filename=strcat('1theta response',num2str(a) ,'-',num2str(b) ,'.jpg');
F = getframe(gcf);
imwrite(F.cdata ,filename);
 
figure(2)
hold on;
subplot(2,2,1),plot(theta1,theta2, 'b');
xlabel('theta1(rad)'),ylabel('theta2(rad)'),title('phase potrait1');
subplot(2,2,2),plot(theta1,dtheta1, 'b');
xlabel('theta1(rad)'),ylabel('dtheta1(rad/s)'),title('phase potrait2');
subplot(2,2,3),plot(theta2,dtheta2, 'b');
xlabel('theta2(rad)'),ylabel('dtheta2(rad/s)'),title('phase potrait3');
subplot(2,2,4),plot(dtheta1,dtheta2, 'b');
xlabel('dtheta1(rad/s)'),ylabel('dtheta2(rad/s)'),title('phase potrait4');
hold off;
filename=strcat('1phase potrait',num2str(a) ,'-',num2str(b) ,'.jpg');
F = getframe(gcf);
imwrite(F.cdata ,filename);

figure(3)
kineticv=0.5*(m1*Lc1^2*dtheta1.^2+m2*(L1^2*dtheta1.^2+Lc2^2*(dtheta1+dtheta2).^2+2*Lc2*L1*dtheta1.*(dtheta1+dtheta2)).*cos(theta2));
kinetict=0.5*(I1*dtheta1.^2+I1*(dtheta1+dtheta2).^2)+kineticv;
position=-g*(m1*Lc1*cos(theta1)+m2*(L1*cos(theta1)+Lc2*cos(theta1+theta2)));
energy=position+kinetict;
hold on
plot(t,position);
plot(t,kinetict);
plot(t,energy);
xlabel('t(s)');ylabel('energy(J)');
title('energy diagram');
legend('U','K','U+K');
hold off 
filename=strcat('1energy diagram',num2str(a) ,'-',num2str(b) ,'.jpg');
F = getframe(gcf);
imwrite(F.cdata ,filename);

figure(4)
   Ncount=0;
   fram=0;
myVideo = VideoWriter('12-90-1'); % open video file 
myVideo.FrameRate = 24*4; % can adjust this
open(myVideo);
%become gif    
     for i=1:length(theta1)
         Ncount=Ncount+1;
         fram=fram+1;
         subplot(1,2,1);
         plot(0, 0,'.','markersize',20);
         hold on
         plot(y1(i),x1(i),'.','markersize',20);
         plot(y2(i),x2(i),'.','markersize',20);
         hold off
         line([0 y1(i)], [0 x1(i)],'Linewidth',2);
         axis([-(L1+L2) L1+L2 -(L1+L2) L1+L2]);
         line([y1(i) y2(i)], [x1(i) x2(i)],'linewidth',2);
         h=gca; 
         get(h,'fontSize') ;
         set(h,'fontSize',12);
         xlabel('Y(unit lenght)','fontSize',12);
         ylabel('X(unit lenght)','fontSize',12);
         title('dual pendulum','fontsize',14)
          frame=getframe(gcf);
     writeVideo(myVideo, frame);
        subplot(1,2,2);
         plot(t,theta1,'r', t,theta2,'b');
        hold on;
        plot(t(i),theta1(i),'ro');  
         plot(t(i),theta2(i),'bo');
        axis([0 20 -3 3]);
        hold off;
         xlabel('t(s)');ylabel('rad');
         fh = figure(4);
         set(fh, 'color', 'white'); 
         F=getframe;
         frame=getframe(gcf);
     writeVideo(myVideo, frame);
     end
close(myVideo)
     
movie(F,fram,20);


function [theta1, theta2,dtheta1,dtheta2] =ode45_ex(theta0,a0,b0,c0,d0,e0,f0,g0,h0,i0,damp)
I1=a0;I2=b0;m1=c0;m2=d0;L1=e0;L2=f0;Lc1=g0;Lc2=h0;g=i0;
t=linspace(0,20,500);
[t,theta] = ode45(@smd,t,theta0);
theta1=theta(:,1);
theta2=theta(:,2);

dtheta1=theta(:,3);
dtheta2=theta(:,4);


function dtheta = smd(t,theta)
dtheta(1,:) = theta(3);
dtheta(2,:) = theta(4);

a1=I1+I2+m1*Lc1^2+m2*(L1^2+Lc2^2+2*L1*Lc2*cos(theta(2)));
b1=I2+m2*(Lc2^2+L1*Lc2*cos(theta(2)));
c1=I2+m2*(Lc2^2+L1*Lc2*cos(theta(2)));
d1=I2+m2*Lc2^2;
e1=m1*g*Lc1*sin(theta(1))+m2*g*(L1*sin(theta(1))+Lc2*sin(theta(1)+theta(2)))-m2*L1*Lc2*sin(theta(2))*theta(3)*(2*theta(3)+theta(4))+damp*theta(3);
f1=m2*g*Lc2*sin(theta(1)+theta(2))+m2*L1*Lc2*sin(theta(2))*theta(3)^2+damp*theta(4);

dtheta(3,:)= (-e1*d1+b1*f1)/(a1*d1-c1*b1) ;
dtheta(4,:)= (-a1*f1+c1*e1)/(a1*d1-c1*b1) ;
end
end
