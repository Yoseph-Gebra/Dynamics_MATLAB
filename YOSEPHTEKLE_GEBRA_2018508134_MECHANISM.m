%CLEAR: Variables and command window in MATLAB

clc,clear
% INPUT: PHYSICAL PARAMETERS of MECHANISM (centemeter)
X=0.5493; a2=0.3; a3=0.6; a4=0.4394 ;a5=0.9558; Y=0.466;
% INPUT: Maximum Iteration Number Nmax
Nmax=100;
% INPUT: INITIAL GUESS VALUES for a6, th3, th4 and, th5 to respectively
x=[245*pi/180,305*pi/180,20*pi/180, 1.05];
% INPUT: ERROR TOLERANCE
xe=0.001*abs(x);
% INPUT: SYSTEM INPUTS (th2,w2,al2)
dth=5*pi/360;
th2=0*pi/180:dth:360*pi/180;

w2=8*ones(1,length(th2));
al2=0*ones(1,length(th2));
%----------------------------------------------
xe=transpose(abs(xe));
kerr=1; %If kerr=1, results are not converged
%% 
for k=1:1:length(th2);
 for n=1:Nmax
 %----------------------------------------------
 %Assign initial guess to unknowns

 th3(k)=x(1); th4(k)=x(2); 
 th5(k)=x(3);a6(k)=x(4);
 
 % INPUT: JACOBIAN Matrix
 J=zeros(4,4);
 J(1,1)=-a3*sin(th3(k));J(1,2)=-a4*sin(th4(k));
 J(2,1)=a3*cos(th3(k));J(2,2)=a4*cos(th4(k));
 J(3,4)=-1;J(3,2)=-a4*sin(th4(k));J(3,3)=-a5*sin(th5(k));
 J(4,2)=a4*cos(th4(k));J(4,3)=a5*cos(th5(k));
 % INPUT: Function f
 f=zeros(4,1);
 f(1,1)= -X - a2*cos(th2(k)) - a3*cos(th3(k)) - a4*cos(th4(k));
 f(2,1)= -Y - a2*sin(th2(k)) - a3*sin(th3(k)) - a4*sin(th4(k));
 f(3,1)= a6(k) - a5*cos(th5(k)) - a4*cos(th4(k));
 f(4,1)= -a5*sin(th5(k)) - a4*sin(th4(k));
 %----------------------------------------------
 eps=inv(J)*f;x=x+transpose(eps);
% if abs(eps)<xe
%  kerr=0;break
%  end
 end
%  if kerr==1
%  print('error')
%  end
 
 th3(k)=x(2); a6(k)=x(1);
 th4(k)=x(3);th5(k)=x(4);
 
 %---velocity---------------------------
 fv(1,1)=-w2(k)*a2.*sin(th2(k));
 fv(2,1)=-w2(k)*a2.*cos(th2(k));
 fv(3,1)=0;
 fv(4,1)=0;
 vel=inv(J)*fv;
 w3(k)=vel(1);w4(k)=vel(2);
 w5(k)=vel(3);V6(k)=vel(4);
%---acceleration---------------------------
 fa(1,1)=-al2(k)*a2*sin(th2(k))-w2(k)^2*a2*cos(th2(k)) - w3(k)^2*a3*cos(th3(k)) - w4(k)^2 *a4*cos(th4(k));
 fa(2,1)= -al2(k)*a2*cos(th2(k)) + w2(k)^2*a2*sin(th2(k)) +w3(k)^2*a3*sin(th3(k)) + w4(k)^2*a4*sin(th4(k));
 fa(3,1)=-w4(k)^2*a4*cos(th4(k))- a5*w5(k)^2*cos(th5(k));
 fa(4,1)= w4(k)^2*a4*sin(th4(k))+ a5*w5(k)^2*sin(th5(k));
 acc=inv(J)*fa;
 al3(k)=acc(1);al4(k)=acc(2);
 al5(k)=acc(3);A6(k)=acc(4);
end
% Angle: radian --> degree
th2d=th2*180/pi;
th3d=th3*180/pi;
th4d=th4*180/pi;
th5d=th5*180/pi;

%--------GRAPHINGS---------------
figure(1),

subplot(4,3,1),plot(th2d,th3d,'g','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\theta_3(^o)'), title('\theta_2(^o) and \theta_3(^o) '), grid on;
subplot(4,3,2),plot(th2d,w3,'k','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\omega_3(r/s)'),title('\theta_2(^o) and \omega_3 '),grid on;
subplot(4,3,3),plot(th2d,al3,'b','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\alpha_3(r/s^2)'),title('\theta_2(^o) and \alpha_3 '),grid on;


subplot(4,3,4),plot(th2d,th4,'g','linewidth',2),xlabel('\theta_2(^o)'),ylabel('\theta_4 (^o)'),title('\theta_2(^o) and \theta_4 '),grid on;
subplot(4,3,5),plot(th2d,w4,'k','linewidth',2),xlabel('\theta_2(^o)'),ylabel('\omega_4 (r/s)'),title('\theta_2(^o) and \omega_4 '),grid on;
subplot(4,3,6),plot(th2d,al4,'b','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\alpha_4(r/s^2)'),title('\theta_2(^o) and \alpha_4 '),grid on;


subplot(4,3,7),plot(th2d,th5,'g','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\theta_5(^o)'),title('\theta_2(^o) and \theta_5 '),grid on;
subplot(4,3,8),plot(th2d,w5,'k','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\omega_5(r/s)'),title('\theta_2(^o) and \omega_5 '),grid on;
subplot(4,3,9),plot(th2d,al5,'b','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('\alpha_5(r/s^2)'),title('\theta_2(^o) and \alpha_5 '),grid on;


subplot(4,3,10),plot(th2d,a6,'g','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('a6(m)'),title('\theta_2(^o) and a6 '),grid on;
subplot(4,3,11),plot(th2d,V6,'k','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('V6(m/s)'),title('\theta_2(^o) and V6 '),grid on;
subplot(4,3,12),plot(th2d,A6,'b','linewidth',2),xlabel('\theta_2 (^o)'),ylabel('A6(m/s^2)'),title('\theta_2(^o) and A6 '),grid on;


sgtitle('YOSEPH TEKLE GEBRA ID-NO....20018508134')
