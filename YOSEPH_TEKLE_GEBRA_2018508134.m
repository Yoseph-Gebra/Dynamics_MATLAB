clc,clear
% GIVEN VALUES
W=[0 -44497.37];
a2=7.25; a31=1.43; a32=2.73; a4=5.79;
ag3=(a31+a32)/2;
a3=a31+a32;
w2=8;
m2=500;m3=350; m4=400;
I3=504.47; I4=1117.47;
G2=[0 -4905]; 
G3=[0 -3433.5];
G4=[0 -3924];
err=[0.001;0.001]; % ERROR VOR THE NETON RAPHSON CALCULATION OF 
% INITIAL GUESS FOR THETA
x=[290*pi/180,237.34*pi/180,127.97*pi/180];
 
%% POSITION ANALYSIS
pr=[-(a4*cos(x(3))+a31*cos(x(2))+a2*cos(x(1)))
    -(a4*sin(x(3))+a31*sin(x(2))+a2*sin(x(1)))];
X_D=pr(1,1);
Y_D=pr(2,1);
%% FORCE ANALYSIS
Nmax=20; %number of iteration
dth=1*pi/360;
for th2=280*pi/180:dth:306*pi/180
    
     for n=1:200
X=[127.97*pi/180; 237.34*pi/180]; % initial guess for theta4 and theta3 respectively
% The jacobian for the neton raphson method 
j(1,1)=-a4*sin(X(1)); j(1,2)=-a31*sin(X(2));
j(2,1)=a4*cos(X(1)); j(2,2)=a31*cos(X(2));
b=[-(X_D + a2*cos(th2)+a4*cos(X(1))+a31*cos(X(2)))
    -(Y_D +a2*sin(th2)+a4*sin(X(1))+a31*sin(X(2)))];
eps=inv(j)*b;
X=X+(eps);
if abs(eps)<=err
    break
elseif n==Nmax
end
     end
n=n+1;
% Final result for the two angles 
th4=X(1,1);
th3=X(2,1);
%% velocity analysis
vx=[-a4*sin(X(1)) -a31*sin(X(2))
    a4*cos(X(1)) a31*cos(X(2))];
vr=[w2*a2*sin(th2)
    -w2*a2*cos(th2)];
VR=inv(vx)*vr;  
w4=VR(1,1);
w3=VR(2,1)
%% acceleration analysis
ax=[-a4*sin(X(1)) -a31*sin(X(2))
    a4*cos(X(1)) a31*cos(X(2))];
ar=[w2^2*a2*cos(th2)+w3^2*a31*cos(X(2))+w4^2*a4*cos(X(1))
    w2^2*a2*sin(th2)+ w3^2*a31*sin(X(2))+ w4^2*a4*sin(X(1))];
AR=inv(ax)*ar;  
alp4=AR(1,1);
alp3=AR(2,1);
%% Trial for Force equations
REo2=[1.45*cos(th2) 1.45*sin(th2)]; RGo2=[(a2/2)*cos(th2) (a2/2)*sin(th2)];RBo2=[a2*cos(th2) a2*sin(th2)];
RGo4=[(a4/2)*cos(th4) (a4/2)*sin(th4)]; RAo4=[a4*cos(th4) a4*sin(th4)];
RAG3=[ag3*cos(th3) ag3*sin(th3)]; RCA=[a3*cos(th3) a3*sin(th3)];
RBA=[a31*cos(th3) a31*sin(th3)];
VA=[-w2*RAo4(1,2) w2*RAo4(1,1)];
M_VA=((VA(1,1))^2 + (VA(1,2))^2)/sqrt((RAo4(1,2)^2+RAo4(1,1)^2));
VG2=[-w2*RGo2(1,2) w2*RGo2(1,1)];
VE=[-w2*REo2(1,2) w2*REo2(1,1)];
M_VG2=((VG2(1,1))^2 + (VG2(1,2))^2 )/sqrt((RGo2(1,2)^2+RGo2(1,1)^2));
aG2n=M_VG2*[cos(th2+pi) sin(th2+pi)];
aG2=aG2n;
VG4=[-w4*RGo4(1,2) w4*RGo4(1,1)];
M_VG4=((VG4(1,1))^2 + (VG4(1,2))^2 )/sqrt((RGo4(1,2)^2+RGo4(1,1)^2));
aG4n=M_VG4*[cos(th4+pi) sin(th4+pi)];
aG4t=[-alp4*RGo4(1,2) alp4*RGo4(1,1)];
aG4=aG4n+aG4t;
aAn=M_VA*[cos(th4+pi) sin(th4+pi)];
aAt=[-alp4*RAo4(1,2) alp4*RAo4(1,1)];
aA=aAn+aAt;
VG3relA=[-w3*RAG3(1,2) w3*RAG3(1,1)];
VG3=[VG3relA(1,1)+VA(1,1) VG3relA(1,2)+VA(1,2)];
M_VG3=((VG3(1,1))^2 + (VG3(1,2))^2 )/sqrt(RAG3(1,2)^2+RAG3(1,1)^2);
aG3relAn=M_VG3*[cos(th3) sin(th3)];
aG3relAt=[-alp4*RAG3(1,2) alp4*RAG3(1,1)];
aG3=aG3relAn+aG3relAt+aA;
VCrelA=[-w3*RCA(1,2) w3*RCA(1,1)];
VC=VA+VCrelA;
Ti3=-I3*alp3;
FG3i=-(m3*aG3);
Ti4=-I4*alp4;
FG4i=-(m4*aG4);
FG2i=-(m2*aG2);

fx=[RBA(1,1) RBA(1,2) 0 0 0 0 0 0 0
    -1 0 1 0 0 0 0 0 0
    0 -1 0 1 0 0 0 0 0
    RBo2(1,1) RBo2(1,2) 0 0 REo2(1,2) 0 0 0 0
    1 0 0 0 0 1 0 0 0 
    0 -1 0 0 1 0 -1 0 0
    0 0 RAo4(1,1) RAo4(1,2) 0 0 0 0 0
    0 0 1 0 0 0 0 1 0
    0 0 0 1 0 0 0 0 1];


fr=[-(W(1,2)*RCA(1,1)+FG3i(1,1)*RAG3(1,2)+ FG3i(1,2)*RAG3(1,1)+Ti3+G3(1,2)*RAG3(1,1))
    -(W(1,2)+FG3i(1,2)+G3(1,2))
     -FG3i(1,1)
     -(G4(1,2)*RAo4(1,1)+FG4i(1,1)*RGo4(1,2)+ FG4i(1,2)*RGo4(1,1)+ Ti4)
     -(G4(1,2)+FG4i(1,2))
     -FG4i(1,1)
     -(G2(1,2)*RGo2(1,1)+FG2i(1,1)*RGo2(1,2)+ FG2i(1,2)*RGo2(1,1))
     -(G2(1,2)+FG2i(1,2))
     -FG2i(1,1)];


FR=inv(fx)*fr;
F23Y=FR(1,1);
F23X=FR(2,1);
F43X=FR(3,1);
F43Y=FR(4,1);
F=FR(5,1);
F12X=FR(6,1);
F12Y=FR(7,1);
F14Y=FR(8,1);
F14X=FR(9,1);
F14=sqrt(F14X^2+F14Y^2);
F23=sqrt(F23X^2+F23Y^2);
F12=sqrt(F12X^2+F12Y^2);
F43=sqrt(F43X^2+F43Y^2);
f=FR(5,1);

figure(1),

subplot(2,3,1),plot(th2,F43,'*'),xlabel('\theta_2 '),ylabel('F43'), title('\theta_2 and F43'), grid on,hold on;
subplot(2,3,2),plot(th2,F14,'*'),xlabel('\theta_2 '),ylabel('F14'), title('\theta_2 and F14'),grid on,hold on;
subplot(2,3,3),plot(th2,F23,'*'),xlabel('\theta_2 '),ylabel('F23'), title('\theta_2 and F23'),grid on,hold on;
subplot(2,3,4),plot(th2,F12,'*'),xlabel('\theta_2 '),ylabel('F12'), title('\theta_2 and F12'),grid on,hold on;
subplot(2,3,[5,6]),plot(th2,f,'*'),xlabel('\theta_2 '),ylabel('F'), title('\theta_2 and F'),grid on,hold on;

end 

