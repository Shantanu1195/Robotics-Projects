s=tf('s');
LA=1;tau=0.001;RA=1;Kp=400;Ki=5;Kd=20;
sys=((Kp*tau+Kd)*s^2+(Kp+Ki*tau)*s+Ki)/(LA*tau*s^3+(LA+RA*tau+Kp*tau+Kd)*s^2+(Kp+Ki*tau+RA)*s+Ki);
bode(sys)
