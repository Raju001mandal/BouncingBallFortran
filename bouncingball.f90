program bouncingball
implicit none

integer::i,n,n1,n2,n3,n4,n5
real::h,d,thetarad,theta,ux0,uy0,u,r,x,y,x0,y0,vxhit1,vyhit1,thit1,t,dt,vxhit2,vyhit2,tf,ux20,uy20,thit2,thetarad2,theta2 ,&
      xs1,ys1,xs2,ys2,vxhit3,vyhit3,ux30,uy30,thit3,thetarad3,theta3,ux40,uy40,xs3,ys3,vxhit4,vyhit4,thit4,thetarad4,theta4 ,&
      ux50,uy50,xs4,ys4,vxhit5,vyhit5,thit5,thetarad5,theta5,xs5,ys5
      
real,parameter::pi=acos(-1.0),t0=0.0,g=9.81,cor=0.8

!real,dimension(1:1000)::x,y

print*,"give u,h,d"
read*,u,h,d

print*,"give the step size dt of time"
read*,dt




r=h/d
thetarad=atan(r)
theta = 180*thetarad/pi
!=========================calculation of time of first hit===============================

ux0=u*cos(thetarad)
uy0=u*sin(thetarad)

vxhit1=ux0
vyhit1=sqrt(uy0**2+2*g*h)
thit1=abs((vyhit1-uy0)/g)
n1=int((thit1-t0)/dt)

!==============time of second hit and initialization of 2nd projection==================================================

ux20=cor*vxhit1
uy20=abs(cor*vyhit1)
thetarad2=abs(uy20/ux20)
theta2=180*thetarad2/pi
thit2=abs(2*uy20*sin(thetarad2)/g)
thit2=thit1+thit2+0.16
n2=int((thit2-t0)/dt)
vxhit2=ux20
vyhit2=-uy20

!==============time of third hit===========================================================================================================

ux30=cor*vxhit2
uy30=abs(cor*vyhit2)
thetarad2=abs(uy30/ux30)
theta3=180*thetarad3/pi
thit3=abs(2*uy20*sin(thetarad3)/g)
thit3=thit2+thit3+0.58
n3=int((thit3-t0)/dt)
vxhit3=ux30
vyhit3=-uy30
!==============================fourth hit================================================================================================

ux40=cor*vxhit3
uy40=abs(cor*vyhit3)
thetarad4=abs(uy40/ux40)
theta4=180*thetarad4/pi
thit4=abs(2*uy20*sin(thetarad4)/g)
thit4=thit4+thit3-0.12
n4=int((thit4-t0)/dt)
vxhit4=ux40
vyhit4=-uy40

!================================5th hit============================================================================================

ux50=cor*vxhit4
uy50=abs(cor*vyhit4)
thetarad5=abs(uy50/ux50)
theta5=180*thetarad5/pi
thit5=abs(2*uy20*sin(thetarad5)/g)
thit5=thit4+thit5-0.21
n5=int((thit5-t0)/dt)
vxhit5=ux50
vyhit5=-uy50

!=================================================================================================================================

x0=0
y0=h

print*,"n1= ",n1,"n2= ",n2,"n3= ",n3,"n4= ",n4,"n5= ",n5

open(1,file="newball.dat")
t=t0
do i=0,(n5+60)

t=i*dt+t0
if(i<=n1)then
x=ux0*t+x0
y=-uy0*t-0.5*g*t**2+y0
xs1=x
ys1=y
else if(i>n1 .and. i<=n2)then
t=(i-n1)*dt+t0
x=xs1+ux20*t
y=ys1+uy20*t-0.5*g*t**2
xs2=x
ys2=y
else if(i>n2 .and. i<=n3)then
t=(i-n2)*dt+t0
x=xs2+ux30*t
y=ys2+uy30*t-0.5*g*t**2
xs3=x
ys3=y
else if(i>n3 .and. i<=n4)then
t=(i-n3)*dt+t0
x=xs3+ux40*t
y=ys3+uy40*t-0.5*g*t**2
xs4=x
ys4=y
else
t=(i-n4)*dt+t0
x=xs4+ux50*t
y=ys4+uy50*t-0.5*g*t**2
xs5=x
ys5=y
end if
write(1,*)x,y,t

end do

end program
