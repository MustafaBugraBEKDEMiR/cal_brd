%Calculation of the GPS satellite position using the navigation message
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spos]=cal_brd(eph, brd)
%calculation toe
a = fix(brd(3,1)/86400)
toe=brd(3,1)-a*86400
%definitions of matrix
crs=brd(1,2)
deltan=brd(1,3)
mo=brd(1,4)
cuc=brd(2,1)
e=brd(2,2)
cus=brd(2,3)
sqrta=brd(2,4)
cic=brd(3,2)
omegao=brd(3,3)
cis=brd(3,4)
io=brd(4,1)
crc=brd(4,2)
omega=brd(4,3)
omegadot=brd(4,4)
idot=brd(5,1)
%formulas
tk=eph-toe
mk=mo+((sqrt(3.986005/power(sqrta*sqrta,3))+deltan))*tk
en  = mk;
	ens = en - (en-e*sin(en)- mk)/(1 - e*cos(en));
	while ( abs(ens-en) > eps )
		en = ens;
		ens = en - (en - e*sin(en) - mk)/(1 - e*cos(en));
	end;
	ek = ens
vk=atan(sqrt(1-power(e,2))*sin(ek)/cos(ek)-e)
uk=omega+vk+cuc*cos(2*omega+vk)+cus*sin(2*omega+vk)
rk=(sqrta*sqrta)*(1-e*cos(ek))+crc*cos(2*omega+vk)+crs*sin(2*omega+vk)
ik=io+idot*tk+cic*cos(2*omega+vk)+cis*sin(2*omega+vk)
hk=omegao+(omegadot-omega)*tk-omega*toe
xkk=rk*cos(uk)
ykk=rk*sin(uk)
xk=xkk*cos(hk)-ykk*cos(ik)*sin(hk)
yk=xkk*sin(hk)+ykk*cos(ik)*cos(hk)
zk=ykk*sin(ik)
[xk,yk,zk]