function [SVpos,SVbias] = findSVpos(eph,time,DF)
% Calculate satellite positions in WGS-84 ECEF frame and satellite clock biases
% Inputs:
%   eph = GPS ephemeris
%   time = time to calculate SV positions [GPSweek GPSsec]
%   DF = dual frequency option (0=L1 only; 1=L1/L2 DF)
% Outputs:
%   SVpos = satellite positions in WGS-84 frame (m)
%   SVbias = satellite clock biases wrt GPS time (s)

% CT:
if isempty(eph) || isempty(time)
    SVpos = [];
    SVbias =[];
    return;
end


% Constants
c=utilities.c;             %speed of light (m/s)
omega_e=7.2921151467e-5;    %rotation rate of earth (rad/s)
mu=utilities.mu;             %Earth gravitational parameter (m^3/s^2)

nsat=length(eph(:,1));      %number of satellites
M0=eph(:,11);               %mean anomaly at reference time (rad)
dn=eph(:,10);               %mean motion difference (rad/s)
e=eph(:,12);                %eccentricity        
A=eph(:,9);                 %semimajor axis (m)
Omega0=eph(:,22);           %right ascension / longitude of ascending node (rad)
i0=eph(:,20);               %inclination angle (rad)
w=eph(:,13);                %argument of perigee (rad)
Omega_dot=eph(:,23);        %rate of right ascension (rad/s)
i_dot=eph(:,21);            %rate of inclination angle (rad/s)
cuc=eph(:,14);              %argument of latitude cos (rad)
cus=eph(:,15);              %argument of latitude sin (rad)
crc=eph(:,16);              %orbit radius cos (m)
crs=eph(:,17);              %orbit radius sin (m)
cic=eph(:,18);              %inclination cos (rad)
cis=eph(:,19);              %inclination sin (rad)
toe=eph(:,8);               %reference time of ephemeris (s)
IODE=eph(:,4);              %issue of data ephemeris
GPSweek=eph(:,6);           %GPS week number
toC=eph(:,25);              %SV clock correction (s)
tGD=eph(:,26);              %ionospheric group delay (s)
af0=eph(:,27);              %clock aging parameter (s)
af1=eph(:,28);              %clock aging parameter (s/s)
af2=eph(:,29);              %clock aging parameter (s/s/s)
URA=eph(:,32);              %user range accuracy variance (m^2)

% Time from ephemeris reference epoch
tk=time(:,2)-toe;
for k=1:nsat
    if tk(k)>302400
        tk(k)=tk(k)-604800;
    elseif tk(k)<-302400
        tk(k)=tk(k)+604800;
    end
end

n0=sqrt(mu./A.^3);
% Corrected mean motion
N=n0+dn;
% Mean anomaly
Mk=M0+N.*tk;
Eo=Mk; %old eccentric anomaly for use in iteration
error=1;
% Iterate to solve for eccentric anomaly
while abs(error)>1e-8
    E=Eo+(Mk-Eo+e.*sin(Eo))./(1-e.*cos(Eo));
    error=E-Eo;
    Eo=E;
end

% True anomaly
nu=atan2(sqrt(1-e.^2).*sin(E),cos(E)-e); 
% Argument of latitude
phi=nu+w;

% Second harmonic perturbations
du=cus.*sin(2*phi)+cuc.*cos(2*phi); %argument of latitude correction
dr=crs.*sin(2*phi)+crc.*cos(2*phi); %radius correction
di=cis.*sin(2*phi)+cic.*cos(2*phi); %inclination correction

% Corrected argument of latitude
muk=phi+du;
% Corrected radius
rk=A.*(1-e.*cos(E))+dr;
% Corrected inclination
ik=i0+di+i_dot.*tk;

% SV positions in orbital plane
xkp=rk.*cos(muk);
ykp=rk.*sin(muk);

% Corrected longitude of ascending node
Omega_k=Omega0+(Omega_dot-omega_e).*tk-omega_e.*toe;

% SV positions in ECEF frame
x1=xkp.*cos(Omega_k)-ykp.*cos(ik).*sin(Omega_k);
y1=xkp.*sin(Omega_k)+ykp.*cos(ik).*cos(Omega_k);
z=ykp.*sin(ik);
SVpos=[eph(:,1) x1 y1 z];

% Compute satellite clock biases
tc=time(:,2)-toC;
for k=1:nsat
    if tc(k)>302400
        tc(k)=tc(k)-604800;
    elseif tc(k)<-302400
        tc(k)=tc(k)+604800;
    end
end
dt=af0+af1.*tc+af2.*tc.^2-(2.*e.*sqrt(mu.*A).*sin(E))./(c^2);
if DF==0
    dt=dt-tGD;
end
SVbias=[eph(:,1) dt];

end