%Gamma transport
clc; clear;
global w_reg air_reg  det_reg void_reg1 void_reg2 boundary_reg

%tank width
 tank_width=20:2:26;


%water level
w_level=120;

%Energy cutoff
E_cutoff=0.0015; %Mev

flag=true;
%save detector response
det_res=zeros(1,length(tank_width));
for j=1:length(tank_width)



% tank boundary
%left, right, top, bottom
bl=0; br=tank_width(j); bt=200; bb=0;


%water region
%left right bottom top
w_reg=[bl br bb w_level];

%Air region
%left right bottom top
air_reg=[bl br w_level bt];

%Detector region
%left right bottom top
det_reg=[br br+3 90 110];

%void region
%left right bottom top
void_reg1=[br br+3 110 200];
void_reg2=[br br+3 0 90];

%outside boundary
%left right bottom top
boundary_reg=[bl br+3 bb bt];





%Particle track save
%X_dim=[]; Y_dim=[];
%----------------------------------------%
%particle initiate

%Number of source particels
Nps=20;
%number of particles detected
detected=0;
for i=1:Nps
    %Particle energy
E_p=0.662; %Mev
%source disersion angle 30 degrees
%initialization of source particle angle
%from -15 to 15 degrees
a_min=-15; a_max=15;
%random uniform angle 
S_theta=(a_min+rand*(a_max-a_min));


x=[]; y=[];
%initial position of the particle
x(1)=bl; y(1)=100;
s=1;
x(2)=x(1)+s*cosd(S_theta);
y(2)=y(1)+s*sind(S_theta);

alive=true;
it=2;
%-----------------------------------%
%           Particle Track          %
    while alive
        Sig_pe=0; Sig_c=0; Sig_s=0; Sig_a=0;
        [Sig_pe,Sig_c,Sig_s,Sig_a]=Xsec(x(it),y(it),E_p);
        
        Sig_t=Sig_pe+Sig_c+Sig_s+Sig_a;
        it=it+1;
        %interaction length
        s=-log(rand)/Sig_t;
        %update position data
        x(it)=x(it-1)+s*cosd(S_theta);
        y(it)=y(it-1)+s*sind(S_theta);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (x(it)<=boundary_reg(1))||(x(it)>=boundary_reg(2))||(y(it)<=boundary_reg(3))||(y(it)>=boundary_reg(4))
           
        alive=false;
        if x(it)<=boundary_reg(1)
            m=(y(it)-y(it-1))/(x(it)-x(it-1));
            x(it)=boundary_reg(1);
            y(it)=y(it-1)+m*(x(it)-x(it-1));
        elseif x(it)>=(boundary_reg(2))
            m=(y(it)-y(it-1))/(x(it)-x(it-1));
            x(it)=boundary_reg(2);
            y(it)=y(it-1)+m*(x(it)-x(it-1));
        elseif y(it)<=boundary_reg(3)
            m=(y(it)-y(it-1))/(x(it)-x(it-1));
            y(it)=boundary_reg(3);
            x(it)=(y(it)-y(it-1))/m;
        elseif y(it)>=boundary_reg(4)
            m=(y(it)-y(it-1))/(x(it)-x(it-1));
            y(it)=boundary_reg(4);
            x(it)=(y(it)-y(it-1))/m;  
        end
         
            
        end
        %detector
        if (x(it)>=det_reg(1))&& (y(it)>=det_reg(3))&&(y(it)<=det_reg(4))
                detected=detected+1; 
                break
        end
        if alive 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %process type
        r=rand;
        if r<Sig_pe/Sig_t
            %photoelectric
            
            alive=false;
        elseif (r>Sig_pe/Sig_t)&&(r<(Sig_pe+Sig_c)/Sig_t)
            %compton
            [angle,E_s]=Comp_scat(E_p);
            E_p=E_s;
            %Angle update after scattering
            S_theta=S_theta+angle;
            
             %Energy cutoff limit
             if E_p<E_cutoff
             alive=false;
             break
             end
             
        elseif (r>(Sig_pe+Sig_c)/Sig_t) && (r<(Sig_pe+Sig_c+Sig_s)/Sig_t)
            %in air or void
            alive=true;
            
        else
            
            
           alive=false;
           break
        end
        else 
            break
        end
    end
%Plot of 10 particle tracks
if flag
  if i>10
      flag=false;
  end
plot(x,y)
hold on
end


end
det_res(j)=detected;
end
%____________________Plots_________________________%


%Geometry plot


%Air region
rectangle('Position',[air_reg(1) air_reg(3) air_reg(2) (air_reg(4)-air_reg(3))])
text((air_reg(1)+air_reg(2))/2,(air_reg(3)+air_reg(4))/2,'Air Region')
%water region
rectangle('Position',[0 0 w_reg(2) w_reg(4)])
text((w_reg(1)+w_reg(2))/2,(w_reg(3)+w_reg(4))/2,'Water Region')
%boundary region
rectangle('Position',[0 0 boundary_reg(2) boundary_reg(4)])


%detector region
rectangle('Position',[det_reg(1) det_reg(3) 3 20])
text((det_reg(1)+det_reg(2))/2,(det_reg(3)+det_reg(4))/2,'Detector Region')

title('Monte Carlo Gamma-ray Transport Particle tracks')
xlabel('Dimension in X (cm)')
ylabel('Dimension in Y (cm)')


%detector response plot
figure(2)
plot(tank_width,det_res)
xlabel('Tank width (cm)'); ylabel('Detector Reponse')

%________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Functions                            %

%Cross section calculation at different energies for water

%Function returns Photoelectric cross section at different energies
%Tabulated PE cross section data for H2O
function x_secPE=PE_x_sec(E)
    load('PE_xsec.mat', 'PE_xsec')
    [ ~, ix ] = min( abs( PE_xsec(1,:)-E ) );
    if ix==1
    x_secPE=interp1(PE_xsec(1,ix:ix+1),PE_xsec(2,ix:ix+1),E);
    elseif ix==length(PE_xsec)
    x_secPE=interp1(PE_xsec(1,ix-1:ix),PE_xsec(2,ix-1:ix),E);
    else
    x_secPE=interp1(PE_xsec(1,ix-1:ix+1),PE_xsec(2,ix-1:ix+1),E);
    end

%x_secPE=interp1(PE_xsec(1,:),PE_xsec(2,:),E);

%cross section

end

%compton scattering cross section
function x_Comp=Comp_x_sec(E)
Z=10;   %for water
ao=E/0.511;
re=2.82E-13; %cm
K=pi*re^2/ao;
eta=1+0.222037*ao;
c1=1.651035;
c2=9.34022;
c3=-8.325004;
d1=12.501332;
d2=-14.200407;
d3=1.699075;
x_Comp=K*ao*(c1*eta^2+c2*eta+c3)/(eta^3+d1*eta^2+d2*eta+d3);
%number density of water
N=33.3679E21; %cm-3
%macroscopic cross section
x_Comp=x_Comp*Z*N;
% sig_c2=2*K*((1+ao)/ao*((2+2*ao)/(1+2*ao)-log(1+2*ao)/ao)+log(1+2*ao)/(2*ao)-(1+3*ao)/(1+2*ao)^2);
% disp(sig_c2)
end

%Function returns Compton cross section at different energies
%Compton scattering angle and Energy 
function [theta, En]=Comp_scat(E)
K=E/0.511; %E in MeV
%PDF for x=cos(theta)
f=@(x) (1+x^2+(K^2*(1-x^2))/(1+K*(1-x)))/(1+K*(1-x))^2;
%Envelop fucntion for rejection sampling
co=2*(2*K^2+2*K+1)/(2*K+1)^3;
b=(1+co/2)/(1-co/2);
a=2*(b-1);
h=@(x) a/(b-x);
%Sampling cos(theta) from h(x) using relation
% x=b-(b+1)*(co/2)^r where r is uniform random (0:1)
fxr=@(r) b-(b+1)*(co/2)^r;

%Rejection sampling

flag=1;
while flag
    u1=rand;
    u2=rand;
    if u2<f(fxr(u1))/h(fxr(u1))
        x1=fxr(u1);
        flag=0;
    else
        continue
    end
end
if rand>0.5
    theta=acosd(x1);
else
    theta=-acosd(x1);
    
end
En=E/(1+K*(1-cosd(theta)));
end


%interaction length function for region
function [Sig_pe,Sig_c,Sig_s,Sig_a]=Xsec(x,y,E_p)
global w_reg air_reg  det_reg void_reg1 void_reg2 boundary_reg

    water_region=((x>=w_reg(1))&&(x<=w_reg(2))&&(y>=w_reg(3))&&(y<=w_reg(4)));
    air_region=(x>=air_reg(1))&&(x<=air_reg(2))&&(y>=air_reg(3))&&(y<=air_reg(4));
    void_region=(x>=void_reg1(1))&&(x<=void_reg1(2))&&(y>=void_reg1(3))&&(y<=void_reg1(4))||(x>=void_reg2(1))&&(x<=void_reg2(2))&&(y>=void_reg2(3))&&(y<=void_reg2(4));
 %   det_region=(x>=det_reg(1))&&(x<=det_reg(2))&&(y>=det_reg(3))&&(y<=det_reg(4));
 %   outside_boundary=(x<=boundary_reg(1))&&(x>=boundary_reg(2))&&(y<=boundary_reg(3))&&(y>=boundary_reg(4));
    if water_region
       % disp('water region')
     %cross section
     Sig_pe=PE_x_sec(E_p); %photoelectric
     Sig_c=Comp_x_sec(E_p); %Compton
     Sig_s=0;
     Sig_a=0;
     
    
    elseif (air_region || void_region)
        Sig_pe=0;
        Sig_c=0;
        Sig_a=0;
        Sig_s=1;
       % disp('air void')
    else
        Sig_pe=0;
        Sig_c=0;
        Sig_a=1E3;
        Sig_s=0;
       % disp('boundary detector')
    end

end
