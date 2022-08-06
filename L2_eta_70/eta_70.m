%function script_eta_PM_scenario

format short

%addpath(genpath('/home/tvertesi/Documents/MATLAB/fq444/qubit4matlab'))

paulixyz 

% eta_planar_model = 2/pi = 0.6366

load W70i.txt

Wint = W70i;

[mA,mB] = size(Wint)


% generate EQ %%%%%%%%%%%%%

format short 


load grassc_3_1_70.txt

v = reshape(grassc_3_1_70,70,3);

for i=1:mA,
v(i,:) = v(i,:)/sqrt(v(i,:)*v(i,:)');
end;

size(v)

E = zeros(mA,mB);
for i1=1:mA,
    for i2=1:mB,
    va = v(i1,:);
    vb = v(i2,:);
    E(i1,i2) = va*vb';
    end;
end;

format short
% target point

% The value computed by L2 code.
L2 = 412667
% L2 = 412667

Q = trace(E*Wint')
% Q = 5.3672e+05

S = sum(sum(Wint))
% S = 194369

eta_crit = (L2-S)/(Q-S)
% eta_crit = 0.6376 


%dlmwrite('W70i.txt',Belli,'delimiter','\t','precision',6)
