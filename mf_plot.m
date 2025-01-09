clc, clear all, close all

x = -0.25:0.1:1.25;

x0 = -0.25:0.01:0;
x1 = 0:0.01:0.25;
x2 = 0.25:0.01:0.5;
x3 = 0.5:0.01:0.75;
x4 = 0.75:0.01:1;
x5 = 1:0.01:1.25;

y1a = 1+4.*x0;
y1b = 1-4.*x1;

y2a = 4.*x1;
y2b = 2*(1-2.*x2);

y3a = 4.*x2-1;
y3b = 3-4.*x3;

y4a = 2*(2.*x3-1);
y4b = 4*(1-x4);

y5a = 4.*x4-3;
y5b = 5-4.*x5;

figure('Name','Membership Functions')
plot(x0,y1a,x1,y1b,x1,y2a,x2,y2b,x2,y3a,x3,y3b,x3,y4a,x4,y4b,x4,y5a,x5,y5b)
title('Equi-partitioned universes of discourse (UD)')
xlabel('input: x')
ylabel('output: \mu')
axis([-0.25 1.25 0 1])

%% Liquid Side Condition (Gas Side is neglected)

% 1.a Pump Set

%% 11.4. System Condition Audit

%% 11.4.1. Filter
cf = 0.46; % rate of pressure drop for a given flow rate
cfp = 0.56; % observed in-service value
f = cf/cfp % filter condition

%% 11.4.2. Pump and Motor
ch = 0.5; % design pump characteristic
chp = 0.6; % observed value
p = (2/pi)*atan(ch/(chp-ch)) % pump condition

g = 0.8; % assumed pump-gland condition

%% 11.4.3. Pump and Filter Condition

F = epud(f) % pump condition
P = epud(p) % filter condition
G = epud(g) % gland condition

PG = rulePG(P,G)

%% 11.4.4. Pump Unit Condition

PU = rulePU(PG,F)

defPU = defuzzify(PU)

%% 11.4.5. By-Pass

l = 1.0; % observed valve gland seal condition

L = epud(l)

e = 0.85; % operational ease condition

E = epud_inv(e)

B = ruleB(E,L)

defB = defuzzify(B)

%% 11.4.6. Control Valve

d = 1.0; % the unit condition observation is "as good as new"

D = epud(d)

defD = defuzzify(D)

%% 11.4.7.Non-Return Valve

r = 0.95; % the observed condition

R = epud(r)

defR = defuzzify(R)

%% 11.4.8. Heat Exchanger Unit

% term = tan(pi*0.8585/2)
% t = (2/pi)*atan(st/(stp-st)) % pump condition
% ratio = 4.4247;
% stp_over_st = 1.226001710710897;
stp_over_st = 1.225506455158971;

t = (2/pi)*atan(1/(stp_over_st-1)) % pump condition

T = epud(t)

defT = defuzzify(T)

%% 11.4.9. System Water-Side Condition Metric
firingPU = find(PU > 0)
firingB = find(B > 0)
firingR = find(R > 0)
firingT = find(T > 0)

for i = 1:length(firingPU)
    outPU(i) = PU(firingPU(i));
end
for i = 1:length(firingB)
    outB(i) = B(firingB(i));
end
for i = 1:length(firingR)
    outR(i) = R(firingR(i));
end
for i = 1:length(firingT)
    outT(i) = T(firingT(i));
end

subscript1 = [4,4,4,4,4,4,4,5];

W2MF1 = ruleW(outPU,outB,outR,outT,subscript1);

W1 = zeros(1,5);
W1(4) = W2MF1(1);
W1(5) = W2MF1(2);

defW1 = defuzzify(W1)

subscript2 = [4,4,4,5,4,4,4,5];

W2MF2 = ruleW(outPU,outB,outR,outT,subscript2);

W2 = zeros(1,5);
W2(4) = W2MF2(1);
W2(5) = W2MF2(2);

defW2 = defuzzify(W2)

100*(defW2-defW1)/defW1

%% Functions

function mu = epud(x)
    mu = zeros(1,5);
    if x < 0
        mu(1) = 1+4*x;
    elseif x < 0.25
        mu(1) = 1-4*x;
        mu(2) = 4*x;
    elseif x < 0.5
        mu(2) = 2*(1-2*x);
        mu(3) = 4*x-1;
    elseif x < 0.75
        mu(3) = 3-4*x;
        mu(4) = 2*(2*x-1);
    elseif x < 1
        mu(4) = 4*(1-x);
        mu(5) = 4*x-3;
    else
        mu(5) = 5-4*x;
    end
end

function mu = epud_inv(x)
    mu = zeros(1,5);
    if x < 0
        mu(5) = 1+4*x;
    elseif x < 0.25
        mu(5) = 1-4*x;
        mu(4) = 4*x;
    elseif x < 0.5
        mu(4) = 2*(1-2*x);
        mu(3) = 4*x-1;
    elseif x < 0.75
        mu(3) = 3-4*x;
        mu(2) = 2*(2*x-1);
    elseif x < 1
        mu(2) = 4*(1-x);
        mu(1) = 4*x-3;
    else
        mu(1) = 5-4*x;
    end
end

function PG = rulePG(P,G)
    R11 = min(P(1),G(1));
    R21 = min(P(1),G(2));
    R12 = min(P(2),G(1));
    R22 = min(P(2),G(2));
    R32 = min(P(2),G(3));
    R23 = min(P(3),G(2));
    R33 = min(P(3),G(3));
    R43 = min(P(3),G(4));
    R34 = min(P(4),G(3));
    R44 = min(P(4),G(4));
    R54 = min(P(4),G(5));
    R45 = min(P(5),G(4));
    R55 = min(P(5),G(5));
    
    PG = [max([R11,R21,R12]),max([R22,R32,R23]),max([R33,R43,R34]),max([R44,R54,R45]),R55];    
end

function PU = rulePU(P,U)
    R11 = min(P(1),U(1));
    R21 = min(P(1),U(2));
    R12 = min(P(2),U(1));
    R22 = min(P(2),U(2));
    R32 = min(P(2),U(3));
    R23 = min(P(3),U(2));
    R33 = min(P(3),U(3));
    R43 = min(P(3),U(4));
    R34 = min(P(4),U(3));
    R44 = min(P(4),U(4));
    R54 = min(P(4),U(5));
    R45 = min(P(5),U(4));
    R55 = min(P(5),U(5));
    
    PU = [max([R11,R21]),max([R12,R22,R32,R23]),max([R33,R43,R34]),max([R44,R54,R45]),R55];    
end

function B = ruleB(E,L)
    R11 = min(E(5),L(1));
    R21 = min(E(5),L(2));
    R31 = min(E(5),L(3));
    R41 = min(E(5),L(4));
    R51 = min(E(5),L(5));
    R12 = min(E(4),L(1));
    R22 = min(E(4),L(2));
    R32 = min(E(4),L(3));
    R42 = min(E(4),L(4));
    R52 = min(E(4),L(5));
    R13 = min(E(3),L(1));
    R23 = min(E(3),L(2));
    R33 = min(E(3),L(3));
    R43 = min(E(3),L(4));
    R53 = min(E(3),L(5));
    R14 = min(E(2),L(1));
    R24 = min(E(2),L(2));
    R34 = min(E(2),L(3));
    R44 = min(E(2),L(4));
    R54 = min(E(2),L(5));
    R15 = min(E(1),L(1));
    R25 = min(E(1),L(2));
    R35 = min(E(1),L(3));
    R45 = min(E(1),L(4));
    R55 = min(E(1),L(5));
    
    B = [max([R11,R21,R31,R41,R51,R12,R13,R14,R15]),max([R22,R32,R42,R52,R23,R24,R25]),max([R33,R43,R53,R34,R35]),max([R44,R54,R45]),R55];    
end

function W = ruleW(outPU,outB,outR,outT,subscript)
    
    maxPU = max(outPU);
    
    R4 = zeros(1,8);
    R5 = zeros(1,8);
    
    if subscript(1) == 4
        R4(1) = min([maxPU,outB(1),outR(1),outT(1)]);
    elseif subscript(1) == 5
        R5(1) = min([maxPU,outB(1),outR(1),outT(1)]);
    end
    if subscript(2) == 4
        R4(2) = min([maxPU,outB(1),outR(1),outT(2)]);
    elseif subscript(2) == 5
        R5(2) = min([maxPU,outB(1),outR(1),outT(2)]);
    end
    if subscript(3) == 4
        R4(3) = min([maxPU,outB(1),outR(2),outT(1)]);
    elseif subscript(3) == 5
        R5(3) = min([maxPU,outB(1),outR(2),outT(1)]);
    end
    if subscript(4) == 4
        R4(4) = min([maxPU,outB(1),outR(2),outT(2)]);
    elseif subscript(4) == 5
        R5(4) = min([maxPU,outB(1),outR(2),outT(2)]);
    end
    if subscript(5) == 4
        R4(5) = min([maxPU,outB(2),outR(1),outT(1)]);
    elseif subscript(5) == 5
        R5(5) = min([maxPU,outB(2),outR(1),outT(1)]);
    end
    if subscript(6) == 4
        R4(6) = min([maxPU,outB(2),outR(1),outT(2)]);
    elseif subscript(6) == 5
        R5(6) = min([maxPU,outB(2),outR(1),outT(2)]);
    end
    if subscript(7) == 4
        R4(7) = min([maxPU,outB(2),outR(2),outT(1)]);
    elseif subscript(7) == 5
        R5(7) = min([maxPU,outB(2),outR(2),outT(1)]);
    end
    if subscript(8) == 4
        R4(8) = min([maxPU,outB(2),outR(2),outT(2)]);
    elseif subscript(8) == 5
        R5(8) = min([maxPU,outB(2),outR(2),outT(2)]);
    end
    
	W = [max(R4),max(R5)];    
end

function out = defuzzify(R)

    out = (R(1)*0 + R(2)*0.25 + R(3)*0.5 + R(4)*0.75 + R(5)*1)/(sum(R));

end
