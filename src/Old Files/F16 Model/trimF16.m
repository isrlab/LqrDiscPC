clear; clc;
d2r = pi/180;
%Trim for steady-level flight at various velocities
Vtrims = [300:100:900]; % ft/s

for i=1:numel(Vtrims)
    Vtrim = Vtrims(i);
    
    X0 = [500 0 0 0]'; % Initial guess
    U0 = [10000 0]'; % Initial guess
    
    %        V            a       th    q    T    de
    UB = [Vtrim    45*d2r   45*d2r   0  19000  25];
    LB = [Vtrim   -20*d2r  -20*d2r   0  1000  -25];
    P0 = [X0;U0]; % Initial guess for trim states and control
    Aeq = [0 -1 1 0 0 0]; % Constraints for gamma = 0
    Beq = 0;
    
    options = optimset('TolFun',1e-16,'TolX',1e-16,'Algorithm','sqp');
    % trim
    [x,fval,exitflag,output] = fmincon(@trimCost,P0,[],[],Aeq,Beq,LB,UB,[],options);
    xtrim(:,i) = x(1:4);
    utrim(:,i) = x(5:6);
    str1 = ['Objective value: ' num2str(fval)];
    str2 = ['xtrim = [' num2str(xtrim(:,i)') ']'];
    str3 = ['utrim = [' num2str(utrim(:,i)') ']'];
    
    disp(str1);
    disp(str2);
    disp(str3);
    optVal(i) = fval;
    
    % linearize
    [A,B,C,D] = linmod('f16longi_sim',xtrim(:,i),utrim(:,i));
    C = eye(4); D = zeros(4,2);
    f16ss(i).sys = ss(A,B,C,D);
    disp(i);
end
save trimmedF16 f16ss Vtrims xtrim utrim
% % Plot frsp
% bodemag(f16ss);
%%
load trimmedF16
scale = [1 1/d2r 1/d2r 1];
for i=1:size(xtrim,2)
    for j=1:size(xtrim,1)
        fprintf(1,'%.3f & ',scale(j)*xtrim(j,i));
    end
    for j=1:size(utrim,1)
        fprintf(1,'%.3f & ',utrim(j,i));
    end
    fprintf(1,'\n');
end
