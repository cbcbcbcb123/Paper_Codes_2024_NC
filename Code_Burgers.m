clear
clc

KK=1;
Ecm_Array_Array=zeros(500,1);
Contraction_Arry_Arry=zeros(500,1);
Vactin_Arry_Arry=zeros(500,1);
N_link_0_Arry_Arry=zeros(500,1);
N_link_1_Arry_Arry=zeros(500,1);
Averageforce_Arry_Arry=zeros(500,1);

for kk=0.9:0.05:0.9 % kk 表示：Klong/Kshort
    kk
    % Parameters
    Klink=1;
    
    k2=10^0; % 初始模量:k2 稳态模量:k1k2/(k1+k2)
    k1=kk*k2/(1-kk);
    
    R1=10^-1;
    R2=10^1;
    
    A=R1/k1+R2/k1+R2/k2;
    B=R1*R2/k1/k2;
    C=R2;
    D=R1*R2/k1;
    
    t=0;
    V0=120;
    Num_links=50;
    Nmyosin=50;
    Fstall=2;
    
    ECMloc=0;
    N_link_1=0;
    Contraction=0;
    ECMloc_old2=0;
    ECMloc_old1=0;
    sumxi_old2=0;
    sumxi_old1=0;
    nc_old2=0;
    nc_old1=0;
    
    Fcr=3;
    dint=300;
    dadd=4;
    dmax=3000;
    dmin=300;
    Punfold=1;
    kont=0.001;
    koff0=0.1;
    kF=2;
    
    % N_links states:
    % 1,[(0 no bound) (1 bound)]
    % 2,[Single clutch length]
    % 3,[the rate of state1 -> state0].
    
    N_links=zeros(Num_links,3);
    
    Ncyc=1000;
    TimeArry=zeros(Ncyc,1);
    N_link_0_Arry=zeros(Ncyc,1);
    N_link_1_Arry=zeros(Ncyc,1);
    Contraction_Arry=zeros(Ncyc,1);
    Vactin_Arry=zeros(Ncyc,1);
    Averageforce_Arry=zeros(Ncyc,1);
    VECAD_Arry=zeros(Ncyc,1);

    for i=1:Ncyc
        %======
        % Caculate on-rate
        Averageforce=max(Contraction/N_link_1,0);

        VECAD = 1 / (Averageforce + 1)

        dint = 10 * R2 / (R2 + 10);

        dint=min(dint,dmax);
        kon=kont*dint;
        %======
        % Caculate Force-dependent off-rate
        for j=1:Num_links
            if N_links(j,1)==1
                N_links(j,3)=koff0*exp((N_links(j,2)-ECMloc)*Klink/kF);
            end
        end
        %======
        % Caculate time
        Tlink=100*ones(Num_links,1);
        for j=1:Num_links
            if N_links(j,1)==0
                Tlink(j,1)=-log(rand)/kon;
            elseif N_links(j,1)==1
                Tlink(j,1)=-log(rand)/N_links(j,3);
            end
        end
        %======
        % Min time
        [tmin,X_min_loc]=min(Tlink(:));
        %======
        % Go on
        t=t+tmin;
        %======
        % Things happen
        if N_links(X_min_loc,1)==0
            N_links(X_min_loc,1)=1;
        elseif N_links(X_min_loc,1)==1
            N_links(X_min_loc,1)=0;
        end
        %======
        % Update clutch strectch
        Vactin=V0*(1-Contraction/(Nmyosin*Fstall));
        Vactin=max(Vactin,0);
        for j=1:Num_links
            if N_links(j,1)==1
                N_links(j,2)=N_links(j,2)+Vactin*tmin;
            end
        end
        %======
        % Sum strectch and bind number
        sumXi=0;
        N_link_0=0;
        N_link_1=0;
        for j=1:Num_links
            if N_links(j,1)==0
                N_link_0=N_link_0+1;
            end
            if N_links(j,1)==1
                N_link_1=N_link_1+1;
                sumXi=sumXi+N_links(j,2);
            end
        end
        
        sumxi=sumXi;
        nc=N_link_1;
        %======
        % Calculation derta ECM
        if N_link_1>0  
            tmin=0.01; %%%%% 注意这里必须=0.01 ！！！
            M0=Klink*sumxi*(1+A/tmin+B/tmin^2)-Klink*sumxi_old2*(A/tmin+2*B/tmin^2)...
                +Klink*B/tmin^2*sumxi_old1;
            M1=Klink*nc_old2*A/tmin+2*Klink*nc_old2*B/tmin^2+C/tmin+2*D/tmin^2;
            M2=Klink*nc_old1*B/tmin^2+D/tmin^2;
            N=Klink*nc+Klink*nc*A/tmin+Klink*nc*B/tmin^2+C/tmin+D/tmin^2;

            ECMloc=(M0+M1*ECMloc_old2-M2*ECMloc_old1)/N;

            ECMloc_old1=ECMloc_old2;
            ECMloc_old2=ECMloc;
            sumxi_old1=sumxi_old2;
            sumxi_old2=sumxi;
            nc_old1=nc_old2;
            nc_old2=nc;
        else
            ECMloc=0;
            ECMloc_old2=0;
            ECMloc_old1=0;
            sumxi_old2=0;
            sumxi_old1=0;
            nc_old2=0;
            nc_old1=0;
        end
        %======
        % Clear data (位置很重要)
        for j=1:Num_links
            if N_links(j,1)==0
                N_links(j,2)=ECMloc; % 与基质位置相同
            end
            N_links(j,3)=0;
        end
        %======
        % Update total ECM location
        Contraction=0;
        for j=1:Num_links
            if N_links(j,1)==1
                Contraction=Contraction+Klink*(N_links(j,2)-ECMloc);
            end
        end

        TimeArry(i,1)=t;
%         N_link_0_Arry(i,1)=N_link_0;
%         N_link_1_Arry(i,1)=N_link_1;
         Contraction_Arry(i,1)=Contraction;
%         Contraction_Arry(i,1)=dint*(N_link_1/Num_links);
%         Vactin_Arry(i,1)=Vactin;
%         Averageforce_Arry(i,1)=Averageforce;

    end
    Ecm_Array_Array(KK,1)=kk;
    Contraction_Arry_Arry(KK,1)=mean(Contraction_Arry);
%     Vactin_Arry_Arry(KK,1)=mean(Vactin_Arry);
%     N_link_0_Arry_Arry(KK,1)=mean(N_link_0_Arry);
%     N_link_1_Arry_Arry(KK,1)=mean(N_link_1_Arry);
%     Averageforce_Arry_Arry(KK,1)=mean(Averageforce_Arry);
    KK=KK+1;
end

% figure (1)
% plot(TimeArry,N_link_0_Arry,'r',TimeArry,N_link_1_Arry,'b')

% figure (2)
% plot(TimeArry,Contraction_Arry)
% hold on

% figure (3)
% plot(TimeArry,Vactin_Arry)
%====================
% figure (1)
% plot(Ecm_Array_Array(1:KK-1),Contraction_Arry_Arry(1:KK-1))
% hold on

% figure (2)
% semilogx(Ecm_Array_Array(1:KK-1),Vactin_Arry_Arry(1:KK-1))

% figure (3)
% semilogx(Ecm_Array_Array(1:KK-1),Averageforce_Arry_Arry(1:KK-1))

Contraction_Arry_Arry(1:KK-1)
if Contraction_Arry_Arry(1:KK-1)<50
    Cad = 300*(1+Contraction_Arry_Arry(1:KK-1)/(10+Contraction_Arry_Arry(1:KK-1)))
else
    Cad = 300*(1-1/(10+Contraction_Arry_Arry(1:KK-1)))
end
