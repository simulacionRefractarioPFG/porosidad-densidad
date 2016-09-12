function [e, rho, mu] = PostProcDump(filename)
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lee archivo con la informacion de las particulas de la simulacion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DUMP     = dlmread(filename,' ',9,0); % skip 9 rows, 0 columns

material   = DUMP(:,2);         % 1->MgO ; 2->AL2O3
x          = DUMP(:,3);         % posicion x del centro de la particula
y          = DUMP(:,4);         % posicion y del centro de la particula
z          = DUMP(:,5);         % posicion z del centro de la particula
r          = DUMP(:,18);        % radio de la particula
z_max      = z + r;             % altura maxima de la particula
z_min      = z - r;             % altura minima de particula
Nparticles = length(x);         % numero total de particulas
Vol        = (4/3)*pi*r.^3;     % volumen de la particula
d          = sqrt(x.^2 + y.^2); % distancia al centro
d_max      = d + r;             % distancia maxima a (0,0) de la particula
d_min      = d - r;             % distancia minima a (0,0) de la particula 

rho_MgO    = 3.5;               % pg*nu_m-3
rho_Al2O3  = 3;                 % pg*nu_m-3

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dimensiones que definen la cama de particulas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Limite superior representativo del disco de particulas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La altura maxima media del 20% de las particulas
z_max_sorted = sort(z_max);
top10percent = z_max_sorted((end-ceil(0.2*Nparticles)):end);
z_top        = sum(top10percent)/(length(top10percent));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Limite inferior representativo del disco de particulas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La altura minima media del 20% de las particulas
z_min_sorted = sort(z_min);
inf10percent = z_min_sorted(1:ceil((0.2*Nparticles)));
z_inf        = sum(inf10percent)/(length(inf10percent));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Radio representativo del disco de particulas           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La distancia al centro maxima media del 30% de las particulas
d_max_sorted    = sort(d_max);
d_max_10percent = d_max_sorted((end-ceil(0.3*Nparticles)):end);
R               = sum(d_max_10percent)/(length(d_max_10percent));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas dentro del disco representativo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = z_top-z_inf;
Vd = pi*R^2*h;

Vp      = sum(Vol((d_max <= R) & (z_max <= z_top) & (z_min >= z_inf)));
V_MgO   = sum(Vol((d_max <= R) & (z_max <= z_top) & (z_min >= z_inf) & (material==1)));
V_Al2O3 = sum(Vol((d_max <= R) & (z_max <= z_top) & (z_min >= z_inf) & (material==2)));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas justo en las fronteras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas que se han visto cortadas por el limite superior de z_top
%%% tienen parte de su geometria dentro del disco. Considera que parte
%%% de su volumen y masa entra dentro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(x)
    if((z_max(i)>z_top) && (z_min(i)<z_top) && (d_max(i)<R))
        h_cap = z_top - z_min(i);
        V_cap = 0;
        if(h_cap<r(i))
            V_cap = pi*h_cap^2*(3*r(i)-h_cap)/3;
        else
            h_cap = 2*r(i) - h_cap;
            V_cap = (4/3)*pi*r(i)^3 - pi*h_cap^2*(3*r(i)-h_cap)/3;
        end
        ratioCap = V_cap/(4/3*pi*r(i)^3);
        % if(ratioCap>0)
        %     ratioCap
        %     viscircles([d(i),z(i)], r(i),'Color','b');
        %     hold on
        % end

        if(material(i)==1)
            V_MgO   = V_MgO + V_cap;
        else
            V_Al2O3 = V_Al2O3 + V_cap;
        end
        Vp    = Vp + V_cap;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas que se han visto cortadas por el limite inferior de z_inf
%%% tienen parte de su geometria dentro del disco. Considera que parte
%%% de su volumen y masa entra dentro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(x)
    if((z_max(i)>z_inf) && (z_min(i)<z_inf) && (d_max(i)<R))
        h_cap = z_max(i) - z_inf;
        V_cap = 0;
        if(h_cap<r(i))
            V_cap = pi*h_cap^2*(3*r(i)-h_cap)/3;
        else
            h_cap = 2*r(i) - h_cap;
            V_cap = (4/3)*pi*r(i)^3 - pi*h_cap^2*(3*r(i)-h_cap)/3;
        end
        ratioCap = V_cap/(4/3*pi*r(i)^3);
        % if(ratioCap>0)
        %     ratioCap
        %     viscircles([d(i),z(i)], r(i),'Color','b');
        %     hold on
        % end

        if(material(i)==1)
            V_MgO   = V_MgO + V_cap;
        else
            V_Al2O3 = V_Al2O3 + V_cap;
        end
        Vp    = Vp + V_cap;
    end
end


function fl1 = fl1(zs,i)
    xs = sqrt(r(i)^2-zs.^2-((r(i)^2+d(i)^2-R^2-zs.^2)/(2*d(i))).^2);
    fl1 = -(d(i).*xs)+...
         (((r(i)^2-zs.^2)/2).*asin((xs)./(sqrt(r(i)^2-zs.^2)))+(((xs)/2).*sqrt((r(i)^2-zs.^2)-(xs).^2)))+...
         (((R^2)/2).*asin((xs)./(R))+(((xs)/2).*sqrt(R^2-(xs).^2)));
end
function fl2 = fl2(zs,i)
    xs = sqrt(r(i)^2-zs.^2-((r(i)^2+d(i)^2-R^2-zs.^2)/(2*d(i))).^2);
    fl2 = (d(i).*xs)+...
         (((r(i)^2-zs.^2)/2).*asin((xs)./(sqrt(r(i)^2-zs.^2)))+(((xs)/2).*sqrt((r(i)^2-zs.^2)-(xs).^2)))-...
         (((R^2)/2).*asin((xs)./(R))+(((xs)/2).*sqrt(R^2-(xs).^2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas que se han visto cortadas por el limite lateral de
%%% Diam_disco tienen parte de su geometria dentro del disco.
%%% Considera que parte de su volumen y masa entra dentro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(x)
   if((d_min(i)<R) && (d_max(i)>R) && (sqrt((d(i)-R)^2+(z(i)-z_top)^2)>r(i)) &&...
      (sqrt((d(i)-R)^2+(z(i)-z_inf)^2)>r(i)) && (z(i)<z_top) && (z(i)>z_inf))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        h_cap = R -(d_min(i));
%        V_cap = 0;
%        if(h_cap<r(i))
%            V_cap = (1/3)*pi*h_cap^2*(3*r(i)-h_cap);
%        else
%            h_cap = 2*r(i) - h_cap;
%            V_cap = (4/3)*pi*r(i)^3 - pi*h_cap^2*(3*r(i)-h_cap)/3;
%        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V_cap = 0;

        % Corte lateral
        if(d(i)>=R)
            % Integral con trapecios
            z1 = (sqrt(r(i)^2-(d(i)-R)^2)); 
            z0 = 0;
            Nintervalos = 100000;
            h = (z1-z0)/Nintervalos;
            j = linspace(z0,z1,Nintervalos);
            % Multiplicado x4 por la simetria de la integral en X y Z
            V_cap = real(4*h*sum((fl1(j(1:end-1),i)+fl1(j(2:end),i))/2));

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap>0.45)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end

        elseif(d(i)<R)
            % Integral con trapecios
            z1 = (sqrt(r(i)^2-(d(i)-R)^2)); 
            z0 = 0;
            Nintervalos = 100000;
            h = (z1-z0)/Nintervalos;
            j = linspace(z0,z1,Nintervalos);
            % Multiplicado x4 por la simetria de la integral en X y Z
            V_cap = 4/3*pi*r(i)^3 - real(4*h*sum((fl2(j(1:end-1),i)+fl2(j(2:end),i))/2));

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap>0.999)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end
        end

        % si tambien corta con la base superior
        if(z_max(i)>z_top && d(i)<R)              
            h_capSuperior = z_max(i) - z_top;
            V_capSuperior= pi*h_capSuperior^2*(3*r(i)-h_capSuperior)/3; 
            V_cap = V_cap - V_capSuperior;

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap>0.96)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end
        end
        % si tambien corta con la base inferior
        if(z_min(i)<z_inf && d(i)<R)          
            h_capInferior = z_inf - z_min(i);
            V_capInferior = pi*h_capInferior^2*(3*r(i)-h_capInferior)/3;
            V_cap = V_cap - V_capInferior;  

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap>0.96)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end
        end

        % ratioCap = V_cap/(4/3*pi*r(i)^3);
        % if(ratioCap<0.51 && ratioCap>0.49)
        %     ratioCap
        %     viscircles([d(i),z(i)], r(i),'Color','b');
        %     hold on
        % end

        if(material(i)==1)
            V_MgO   = V_MgO + V_cap;
        else
            V_Al2O3 = V_Al2O3 + V_cap;
        end
        Vp    = Vp + V_cap;
    end
end


function f1 = f1(xs,i,a_2,z_top)
    f1 = (z_top-z(i))*((sqrt(a_2^2-xs.^2))-(d(i)-sqrt(R^2-xs.^2)))+...
        (((r(i)^2-xs.^2)/2).*asin((sqrt(a_2^2-xs.^2))/(sqrt(r(i)^2-xs.^2)))+...
        (((sqrt(a_2^2-xs.^2))/2).*(sqrt((r(i)^2-xs.^2)-(a_2^2-xs.^2)))))-...
        (((r(i)^2-xs.^2)/2).*asin((d(i)-sqrt(R^2-xs.^2))/(sqrt(r(i)^2-xs.^2)))+...
        (((d(i)-sqrt(R^2-xs.^2))/2).*(sqrt((r(i)^2-xs.^2)-(d(i)-sqrt(R^2-xs.^2)).^2))));
end
function f2 = f2(z,i)
    xs = sqrt(r(i)^2-z.^2-((r(i)^2+d(i)^2-R^2-z.^2)./(2*d(i))).^2);
    f2 = -(d(i).*xs)+...
         (((r(i)^2-z.^2)/2).*asin((xs)./(sqrt(r(i)^2-z.^2)))+(((xs)/2).*sqrt((r(i)^2-z.^2)-(xs).^2)))+...
         (((R^2)/2).*asin((xs)./(R))+(((xs)/2).*sqrt(R^2-(xs).^2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas que se han visto cortadas por un limite superior y el
%%% limite lateral de Diam_disco tienen parte de su geometria dentro del disco.
%%% Considera que parte de su volumen y masa entra dentro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(x)
    if(sqrt((d(i)-R)^2+(z(i)-z_top)^2)<r(i))
        V_cap = 0;
        if(z_top<z(i))
            h_2 = z_top - z_min(i);
            a_2 = sqrt(2*r(i)*h_2 - h_2^2); 
            % Integral con trapecios
            x_sup = sqrt(a_2^2-((d(i)^2+a_2^2-R^2)/(2*d(i)))^2); 
            x_inf = 0;
            Nintervalos = 100000;
            h = (x_sup-x_inf)/Nintervalos;
            j = linspace(x_inf,x_sup,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = 2*h*sum((f1(j(1:end-1),i,a_2,z_top)+f1(j(2:end),i,a_2,z_top))/2);

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap>0.2)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end 
        elseif(z_top>=z(i))
            %%%%%%%%%%%%%%%%%%%
            % Volumen hasta z=0
            %%%%%%%%%%%%%%%%%%%
            % Integral con trapecios
            z1 = (z_top-z(i)); 
            z0 = 0;
            Nintervalos = 100000;
            h = (z1-z0)/Nintervalos;
            j = linspace(z0,z1,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = 2*h*sum((f2(j(1:end-1),i)+f2(j(2:end),i))/2);

            %%%%%%%%%%%%%%%%%%
            % Volumen restante 
            %%%%%%%%%%%%%%%%%%
            a_2 = r(i)*0.999999; % exactamente r(i) introduce imaginarios
            % Integral con trapecios
            x_sup = sqrt(a_2^2-((d(i)^2+a_2^2-R^2)/(2*d(i)))^2); 
            x_inf = 0;
            Nintervalos = 100000;
            h = (x_sup-x_inf)/Nintervalos;
            j = linspace(x_inf,x_sup,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = V_cap + 2*h*sum((f1(j(1:end-1),i,a_2,z(i))+f1(j(2:end),i,a_2,z(i)))/2);

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(sqrt((d(i)-R)^2+(z(i)-z_top)^2)>0.99*r(i))
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end 
        end

        % if(ratioCap<1 && ratioCap>0.6)
        %     ratioCap
        %     viscircles([d(i),z(i)], r(i),'Color','b');
        %     hold on
        % end 

        if(material(i)==1)
            V_MgO   = V_MgO + V_cap;
        else
            V_Al2O3 = V_Al2O3 + V_cap;
        end

        Vp    = Vp + V_cap;
    end
end


function f3 = f3(xs,i,a_2,z_inf)
    f3 = -(z_inf-z(i))*((sqrt(a_2^2-xs.^2))-(d(i)-sqrt(R^2-xs.^2)))+...
        (((r(i)^2-xs.^2)/2).*asin((sqrt(a_2^2-xs.^2))/(sqrt(r(i)^2-xs.^2)))+...
        (((sqrt(a_2^2-xs.^2))/2).*(sqrt((r(i)^2-xs.^2)-(a_2^2-xs.^2)))))-...
        (((r(i)^2-xs.^2)/2).*asin((d(i)-sqrt(R^2-xs.^2))/(sqrt(r(i)^2-xs.^2)))+...
        (((d(i)-sqrt(R^2-xs.^2))/2).*(sqrt((r(i)^2-xs.^2)-(d(i)-sqrt(R^2-xs.^2)).^2))));
end
function f4 = f4(zs,i)
    xs = sqrt(r(i)^2-zs.^2-((r(i)^2+d(i)^2-R^2-zs.^2)./(2*d(i))).^2);
    f4 = -(d(i).*xs)+...
         (((r(i)^2-zs.^2)/2).*asin((xs)./(sqrt(r(i)^2-zs.^2)))+(((xs)/2).*sqrt((r(i)^2-zs.^2)-(xs).^2)))+...
         (((R^2)/2).*asin((xs)./(R))+(((xs)/2).*sqrt(R^2-(xs).^2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particulas que se han visto cortadas por un limite inferior y el
%%% limite lateral de Diam_disco tienen parte de su geometria dentro del disco.
%%% Considera que parte de su volumen y masa entra dentro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(x)
    if(sqrt((d(i)-R)^2+(z(i)-z_inf)^2)<r(i))
        V_cap = 0;
        if(z_inf>z(i))
            h_2 = z_max(i) - z_inf;
            a_2 = sqrt(2*r(i)*h_2 - h_2^2); 
            % Integral con trapecios
            x_sup = sqrt(a_2^2-((d(i)^2+a_2^2-R^2)/(2*d(i)))^2); 
            x_inf = 0;
            Nintervalos = 100000;
            h = (x_sup-x_inf)/Nintervalos;
            j = linspace(x_inf,x_sup,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = 2*h*sum((f3(j(1:end-1),i,a_2,z_inf)+f3(j(2:end),i,a_2,z_inf))/2);

            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap<0.01)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end 
        elseif(z_inf<=z(i))
            %%%%%%%%%%%%%%%%%%%
            % Volumen hasta z=0
            %%%%%%%%%%%%%%%%%%%
            % Integral con trapecios
            z1 = 0; 
            z0 = (z_inf - z(i));
            Nintervalos = 100000;
            h = (z1-z0)/Nintervalos;
            j = linspace(z0,z1,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = 2*h*sum((f4(j(1:end-1),i)+f4(j(2:end),i))/2);

            %%%%%%%%%%%%%%%%%%
            % Volumen restante 
            %%%%%%%%%%%%%%%%%%
            a_2 = r(i)*0.999999; % exactamente r(i) introduce imaginarios
            % Integral con trapecios
            x_sup = sqrt(a_2^2-((d(i)^2+a_2^2-R^2)/(2*d(i)))^2); 
            x_inf = 0;
            Nintervalos = 100000;
            h = (x_sup-x_inf)/Nintervalos;
            j = linspace(x_inf,x_sup,Nintervalos);
            % Multiplicado x2 por la simetria de la integral en x
            V_cap = V_cap + 2*h*sum((f3(j(1:end-1),i,a_2,z(i))+f3(j(2:end),i,a_2,z(i)))/2);
            
            % ratioCap = V_cap/(4/3*pi*r(i)^3);
            % if(ratioCap<0.02)
            %     ratioCap
            %     viscircles([d(i),z(i)], r(i),'Color','b');
            %     hold on
            % end 
        end

        % if(ratioCap>0)
        %     ratioCap
        %     viscircles([d(i),z(i)], r(i),'Color','b');
        %     hold on
        % end 

        if(material(i)==1)
            V_MgO   = V_MgO + V_cap;
        else
            V_Al2O3 = V_Al2O3 + V_cap;
        end

        Vp    = Vp + V_cap;
    end
end

% plot(0:12000,z_top*ones(12001,1),'b')
% plot(0:12000,z_inf*ones(12001,1),'b')
% hold on
% plot(R*ones(10001,1),0:10000,'g')
% hold on
% axis([R-1000 R+1000 z_top-1000 z_top+1000])
% % axis([R-1000 R+1000 z_inf-1000 z_inf+1000])
% xlabel('distancia al centro')
% ylabel('altura')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Porosidad, densidad y permeabilidad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Porosidad, fraccion de huecos
e   = (Vd-Vp)/Vd; 
% Densidad del disco; multiplicando por 1000 para pasar de micro a si
M   = V_MgO*rho_MgO + V_Al2O3*rho_Al2O3;
rho = (M/Vd)*1000; 
% Permeabilidad
mu  = 1;           

end