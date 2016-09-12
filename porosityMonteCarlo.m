function porosity = porosityMonteCarlo(filename, Ntrials)

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


% Porosidad
Nvoid = 0;
voidRatio = zeros(Ntrials,1);
puntosInteriores = [];
puntosExteriores = [];
for trials = 1:Ntrials
	theta = (2*pi)*rand(1);
	% importante el siguiente sqrt() para distribucion uniforme
	r_sam = R*sqrt(rand(1)); 
	x_sam = r_sam*cos(theta);
	y_sam = r_sam*sin(theta);
	z_sam = (z_top-z_inf)*rand(1) + z_inf;

	% el punto esta a menos de un radio de alguna particula?
	touchVec = sqrt((x_sam-x).^2+(y_sam-y).^2+(z_sam-z).^2) < r;
	if(any(touchVec))
		% el punto ha caido en alguna particula
		puntosInteriores = [puntosInteriores;x_sam,y_sam,z_sam];
	else
		Nvoid = Nvoid+1;
		puntosExteriores = [puntosExteriores;x_sam,y_sam,z_sam];
	end
	voidRatio(trials) = Nvoid/trials;
end

plot3(puntosInteriores(:,1),puntosInteriores(:,2),puntosInteriores(:,3),'b*');
hold on
plot3(puntosExteriores(:,1),puntosExteriores(:,2),puntosExteriores(:,3),'r*');
axis equal
axis off; 
set(gca,'color','none');
figure
plot(1:Ntrials,voidRatio);
xlabel('#samples')
ylabel('\epsilon')
porosity = voidRatio(end);
end