disp('The only input needed is a distance matrix file')
disp('The format of this file should be: ')
disp('Column 1: id of element i')
disp('Column 2: id of element j')
disp('Column 3: dist(i,j)')
%Change the line below to match the actual file name
mdist='s2cPROT_kept.dmtx';
xx=load(mdist);
%Set upper and lower bounds
ND=max(xx(:,2));
NL=max(xx(:,1));
if (NL>ND)
  ND=NL;
end
N=size(xx,1);
%Preallocate memory for pairwise distance matrix and density contrib matrix
dist = zeros(ND,ND,'single');
densit = zeros(ND,'single');
%Populate pairwise distances and local density contributions ("weights")
for i=1:N
  ii=xx(i,1);
  jj=xx(i,2);
  dist(ii,jj)=xx(i,3);
  dist(jj,ii)=xx(i,3);
  densit(ii)=xx(i,4);
  densit(jj)=xx(i,5);
end
percent=2.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
%Set kernel density estimator cutoff to the second percentile of the 
%distance matrix
position=round(N*percent/100);
sda=sort(xx(:,3));
dc=sda(position);

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
%Initialize local densities as density contributions from the input
rho = zeros(ND,1,'single');
for i=1:ND
  rho(i)=densit(i);
end
%
% Gaussian kernel
%Idea: For each point, compute distance to n-th nearest neighbor
%Sigma scales with number of samples (segments) - should decrease with
%segment number
%Consistency of kernel density estimates
%
for i=1:ND-1
  for j=i+1:ND
     rho(i)=rho(i)+densit(j)*exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+densit(j)*exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end
%
% "Cut off" kernel
%
%for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
%end

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
%Figure out distance to nearest neighbor of higher density
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH M(odified)')
disp('column 1:Density')
disp('column 2:Delta')
fid = fopen('DECISION_GRAPH_M', 'w');
for i=1:ND
   fprintf(fid, '%15.6f %15.6f %20.1f %20.1f\n', rho(i),delta(i),ordrho(i),nneigh(i));
end

fid = fopen('dc_value','w');
fprintf(fid, '%10.6f \n',dc);


dgraph='DECISION_GRAPH_M';
yy=load(dgraph);
for i=1:ND
  rho(i)=yy(i,1);
  delta(i)=yy(i,2);
  ordrho(i)=yy(i,3);
  nneigh(i)=yy(i,4);
end

dc=load('dc_value');

disp('Select a polygon enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end
subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')
ax = gca;

%Uncomment these lines if you'd like rho or delta on a log scale
%ax.XAxis.Scale = 'log';
%ax.YAxis.Scale = 'log';

sortdelta = sort(delta,'descend');
fig=subplot(2,1,1)
roi = drawpolygon(fig,'Color','r');
NCLUST=0;
for i=1:ND
 cl(i)=-1;
end
for i=1:ND
 [in,on]=inpolygon(rho(i),delta(i),roi.Position(:,1),roi.Position(:,2));
 if (in || on)
    NCLUST=NCLUST+1;
    cl(i)=NCLUST;
    icl(NCLUST)=i;
 end
end
%Points will be assigned as centroids if present on or in the polygon
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end
%halo
for i=1:ND
  halo(i)=cl(i);
end
if (NCLUST>1)
  for i=1:NCLUST
    bord_rho(i)=0.;
  end
  for i=1:ND-1
    for j=i+1:ND
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end
for i=1:NCLUST
  nc=0;
  nh=0;
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end

%cmap=colormap;
cmap=colormap(jet);
cmap=flip(cmap,1);
for i=1:NCLUST
   ic=int16(256. - ((255.*(i-1))/(NCLUST-1.))); %two times (edited)
   if ic == 0
       ic = 1;
   end
   %ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   rhoi = 0;
   deltai = 0;
   rhoi = double(rho(icl(i)));
   deltai = double(delta(icl(i)));
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
   text(rhoi,deltai,num2str(i),'FontSize',12,'FontWeight','bold');
end
subplot(2,1,2)



%Not super important if you don't care about the 2D projection view
%Feel free to comment the entire section out :)
disp('Performing 2D nonclassical multidimensional scaling')
%Y1 = mdscale(dist, 2, 'criterion','metricstress');
Y1 = mdscale(dist, 2, 'criterion','metricstress','start','random');
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
A(i,1)=0.;
A(i,2)=0.;
end
for i=1:NCLUST
 nn=0;
 ic=int16(256. - ((255.*(i-1))/(NCLUST-1.))); %two times (edited) 
  ic=int8((i*64.)/(NCLUST*1.));
 for j=1:ND
   if (cl(j)==i)
     nn=nn+1;
     A(nn,1)=Y1(j,1);
     A(nn,2)=Y1(j,2);
   end
 end
 hold on
 plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end
faa = fopen('CLUSTER_ASSIGNATION', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end

