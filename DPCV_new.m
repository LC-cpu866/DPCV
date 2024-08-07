% 基于数据间方差的密度峰值聚类算法
% Density peak clustering algorithm based on variance between data
% DPCV下，无需归一化的数据：
%         Pathbased, Jain, Aggregation
% MDDPC下，无需归一化的数据：
%         Pathbased, Jain, Aggregation
clear
clc

t0=cputime;

% load D:\Xnewm\datasets\synthetic\forty.mat;
data=xlsread('D:\Xm\UCI\Letter\Letter.xlsx');
labels=xlsread('D:\Xm\UCI\Letter\label.xlsx');

data = normalize(data, "range");

e0=cputime-t0;
t1=cputime;

ND=size(data,1);
N=ND*(ND-1)/2;

percent = 0.002;
K = 1;

% distM=zeros(ND);
distK=zeros(ND,K);


% for i=1:ND-1
%     for j=i+1:ND
%         distM(i,j)=norm(data(i,:)-data(j,:));
%         distM(j,i)=distM(i,j);
%     end
% end
distM = pdist2(data, data);

maxd=max(max(distM));

for i=1:ND
    [~,orddist]=sort(distM(i,:));
    distK(i,:)=orddist(2:K+1);
end

e1=cputime-t1;
t2=cputime;

position=round(N*percent/100);
sda=sort((distM(triu(true(size(distM)),1)))');
dc=sda(position);
clear sda;

for i=1:ND
    rho(i)=1;
end

for i=1:ND-1
    for j=i+1:ND
        %stdM=std(data(i,:)-data(j,:))/2+1; % DPCV
        stdM=norm(data(i,:)-data(j,:),1)/2+1; % MDDPC
        rho(i)=rho(i)+exp(-(distM(i,j)/(stdM*dc))^2);
        rho(j)=rho(j)+exp(-(distM(i,j)/(stdM*dc))^2);
    end
end

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii=2:ND
    delta(ordrho(ii))=maxd;
    for jj=1:ii-1
        if(distM(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
            delta(ordrho(ii))=distM(ordrho(ii),ordrho(jj));
            nneigh(ordrho(ii))=ordrho(jj);
        end
    end
end
delta(ordrho(1))=max(delta(:));

e2=cputime-t2;

figure(1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('Decision Graph','FontSize',15.0);
xlabel('\rho');
ylabel('\delta');

NCLUST=length(unique(labels));

t3=cputime;

for i=1:ND
    cl(i)=-1;
end

r=rho.*delta;
[r_sorted,ordr]=sort(r,'descend');
r(ordr(NCLUST));
k=0;
for i=1:ND
    if (rho(i)*delta(i)>=r(ordr(NCLUST)))
        k=k+1;
        cl(i)=k;
        icl(k)=i;
    end
end

R=icl;
for i=1:length(icl)
    R(i)=distM(icl(i),distK(icl(i),K));
end

%进行NCLUST轮扫描
a=0.4;
[~,ordicl]=sort(rho(icl),'descend');
for i=1:NCLUST
    Ri=R(ordicl(i));
    Ri=Ri+0.0001;
    Ci=icl(ordicl(i));
    index=find(distM(Ci,:)<=Ri);
    judge=rho(Ci);
    index=intersect(index,find(cl<0));
    queue=index;
    if length(queue)==0
        continue;
    end
    queue(1)=[];
    cl(index)=cl(Ci);
    while(~isempty(queue))
        Ci=queue(1);
        queue(1)=[];
        if (rho(Ci)<=a*judge)
            continue;
        end
        index=find(distM(Ci,:)<=Ri);
        index=intersect(index,find(cl<0));
        cl(index)=cl(Ci);
        queue=[queue index];
    end
end

for i=1:ND
    if(cl(ordrho(i))==-1)
        cl(ordrho(i))=cl(nneigh(ordrho(i)));
    end
end

halo=cl;

e3=cputime-t3;

cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   figure(1);
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

for i=1:ND
    A(i,1)=0;
    A(i,2)=0;
end

nn=0;
ic=int8((3*64.)/((NCLUST+1)*1.));
for j=1:ND
    if (labels(j)==1)
        nn=nn+1;
        A(nn,1)=rho(j);
        A(nn,2)=delta(j);
    end
end
hold on
plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));

figure(2);

for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=data(j,1);
      A(nn,2)=data(j,2);
    end
  end
  hold on
  switch mod(i,14)+1
      case 1
          plot(A(1:nn,1),A(1:nn,2),'c<');%,'MarkerFaceColor','c');
      case 2
          plot(A(1:nn,1),A(1:nn,2),'ms');%,'MarkerFaceColor','m');
      case 3
          plot(A(1:nn,1),A(1:nn,2),'b^');%,'MarkerFaceColor','y');
      case 4
          plot(A(1:nn,1),A(1:nn,2),'g<');%,'MarkerFaceColor','r');
      case 5
          plot(A(1:nn,1),A(1:nn,2),'gs');%,'MarkerFaceColor','g');
      case 6
          plot(A(1:nn,1),A(1:nn,2),'bp');%,'MarkerFaceColor','b');
      case 7
          plot(A(1:nn,1),A(1:nn,2),'kh');%,'MarkerFaceColor','w');
      case 8
          plot(A(1:nn,1),A(1:nn,2),'kp');%,'MarkerFaceColor','c');
      case 9
          plot(A(1:nn,1),A(1:nn,2),'bh');%,'MarkerFaceColor','m');
      case 10
          plot(A(1:nn,1),A(1:nn,2),'g+');%,'MarkerFaceColor','y');
      case 11
          plot(A(1:nn,1),A(1:nn,2),'ro');%,'MarkerFaceColor','r');
      case 12
          plot(A(1:nn,1),A(1:nn,2),'y^');%,'MarkerFaceColor','g');
      case 13
          plot(A(1:nn,1),A(1:nn,2),'mx');%,'MarkerFaceColor','b');
      case 14
          plot(A(1:nn,1),A(1:nn,2),'cs');%,'MarkerFaceColor','w');
  end
  %plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
  
  plot(data(icl,1),data(icl,2),'*','MarkerSize',15);
end

xlabel('属性1');
ylabel('属性2');

(ND*ND-ND)/2
e1
tt=e1+e2+e3;
evaluation_index(labels, cl);