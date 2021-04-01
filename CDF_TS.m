function [ data ] = CDF_TS( data,eta,thre ) 

gap=99999;
t=0;
[n,d]=size(data);
while gap>thre
    odata=normalise(data);
    [ ndata ] = CDF_TS1( odata,eta ); %FScaling CDF_TS
    data=normalise(ndata);
    gap=sum(sum(abs(data-odata)))/n/d;
    t=t+1;
end
X = ['Total iteration time is ',num2str(t)];
disp(X)
end

function [ Ndata ] = CDF_TS1( data,eta ) 
data = normalise(data);
dis = pdist2(data,data,'minkowski',2); 
[n,dim] = size(data);
[ Ndis, dis ] = DScale( dis, eta, dim );
R = Ndis./dis;
R(dis==0) = 1;
Ndata = data .* repmat(sum(R, 1)', [1, dim])...
- reshape(sum(reshape(data, [n, 1, dim]).*reshape(R, [n, n, 1]), 1), [n, dim])...
    + repmat(sum(data, 1), [n, 1]); 
Ndata = normalise(Ndata);
end

function [ Ndis, dis ] = DScale( dis, eta, dim )
n = size(dis, 1);
dis = dis./max(max(dis)); 
Ndis = dis;
Lb = Ndis <= eta;
Rb = ~Lb;
den = sum(Lb, 2);
Rate = (den./n).^(1/dim)./eta;
rep = repmat(Rate, [1, n]);
Ndis(Lb) = Ndis(Lb) .* rep(Lb);
rep = rep .* eta;  
Ndis(Rb) = (Ndis(Rb) - eta).* (1 - rep(Rb))./(1-eta) + rep(Rb);
end

function [ ndata ] = normalise( data )
n = size(data, 1);
l_min = min(data);
ndata = (data - repmat(l_min, [n, 1]))./repmat(max(data) - l_min, [n, 1]);
end