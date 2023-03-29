function RRR=Optivalue(RR,k)

[m, n, c]=size(RR);
RR=floor(RR(:));
U=unique(RR);
H=histc(RR,U);
HU=sortrows([H U],1);

HHU=HU(end-k+1:end,:);


for i=1:length(RR)
    [p,q]=min(abs(RR(i)-HHU(:,2)));
    RR(i)=HHU(q,2);
end

RRR=reshape(RR, m, n, c);