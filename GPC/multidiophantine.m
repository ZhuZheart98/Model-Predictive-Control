function  [E,F,G] = multidiophantine(a,b,c,N)  %d = 1

na = length(a)-1;nb = length(b)-1;nc = length(c)-1;%A、B、C的阶次?

%E、F、G的初值?
E = zeros(N);E(1,1) = 1;F(1,:) = conv(b,E(1,:));
if na>=nc
    G(1,:) = [c(2:nc+1) zeros(1,na-nc)]-a(2:na+1); %令c(nc+2)=c(nc+3)=...=0?
else
    G(1,:) = c(2:nc+1) -[a(2:na+1)-zeros(1,nc-na)]; %令a(nc+2)=a(nc+3)=...=0          
end
%求E、F、G
for j = 1:N-1
    for i = 1:j
        E(j+1,i) = E(j,i);
    end
    E(j+1,i+1) = G(j,1);
        for i = 2:na
            G(j+1,i-1) = G(j,i)-G(j,1)*a(i);
        end
        G(j+1,na) = -G(j,1)*a(na+1);
        F(j+1,:) = conv(b,E(j+1,:));
end