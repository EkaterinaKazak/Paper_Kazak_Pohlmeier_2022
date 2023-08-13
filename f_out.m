function res = f_out( rx,g,a0,a00,sl,c,a_sh,nu)
if (isempty(a00))
    a00 = a0;
end
A = (0.01:0.01:0.5)'; la = length(A);
L = 0:0.01:1; ll = length(L);

GMVP=@(x) (x\ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));
MV = @(S,m)  (inv(S)*m')/(ones(1,size(S,2))*inv(S)*m');
CE = @(x,gamma) (mean(x)-gamma*var(x));

T = size(rx,1)/2;
N = size(rx,2); wn = ones(N,1)./N;

rnx = zeros(T,1); rgx = zeros(T,1); rmx = zeros(T,1); rdx = zeros(T,1);
for j =1:T
    sx = rx(j:T + j-1,:);
    SS = cov(sx); 
    SSi = inv(SS);
    
    wg =  GMVP(SS);
    wx = ((1+rx(T+j-1,:))'.*wg)./sum(((1+rx(T + j-1,:))'.*wg));
    to_ = sum(abs(wg-wx)); rg = rx(T+j,:)*wg; rgx(j) = rg - c*to_;
    
    Sig = (1-a_sh)*SSi + a_sh*nu*eye(N);
    wd =GMVP(Sig);
    wx = ((1+rx(T+j-1,:))'.*wd)./sum(((1+rx(T + j-1,:))'.*wd));
    to_ = sum(abs(wd-wx)); rd = rx(T+j,:)*wd; rdx(j) = rd - c*to_;
    
    Sig = cov(sx) + sl*eye(N);
    wm = MV(Sig,mean(sx));
    wx = ((1+rx(T+j-1,:))'.*wm)./sum(((1+rx(T + j-1,:))'.*wm));
    to_ = sum(abs(wm-wx)); rm = rx(T+j,:)*wm; rmx(j) = rm - c*to_;
    
    wx = ((1+rx(T+j-1,:))'.*wn)./sum(((1+rx(T + j-1,:))'.*wn));
    to_ = sum(abs(wn-wx)); rn = rx(T+j,:)*wn; rnx(j) = rn - c*to_;
end
CEno = CE(rnx,g); CEgo = CE(rgx,g); CEmo = CE(rmx,g);CEdo = CE(rdx,g); MCE = [CEgo,CEdo,CEmo,CEno];
t_1e = (CEgo - CEno)/deltalw_ce(g,[rgx,rnx]); 
t_2e = (CEdo - CEno)/deltalw_ce(g,[rdx,rnx]); 
t_3e = (CEmo - CEno)/deltalw_ce(g,[rmx,rnx]); 


%Pretest with pseudo out-of-sample CE maximizing alpha: 
[~,tmax] = max([repmat(t_1e,la,1),repmat(t_2e,la,1),repmat(t_3e,la,1),norminv(1-A)],[],2);
d = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
ce_grid = MCE*d'; [~, idx] = max(ce_grid); a = A(idx); 

%Pretest with adaptively smoothed alpha maximizing pseudo out-of-sample CE: s24
ax = (a.*(1-L) + L*a00)';
[~,tmax] = max([repmat(t_1e,ll,1),repmat(t_2e,ll,1),repmat(t_3e,ll,1),norminv(1-ax)],[],2);
d = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
ce_grid = MCE*d'; [~, idx] = max(ce_grid); 
res.aL=ax(idx);
[~,tmax] = max([t_1e,t_2e,t_3e,norminv(1-ax(idx))]);
res.d7 = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];


end

