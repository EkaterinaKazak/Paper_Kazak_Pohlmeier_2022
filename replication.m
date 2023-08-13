close all; clear all; clc; warning off;
spmd
    warning('off');
end

GMVP=@(x) (x\ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));
CE = @(x,gamma) (mean(x)-gamma*var(x));
CEin  = @(sample,w,g) (mean(sample)*w-g*w'*cov(sample)*w);
MV = @(S,m)  (inv(S)*m')/(ones(1,size(S,2))*inv(S)*m');

load('100.mat');
data = data(:,2:end);
Ntot = size(data,2);

g = 0.5; c = 0.005; a0 = 0.05;
T = 180; H = 500; W = size(data,1);

MC = 500;
B = 200;
NS = 12;

Ng = 10;

tbl = zeros(4,9);
for xn =1:1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = Ng(xn); 
l = ones(N,1); wn = (1/N).*l;
sl = 0.05*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ce = zeros(MC,NS); cen = ce; TO = ce; srn = ce; sr = ce;
mun = ce; mu = ce; varn = ce; va = ce;
dm = zeros(MC,H); %DeMiguel2013 shrinkage parameter
parfor m = 1:MC
    
    rng(m)
    idx = randperm(Ntot,N); r = data(:,idx);
    S = zeros(H,3,NS);
    
    a00=[]; a000 = []; d00 = ones(1,4)./4; d000 = ones(1,4)./4;
    dm_a = zeros(H,1);
    for i = 1:H
        t = i+W-T-H;
        sample = r(t:t+T-1,:);
        SS = cov(sample); %sample variance 
        SSi = inv(SS); %inverse of sample variance
        %1.  GMVP
        wg = GMVP(SS);
        r1 = r(t+T,:)*wg; wgx = ((1+r(t+T-1,:))'.*wg)./sum(((1+r(t+T-1,:))'.*wg));
        to1 = sum(abs(wg-wgx)); r1n = r1 - c*to1;
        S(i,:,1) = [r1,r1n,to1];
        
        %2. DeMiguel shrunken covariance matrix + GMVP
        Omega = (T-N-2)/((T-N-1)*(T-N-4))*(trace(SSi)*SSi +(T-N-2)*SSi*SSi);
        nu= mean(diag(SS).^(-2));
        n1= trace(Omega)-trace(SSi*SSi);
        n2 = norm(nu*eye(N)-SSi,'fro')^2;
        a_sh = n1/(n1+n2); dm_a(i) = a_sh; %shrinkage intensity
        Sig = (1-a_sh)*SSi + a_sh*nu*eye(N);
        wd =GMVP(Sig);
        r3 = r(t+T,:)*wd; wdx = ((1+r(t+T-1,:))'.*wd)./sum(((1+r(t+T-1,:))'.*wd));
        to3 = sum(abs(wd-wdx)); r3n = r3 - c*to3;
        S(i,:,2) = [r3,r3n,to3];
        
        %3. MV + shrunken covariance matrix
        Sig = SS + sl*eye(N);
        wm = MV(Sig,mean(sample));
        r3 = r(t+T,:)*wm; wmx = ((1+r(t+T-1,:))'.*wm)./sum(((1+r(t+T-1,:))'.*wm));
        to3 = sum(abs(wm-wmx)); r3n = r3 - c*to3;
        S(i,:,3) = [r3,r3n,to3];
        
        %4.  1/N
        r2 = r(t+T,:)*wn; wnx = ((1+r(t+T-1,:))'.*wn)./sum(((1+r(t+T-1,:))'.*wn));
        to2 = sum(abs(wn-wnx)); r2n = r2 - c*to2;
        S(i,:,4) = [r2,r2n,to2];
        
        ww =[wg';wd';wm';wn'];
        
        fitin = f_in( sample,ww,g,a0,a00,d00 );
        d00 = fitin.d11;
        a00 = fitin.aL;
        fitout = f_out( sample,g,a0,a000,sl,c,a_sh,nu);
        a000 = fitout.aL;
        
        %5. pretest 5%: statistical pretest rule
        d = fitin.d5;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,5) = [rp,rpn,to_];       
        
        %6. pretest a in smooth: data drive time-adaptive rule
        d = fitin.d6;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,6) = [rp,rpn,to_];
        
        %8. pretest a out smooth: data-driven time adaptive pseudo
        %out-of-sample with transaction costs
        d = fitout.d7;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,8) = [rp,rpn,to_];
        
        %12. Sequential relative performance
        d  = d00;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,12) = [rp,rpn,to_];
        
        fbag = f_bag( sample,B,ww,a00,a000,g,c);
        %7. pretest a in smooth BAGGING
        d = fbag.d15;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,7) = [rp,rpn,to_];
        
        %9. pretest a out smooth BAGGING
        d = fbag.d16;
        w = (d*ww)';
        rp = r(t+T,:)*w; wx = ((1+r(t+T-1,:))'.*w)./sum(((1+r(t+T-1,:))'.*w));
        to_ = sum(abs(w-wx)); rpn = rp - c*to_;
        S(i,:,9) = [rp,rpn,to_];
        
        %10. Kirby Ostdiek VT
        wm =(var(sample)./sum(var(sample)))';
        r3 = r(t+T,:)*wm; wmx = ((1+r(t+T-1,:))'.*wm)./sum(((1+r(t+T-1,:))'.*wm));
        to3 = sum(abs(wm-wmx)); r3n = r3 - c*to3;
        S(i,:,10) = [r3,r3n,to3];
        
        %11. Kirby Ostdiek RRT
        wm =(mean(sample)./var(sample)/(sum(mean(sample)./var(sample))))';
        r3 = r(t+T,:)*wm; wmx = ((1+r(t+T-1,:))'.*wm)./sum(((1+r(t+T-1,:))'.*wm));
        to3 = sum(abs(wm-wmx)); r3n = r3 - c*to3;
        S(i,:,11) = [r3,r3n,to3];
        
       
       
    end
    R = reshape(S(:,1,:),H,NS); R_net = reshape(S(:,2,:),H,NS); tover = reshape(S(:,3,:),H,NS);
    ce(m,:) = CE(R,g); cen(m,:) =  CE(R_net,g); TO(m,:) = mean(tover,1);
    srn(m,:) = mean(R_net)./std(R_net);
    sr(m,:) = mean(R)./std(R);
    mun(m,:) = mean(R_net); varn(m,:) = var(R_net);
    mu(m,:) = mean(R); va(m,:) = var(R);
    dm(m,:) = dm_a;
    
end

rowNames = {'gmvp','D2013+GMVP','mv-sh','1/N','5% in', 'a in sm','a in bag','a out sm','a out bag','VT','RRT','SP'};
colNames = {'CEn','SRn','CE','TO','sr','mu_net','var_net','mu','var'};

MTX = [mean(cen,1)',mean(srn,1)',mean(ce)',mean(TO)',mean(sr)',mean(mun)',mean(varn)',mean(mu)',mean(va)'];
TAB = array2table(MTX,'RowNames',rowNames,'VariableNames',colNames);

nm = {['N',num2str(N),'T',num2str(T),'.mat']};
save(char(nm))
clc
end








