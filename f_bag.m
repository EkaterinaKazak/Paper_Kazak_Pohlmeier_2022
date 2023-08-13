function res = f_bag( sample,B,ww,a_in,a_out,g,c)
CEin  = @(sample,w,g) (mean(sample)*w-g*w'*cov(sample)*w);


CE = @(x,gamma) (mean(x)-gamma*var(x));

T = size(sample,1); T2 = T/2; 
wgb = ww(1,:)'; wdb = ww(2,:)'; wmb = ww(3,:)'; wn = ww(3,:)';

d15b = zeros(B,4); d16b = d15b;
for b = 1:B
    rng(b)
    J = ceil(rand(T,1)*T);
    sb = sample(J,:);
    rp_inb = [sb*wgb,sb*wdb, sb*wmb,sb*wn];
    CEginb = CEin(sb,wgb,g);  CEdinb = CEin(sb,wdb,g); CEminb = CEin(sb,wmb,g); CEninb = CEin(sb,wn,g);
    t_1eb = (CEginb - CEninb)/deltalw_ce(g,rp_inb(:,[1,4]));
    t_2eb = (CEdinb - CEninb)/deltalw_ce(g,rp_inb(:,[2,4]));
    t_3eb = (CEminb - CEninb)/deltalw_ce(g,rp_inb(:,[3,4]));
    
    [~,tmax] = max([t_1eb,t_2eb,t_3eb,norminv(1-a_in)]);
    d15b(b,:) = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
    
    rnx = zeros(T2,1); rgx = zeros(T2,1); rmx = zeros(T2,1); rdx = zeros(T2,1);
    for j =1:T2
        wg =  wgb;
        wx = ((1+sb(T2+j-1,:))'.*wg)./sum(((1+sb(T2 + j-1,:))'.*wg));
        to_ = sum(abs(wg-wx)); rg = sb(T2+j,:)*wg; rgx(j) = rg - c*to_;
        
        wd = wdb;
        wx = ((1+sb(T2+j-1,:))'.*wd)./sum(((1+sb(T2 + j-1,:))'.*wd));
        to_ = sum(abs(wd-wx)); rd = sb(T2+j,:)*wd; rdx(j) = rd - c*to_;

        wm = wmb;
        wx = ((1+sb(T2+j-1,:))'.*wm)./sum(((1+sb(T2 + j-1,:))'.*wm));
        to_ = sum(abs(wm-wx)); rm = sb(T2+j,:)*wm; rmx(j) = rm - c*to_;
        
        wx = ((1+sb(T2+j-1,:))'.*wn)./sum(((1+sb(T2 + j-1,:))'.*wn));
        to_ = sum(abs(wn-wx)); rn = sb(T2+j,:)*wn; rnx(j) = rn - c*to_;
    end
    CEno = CE(rnx,g); CEgo = CE(rgx,g); CEdo = CE(rmx,g);CEmo = CE(rmx,g); MCE = [CEgo,CEdo,CEmo,CEno];
    t_1e = (CEgo - CEno)/deltalw_ce(g,[rgx,rnx]);
    t_2e = (CEdo - CEno)/deltalw_ce(g,[rdx,rnx]);
    t_3e = (CEmo - CEno)/deltalw_ce(g,[rmx,rnx]);
    
    [~,tmax] = max([t_1e,t_2e,t_3e,norminv(1-a_out)]);
    d16b(b,:) = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
    
end

res.d15 = mean(d15b);
res.d16 = mean(d16b);
end

