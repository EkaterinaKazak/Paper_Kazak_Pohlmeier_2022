function res = f_in( sample,ww,g,a0,a00,d00 )
if (isempty(a00))
    a00 = a0;
end
A = (0.01:0.01:0.5)'; la = length(A);
L = 0:0.01:1; ll = length(L);

CEin  = @(sample,w,g) (mean(sample)*w-g*w'*cov(sample)*w);
wg = ww(1,:)'; wd = ww(2,:)'; wm = ww(3,:)'; wn = ww(4,:)';

rp_in = [sample*wg,sample*wd,sample*wm,sample*wn];
CEgin = CEin(sample,wg,g); CEdin = CEin(sample,wd,g); CEmin = CEin(sample,wm,g);CEnin = CEin(sample,wn,g);
t_1e = (CEgin - CEnin)/deltalw_ce(g,rp_in(:,[1,4])); 
t_2e = (CEdin - CEnin)/deltalw_ce(g,rp_in(:,[2,4]));
t_3e = (CEmin - CEnin)/deltalw_ce(g,rp_in(:,[3,4])); 


% Pretest 5%
[~,tmax] = max([t_1e,t_2e,t_3e,norminv(1-a0)]);
res.d5 = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];

% Pretest a in
[~,tmax] = max([repmat(t_1e,la,1),repmat(t_2e,la,1),repmat(t_3e,la,1),norminv(1-A)],[],2);
d = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
w = (d*ww)';
ce_grid = (mean(sample)*w)'-g*diag(w'*cov(sample)*w);
[~, idx] = max(ce_grid);  a = A(idx);

% Pretest a in smooth
ax = (a.*(1-L) + L*a00)';
[~,tmax] = max([repmat(t_1e,ll,1),repmat(t_2e,ll,1),repmat(t_3e,ll,1),norminv(1-ax)],[],2);
d = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];
w = (d*ww)';
ce_grid = (mean(sample)*w)'-g*diag(w'*cov(sample)*w);
[~, idx] = max(ce_grid); 
res.aL=ax(idx);
[~,tmax] = max([t_1e,t_2e,t_3e,norminv(1-ax(idx))]);
res.d6 = [(tmax==1),(tmax==2),(tmax==3),(tmax==4)];

% sequential performance
res.d11 = d00.*exp([CEgin,CEdin,CEmin,CEnin])./sum(d00.*exp([CEgin,CEdin,CEmin,CEnin]));


end

