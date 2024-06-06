rm(list=ls())

f.w <- function(p,g,d){
    (d*p^g)/(d*p^g+(1-p)^g);
};

fcd <- function(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc) {
    pre.a   <- matrix(v[1:nc],k.u,nc);
    pre.b   <- matrix(v[(nc+1):(2*nc)],1,nc);
    pre.g.g <- matrix(v[(2*nc+1):(3*nc)],k.v,nc);
    pre.l.g <- matrix(v[(3*nc+1):(4*nc)],1,nc);
    pre.g.d <- matrix(v[(4*nc+1):(5*nc)],k.w,nc);
    pre.l.d <- matrix(v[(5*nc+1):(6*nc)],1,nc);
    g.sigma  <- matrix(v[(6*nc+1):(6*nc+j)], g.n, j, byrow=T);
    l.sigma  <- matrix(v[(6*nc+j+1):(6*nc+2*j)], l.n, j, byrow=T);
    ret <- matrix(0,j,nc);
    for (c in 1:nc){
	a   <- matrix(pre.a[,c],g.n,j);
	b   <- matrix(pre.b[,c],l.n,j);
	g.g <- matrix(pre.g.g[,c],g.n,j);
	l.g <- matrix(pre.l.g[,c],l.n,j);
	g.d <- matrix(pre.g.d[,c],g.n,j);
	l.d <- matrix(pre.l.d[,c],l.n,j);
	g.cehat <- (f.w(g.p,g.g,g.d)*g.z1^a + (1-f.w(g.p,g.g,g.d))*g.z2^a)^(1/a);
	g.ret   <- dnorm(g.ce, g.cehat, g.sigma*(g.z1-g.z2))^g.sel;
	l.cehat <- -((1-f.w(1-l.p,l.g,l.d))*(-l.z1)^b + f.w(1-l.p,l.g,l.d)*(-l.z2)^b)^(1/b);
	l.ret   <- dnorm(l.ce, l.cehat, l.sigma*(l.z1-l.z2))^l.sel;
	pre.ret <- rbind(g.ret,l.ret);
	for (i in 1:j){
		ret[i,c] <- prod(pre.ret[,i]);
	};
    };
    ret;
};

ftl <- function(vin,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc){
	lvin <- length(vin);
	v    <- vin[1:(lvin-nc+1)]
	mix  <- vin[(lvin-nc+2):lvin];
	mix  <- c(mix,1-sum(mix));
	m <- matrix(mix,j,nc, byrow=T);
	ret <- sum(log(rowSums(m*fcd(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc))));
	ret;
};

ftg <- function(vin,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc) {
    lvin <- length(vin);
    v    <- vin[1:(lvin-nc+1)]
    mix  <- vin[(lvin-nc+2):lvin];
    mix  <- c(mix,1-sum(mix));
    pre.a   <- matrix(v[1:nc],k.u,nc);
    pre.b   <- matrix(v[(nc+1):(2*nc)],1,nc);
    pre.g.g <- matrix(v[(2*nc+1):(3*nc)],k.v,nc);
    pre.l.g <- matrix(v[(3*nc+1):(4*nc)],1,nc);
    pre.g.d <- matrix(v[(4*nc+1):(5*nc)],k.w,nc);
    pre.l.d <- matrix(v[(5*nc+1):(6*nc)],1,nc);
    g.sigma  <- matrix(v[(6*nc+1):(6*nc+j)], g.n, j, byrow=T);
    l.sigma  <- matrix(v[(6*nc+j+1):(6*nc+2*j)], l.n, j, byrow=T);
    g.sigmas <- g.sigma*(g.z1-g.z2);
    l.sigmas <- l.sigma*(l.z1-l.z2);
    gnj <- g.n*j;
    lnj <- l.n*j;
    g.sel.u    <- matrix(as.vector(g.sel)   ,gnj,k.u);
    g.z1.u     <- matrix(as.vector(g.z1)    ,gnj,k.u);
    g.z2.u     <- matrix(as.vector(g.z2)    ,gnj,k.u);
    g.p.u      <- matrix(as.vector(g.p)     ,gnj,k.u);
    g.ce.u     <- matrix(as.vector(g.ce)    ,gnj,k.u);
    g.sigmas.u <- matrix(as.vector(g.sigmas),gnj,k.u);
    l.sel.u    <- matrix(as.vector(l.sel)   ,lnj,k.u);
    l.z1.u     <- matrix(as.vector(l.z1)    ,lnj,k.u);
    l.z2.u     <- matrix(as.vector(l.z2)    ,lnj,k.u);
    l.p.u      <- matrix(as.vector(l.p)     ,lnj,k.u);
    l.ce.u     <- matrix(as.vector(l.ce)    ,lnj,k.u);
    l.sigmas.u <- matrix(as.vector(l.sigmas),lnj,k.u);
    g.sel.v    <- matrix(as.vector(g.sel)   ,gnj,k.v);
    g.z1.v     <- matrix(as.vector(g.z1)    ,gnj,k.v);
    g.z2.v     <- matrix(as.vector(g.z2)    ,gnj,k.v);
    g.p.v      <- matrix(as.vector(g.p)     ,gnj,k.v);
    g.ce.v     <- matrix(as.vector(g.ce)    ,gnj,k.v);
    g.sigmas.v <- matrix(as.vector(g.sigmas),gnj,k.v);
    l.sel.v    <- matrix(as.vector(l.sel)   ,lnj,k.v);
    l.z1.v     <- matrix(as.vector(l.z1)    ,lnj,k.v);
    l.z2.v     <- matrix(as.vector(l.z2)    ,lnj,k.v);
    l.p.v      <- matrix(as.vector(l.p)     ,lnj,k.v);
    l.ce.v     <- matrix(as.vector(l.ce)    ,lnj,k.v);
    l.sigmas.v <- matrix(as.vector(l.sigmas),lnj,k.v);
    g.sel.w    <- matrix(as.vector(g.sel)   ,gnj,k.w);
    g.z1.w     <- matrix(as.vector(g.z1)    ,gnj,k.w);
    g.z2.w     <- matrix(as.vector(g.z2)    ,gnj,k.w);
    g.p.w      <- matrix(as.vector(g.p)     ,gnj,k.w);
    g.ce.w     <- matrix(as.vector(g.ce)    ,gnj,k.w);
    g.sigmas.w <- matrix(as.vector(g.sigmas),gnj,k.w);
    l.sel.w    <- matrix(as.vector(l.sel)   ,lnj,k.w);
    l.z1.w     <- matrix(as.vector(l.z1)    ,lnj,k.w);
    l.z2.w     <- matrix(as.vector(l.z2)    ,lnj,k.w);
    l.p.w      <- matrix(as.vector(l.p)     ,lnj,k.w);
    l.ce.w     <- matrix(as.vector(l.ce)    ,lnj,k.w);
    l.sigmas.w <- matrix(as.vector(l.sigmas),lnj,k.w);
    ret.p <- fcd(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc)
    pre.ret.a      <- rep(0,nc*k.u);
    pre.ret.b      <- rep(0,nc*k.u);
    pre.ret.g.g    <- rep(0,nc*k.v);
    pre.ret.l.g    <- rep(0,nc*k.v);
    pre.ret.g.d    <- rep(0,nc*k.w);
    pre.ret.l.d    <- rep(0,nc*k.w);
    pre.ret.sigma  <- rep(0,2*j);
    pre.ret.mix    <- rep(0,nc);
    for (c in 1:nc){
	a   <- matrix(pre.a[,c],g.n,j);
	b   <- matrix(pre.b[,c],l.n,j);
	g.g <- matrix(pre.g.g[,c],g.n,j);
	l.g <- matrix(pre.l.g[,c],l.n,j);
	g.d <- matrix(pre.g.d[,c],g.n,j);
	l.d <- matrix(pre.l.d[,c],l.n,j);
	a.u        <- matrix(pre.a[,c],gnj,k.u);
	b.u        <- matrix(pre.b[,c],lnj,k.u);
	a.v        <- matrix(pre.a[,c],gnj,k.v);
	b.v        <- matrix(pre.b[,c],lnj,k.v);
	a.w        <- matrix(pre.a[,c],gnj,k.w);
	b.w        <- matrix(pre.b[,c],lnj,k.w);
	g.g.u        <- matrix(g.g,gnj,k.u);
	l.g.u        <- matrix(l.g,lnj,k.u);
	g.d.u        <- matrix(g.d,gnj,k.u);
	l.d.u        <- matrix(l.d,lnj,k.u);
	g.g.v        <- matrix(g.g,gnj,k.v);
	l.g.v        <- matrix(l.g,lnj,k.v);
	g.d.v        <- matrix(g.d,gnj,k.v);
	l.d.v        <- matrix(l.d,lnj,k.v);
	g.g.w        <- matrix(g.g,gnj,k.w);
	l.g.w        <- matrix(l.g,lnj,k.w);
	g.d.w        <- matrix(g.d,gnj,k.w);
	l.d.w        <- matrix(l.d,lnj,k.w);
	g.cehat <- (f.w(g.p,g.g,g.d)*g.z1^a + (1-f.w(g.p,g.g,g.d))*g.z2^a)^(1/a);
	l.cehat <- -((1-f.w(1-l.p,l.g,l.d))*(-l.z1)^b + f.w(1-l.p,l.g,l.d)*(-l.z2)^b)^(1/b);
	g.ret   <- dnorm(g.ce, g.cehat, g.sigmas)^g.sel;
	l.ret   <- dnorm(l.ce, l.cehat, l.sigmas)^l.sel;
	g.ret.p.sm         <- matrix(ret.p[,c],g.n,j,byrow=T);
	l.ret.p.sm         <- matrix(ret.p[,c],l.n,j,byrow=T);
	g.ret.p.sum	   <- matrix(1/rowSums(matrix(mix,j,nc,byrow=T)*ret.p),g.n,j,byrow=T);
	l.ret.p.sum        <- matrix(1/rowSums(matrix(mix,j,nc,byrow=T)*ret.p),l.n,j,byrow=T);
	g.ret.p.sum.u      <- matrix(as.vector(g.ret.p.sum) ,gnj,k.u);
	l.ret.p.sum.u      <- matrix(as.vector(l.ret.p.sum) ,lnj,k.u);
	g.ret.p.sum.v      <- matrix(as.vector(g.ret.p.sum) ,gnj,k.v);
	l.ret.p.sum.v      <- matrix(as.vector(l.ret.p.sum) ,lnj,k.v);
	g.ret.p.sum.w      <- matrix(as.vector(g.ret.p.sum) ,gnj,k.w);
	l.ret.p.sum.w      <- matrix(as.vector(l.ret.p.sum) ,lnj,k.w);
	g.ret.p.m          <- matrix(ret.p[,c],g.n,j,byrow=T)/g.ret;
	l.ret.p.m          <- matrix(ret.p[,c],l.n,j,byrow=T)/l.ret;
	g.ret.p.m.u        <- matrix(as.vector(g.ret.p.m) ,gnj,k.u);
	l.ret.p.m.u        <- matrix(as.vector(l.ret.p.m) ,lnj,k.u);
	g.ret.p.m.v        <- matrix(as.vector(g.ret.p.m) ,gnj,k.v);
	l.ret.p.m.v        <- matrix(as.vector(l.ret.p.m) ,lnj,k.v);
	g.ret.p.m.w        <- matrix(as.vector(g.ret.p.m) ,gnj,k.w);
	l.ret.p.m.w        <- matrix(as.vector(l.ret.p.m) ,lnj,k.w);
	a.z.grad <- g.sel.u*g.ret.p.m.u*mix[c]*g.ret.p.sum.u* (-1/2*g.un*2^(1/2)*exp(-1/2*(g.ce.u^2-2*((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(1/a.u)*g.ce.u+((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(2/a.u))/g.sigmas.u^2)*(((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(1/a.u)*g.ce.u*log((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))*g.d.u*g.p.u^g.g.u*g.z1.u^a.u+((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(1/a.u)*g.ce.u*log((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))*(1-g.p.u)^g.g.u*g.z2.u^a.u-((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(1/a.u)*g.ce.u*a.u*g.d.u*g.p.u^g.g.u*g.z1.u^a.u*log(g.z1.u)-((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(1/a.u)*g.ce.u*a.u*(1-g.p.u)^g.g.u*g.z2.u^a.u*log(g.z2.u)-((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(2/a.u)*log((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))*g.d.u*g.p.u^g.g.u*g.z1.u^a.u-((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(2/a.u)*log((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))*(1-g.p.u)^g.g.u*g.z2.u^a.u+((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(2/a.u)*a.u*g.d.u*g.p.u^g.g.u*g.z1.u^a.u*log(g.z1.u)+((g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u)/(g.d.u*g.p.u^g.g.u+(1-g.p.u)^g.g.u))^(2/a.u)*a.u*(1-g.p.u)^g.g.u*g.z2.u^a.u*log(g.z2.u))/g.sigmas.u^3/a.u^2/pi^(1/2)/(g.d.u*g.p.u^g.g.u*g.z1.u^a.u+g.z2.u^a.u*(1-g.p.u)^g.g.u));
	g.g.z.grad <- g.sel.v*g.ret.p.m.v*mix[c]*g.ret.p.sum.v* (1/2*g.vn*g.d.v*2^(1/2)*(1-g.p.v)^g.g.v*g.p.v^g.g.v*exp(-1/2*(g.ce.v^2-2*((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(1/a.v)*g.ce.v+((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(2/a.v))/g.sigmas.v^2)*(((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(1/a.v)*g.ce.v*log(g.p.v)*g.z1.v^a.v-((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(1/a.v)*g.ce.v*g.z1.v^a.v*log(1-g.p.v)-((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(1/a.v)*g.ce.v*g.z2.v^a.v*log(g.p.v)+((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(1/a.v)*g.ce.v*g.z2.v^a.v*log(1-g.p.v)-((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(2/a.v)*log(g.p.v)*g.z1.v^a.v+((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(2/a.v)*g.z1.v^a.v*log(1-g.p.v)+((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(2/a.v)*g.z2.v^a.v*log(g.p.v)-((g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v))^(2/a.v)*g.z2.v^a.v*log(1-g.p.v))/(g.d.v*g.p.v^g.g.v*g.z1.v^a.v+g.z2.v^a.v*(1-g.p.v)^g.g.v)/a.v/g.sigmas.v^3/(g.d.v*g.p.v^g.g.v+(1-g.p.v)^g.g.v)/pi^(1/2));
	g.d.z.grad <- g.sel.w*g.ret.p.m.w*mix[c]*g.ret.p.sum.w* (1/2*g.wn*2^(1/2)*(1-g.p.w)^g.g.w*g.p.w^g.g.w*exp(-1/2*(g.ce.w^2-2*((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(1/a.w)*g.ce.w+((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(2/a.w))/g.sigmas.w^2)*(((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(1/a.w)*g.ce.w*g.z1.w^a.w-((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(1/a.w)*g.ce.w*g.z2.w^a.w-((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(2/a.w)*g.z1.w^a.w+((g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w))^(2/a.w)*g.z2.w^a.w)/(g.d.w*g.p.w^g.g.w*g.z1.w^a.w+g.z2.w^a.w*(1-g.p.w)^g.g.w)/a.w/g.sigmas.w^3/(g.d.w*g.p.w^g.g.w+(1-g.p.w)^g.g.w)/pi^(1/2));
	g.sigmas.grad <- (g.z1-g.z2)*g.sel*g.ret.p.m*mix[c]*g.ret.p.sum* (-1/2*2^(1/2)*exp(-1/2*(g.ce^2-2*((g.d*g.p^g.g*g.z1^a+g.z2^a*(1-g.p)^g.g)/(g.d*g.p^g.g+(1-g.p)^g.g))^(1/a)*g.ce+((g.d*g.p^g.g*g.z1^a+g.z2^a*(1-g.p)^g.g)/(g.d*g.p^g.g+(1-g.p)^g.g))^(2/a))/g.sigmas^2)*(g.sigmas^2-g.ce^2+2*((g.d*g.p^g.g*g.z1^a+g.z2^a*(1-g.p)^g.g)/(g.d*g.p^g.g+(1-g.p)^g.g))^(1/a)*g.ce-((g.d*g.p^g.g*g.z1^a+g.z2^a*(1-g.p)^g.g)/(g.d*g.p^g.g+(1-g.p)^g.g))^(2/a))/pi^(1/2)/g.sigmas^4);
	b.z.grad <- l.sel.u*l.ret.p.m.u*mix[c]*l.ret.p.sum.u* (1/2*l.un*2^(1/2)*exp(-1/2*(l.ce.u^2+2*(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(1/b.u)*l.ce.u+(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(2/b.u))/l.sigmas.u^2)*((((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(1/b.u)*l.ce.u*log(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))*(l.p.u)^l.g.u*(-l.z1.u)^b.u+(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(1/b.u)*l.ce.u*log(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))*l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u-(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(1/b.u)*l.ce.u*b.u*(l.p.u)^l.g.u*(-l.z1.u)^b.u*log(-l.z1.u)-(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(1/b.u)*l.ce.u*b.u*l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u*log(-l.z2.u)+(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(2/b.u)*log(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))*(l.p.u)^l.g.u*(-l.z1.u)^b.u+(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(2/b.u)*log(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))*l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u-(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(2/b.u)*b.u*(l.p.u)^l.g.u*(-l.z1.u)^b.u*log(-l.z1.u)-(((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/(l.d.u*(1-l.p.u)^l.g.u+(l.p.u)^l.g.u))^(2/b.u)*b.u*l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u*log(-l.z2.u))/l.sigmas.u^3/((-l.z1.u)^b.u*(l.p.u)^l.g.u+l.d.u*(1-l.p.u)^l.g.u*(-l.z2.u)^b.u)/b.u^2/pi^(1/2));
	l.g.z.grad <- l.sel.v*l.ret.p.m.v*mix[c]*l.ret.p.sum.v* (-1/2*l.vn*l.d.v*2^(1/2)*(l.p.v)^l.g.v*(1-l.p.v)^l.g.v*exp(-1/2*(l.ce.v^2+2*(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(1/b.v)*l.ce.v+(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(2/b.v))/l.sigmas.v^2)*(-(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(1/b.v)*l.ce.v*(-l.z1.v)^b.v*log((1-l.p.v))+(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(1/b.v)*l.ce.v*(-l.z1.v)^b.v*log(l.p.v)+(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(1/b.v)*l.ce.v*log((1-l.p.v))*(-l.z2.v)^b.v-(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(1/b.v)*l.ce.v*(-l.z2.v)^b.v*log(l.p.v)-(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(2/b.v)*(-l.z1.v)^b.v*log((1-l.p.v))+(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(2/b.v)*(-l.z1.v)^b.v*log(l.p.v)+(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(2/b.v)*log((1-l.p.v))*(-l.z2.v)^b.v-(((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v))^(2/b.v)*(-l.z2.v)^b.v*log(l.p.v))/b.v/l.sigmas.v^3/pi^(1/2)/((-l.z1.v)^b.v*(l.p.v)^l.g.v+l.d.v*(1-l.p.v)^l.g.v*(-l.z2.v)^b.v)/(l.d.v*(1-l.p.v)^l.g.v+(l.p.v)^l.g.v));
	l.d.z.grad <- l.sel.w*l.ret.p.m.w*mix[c]*l.ret.p.sum.w* (1/2*l.wn*2^(1/2)*(l.p.w)^l.g.w*(1-l.p.w)^l.g.w*exp(-1/2*(l.ce.w^2+2*(((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(1/b.w)*l.ce.w+(((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(2/b.w))/l.sigmas.w^2)*((((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(1/b.w)*l.ce.w*(-l.z1.w)^b.w-(((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(1/b.w)*l.ce.w*(-l.z2.w)^b.w+(((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(2/b.w)*(-l.z1.w)^b.w-(((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w))^(2/b.w)*(-l.z2.w)^b.w)/b.w/l.sigmas.w^3/pi^(1/2)/((-l.z1.w)^b.w*(l.p.w)^l.g.w+l.d.w*(1-l.p.w)^l.g.w*(-l.z2.w)^b.w)/(l.d.w*(1-l.p.w)^l.g.w+(l.p.w)^l.g.w));
	l.sigmas.grad <- (l.z1-l.z2)*l.sel*l.ret.p.m*mix[c]*l.ret.p.sum* (1/2*2^(1/2)*exp(-1/2*(l.ce^2+2*(((-l.z1)^b*(l.p)^l.g+l.d*(1-l.p)^l.g*(-l.z2)^b)/(l.d*(1-l.p)^l.g+(l.p)^l.g))^(1/b)*l.ce+(((-l.z1)^b*(l.p)^l.g+l.d*(1-l.p)^l.g*(-l.z2)^b)/(l.d*(1-l.p)^l.g+(l.p)^l.g))^(2/b))/l.sigmas^2)*(-l.sigmas^2+l.ce^2+2*(((-l.z1)^b*(l.p)^l.g+l.d*(1-l.p)^l.g*(-l.z2)^b)/(l.d*(1-l.p)^l.g+(l.p)^l.g))^(1/b)*l.ce+(((-l.z1)^b*(l.p)^l.g+l.d*(1-l.p)^l.g*(-l.z2)^b)/(l.d*(1-l.p)^l.g+(l.p)^l.g))^(2/b))/pi^(1/2)/l.sigmas^4);
	mix.grad <- ret.p[,c]*g.ret.p.sum[1,];
	pre.ret.a[((c-1)*k.u+1):(c*k.u)]   <- colSums(a.z.grad);
	pre.ret.b[((c-1)*k.u+1):(c*k.u)]   <- colSums(b.z.grad);
	pre.ret.g.g[((c-1)*k.v+1):(c*k.v)] <- colSums(g.g.z.grad);
	pre.ret.l.g[((c-1)*k.v+1):(c*k.v)] <- colSums(l.g.z.grad);
	pre.ret.g.d[((c-1)*k.w+1):(c*k.w)] <- colSums(g.d.z.grad);
	pre.ret.l.d[((c-1)*k.w+1):(c*k.w)] <- colSums(l.d.z.grad);
	pre.ret.sigma                      <- pre.ret.sigma + c(colSums(g.sigmas.grad),colSums(l.sigmas.grad));
	pre.ret.mix[c]                     <- sum(mix.grad);
    }; 
    pre.ret.mix <- pre.ret.mix-pre.ret.mix[nc];
    grad.ret    <- c(pre.ret.a, pre.ret.b, pre.ret.g.g, pre.ret.l.g, pre.ret.g.d, pre.ret.l.d, pre.ret.sigma,pre.ret.mix[1:(nc-1)]);
    grad.ret;
};

fml <- function(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau){
	sum(tau*log(fcd(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc)));
};

fmp <- function(v,weight,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn) {
    a   <- v[1];
    b   <- v[2];
    g.g <- v[3];
    l.g <- v[4];
    g.d <- v[5];
    l.d <- v[6];
    g.sigma  <- matrix(v[7:(6+j)], g.n, j, byrow=T);
    l.sigma  <- matrix(v[(7+j):(2*(k.u+k.v+k.w+j))], l.n, j, byrow=T);
    g.sigmas <- g.sigma*(g.z1-g.z2);
    l.sigmas <- l.sigma*(l.z1-l.z2);
    gnj <- g.n*j;
    lnj <- l.n*j;
    g.cehat  <- (f.w(g.p,g.g,g.d)*g.z1^a + (1-f.w(g.p,g.g,g.d))*g.z2^a)^(1/a);
    l.cehat  <- -((1-f.w(1-l.p,l.g,l.d))*(-l.z1)^b + f.w(1-l.p,l.g,l.d)*(-l.z2)^b)^(1/b);
    a.u        <- matrix(a       ,gnj,k.u);
    b.u        <- matrix(b       ,lnj,k.u);
    g.g.u      <- matrix(g.g     ,gnj,k.u);
    l.g.u      <- matrix(l.g     ,lnj,k.u);
    g.d.u      <- matrix(g.d     ,gnj,k.u);
    l.d.u      <- matrix(l.d     ,lnj,k.u);
    g.sel.u    <- matrix(as.vector(g.sel)   ,gnj,k.u);
    g.z1.u     <- matrix(as.vector(g.z1)    ,gnj,k.u);
    g.z2.u     <- matrix(as.vector(g.z2)    ,gnj,k.u);
    g.p.u      <- matrix(as.vector(g.p)     ,gnj,k.u);
    g.ce.u     <- matrix(as.vector(g.ce)    ,gnj,k.u);
    g.sigmas.u <- matrix(as.vector(g.sigmas),gnj,k.u);
    l.sel.u    <- matrix(as.vector(l.sel)   ,lnj,k.u);
    l.z1.u     <- matrix(as.vector(l.z1)    ,lnj,k.u);
    l.z2.u     <- matrix(as.vector(l.z2)    ,lnj,k.u);
    l.p.u      <- matrix(as.vector(l.p)     ,lnj,k.u);
    l.ce.u     <- matrix(as.vector(l.ce)    ,lnj,k.u);
    l.sigmas.u <- matrix(as.vector(l.sigmas),lnj,k.u);
    a.v        <- matrix(a       ,gnj,k.v);
    b.v        <- matrix(b       ,lnj,k.v);
    g.g.v      <- matrix(g.g     ,gnj,k.v);
    l.g.v      <- matrix(l.g     ,lnj,k.v);
    g.d.v      <- matrix(g.d     ,gnj,k.v);
    l.d.v      <- matrix(l.d     ,lnj,k.v);
    g.sel.v    <- matrix(as.vector(g.sel)   ,gnj,k.v);
    g.z1.v     <- matrix(as.vector(g.z1)    ,gnj,k.v);
    g.z2.v     <- matrix(as.vector(g.z2)    ,gnj,k.v);
    g.p.v      <- matrix(as.vector(g.p)     ,gnj,k.v);
    g.ce.v     <- matrix(as.vector(g.ce)    ,gnj,k.v);
    g.sigmas.v <- matrix(as.vector(g.sigmas),gnj,k.v);
    l.sel.v    <- matrix(as.vector(l.sel)   ,lnj,k.v);
    l.z1.v     <- matrix(as.vector(l.z1)    ,lnj,k.v);
    l.z2.v     <- matrix(as.vector(l.z2)    ,lnj,k.v);
    l.p.v      <- matrix(as.vector(l.p)     ,lnj,k.v);
    l.ce.v     <- matrix(as.vector(l.ce)    ,lnj,k.v);
    l.sigmas.v <- matrix(as.vector(l.sigmas),lnj,k.v);
    a.w        <- matrix(a       ,gnj,k.w);
    b.w        <- matrix(b       ,lnj,k.w);
    g.g.w      <- matrix(g.g     ,gnj,k.w);
    l.g.w      <- matrix(l.g     ,lnj,k.w);
    g.d.w      <- matrix(g.d     ,gnj,k.w);
    l.d.w      <- matrix(l.d     ,lnj,k.w);
    g.sel.w    <- matrix(as.vector(g.sel)   ,gnj,k.w);
    g.z1.w     <- matrix(as.vector(g.z1)    ,gnj,k.w);
    g.z2.w     <- matrix(as.vector(g.z2)    ,gnj,k.w);
    g.p.w      <- matrix(as.vector(g.p)     ,gnj,k.w);
    g.ce.w     <- matrix(as.vector(g.ce)    ,gnj,k.w);
    g.sigmas.w <- matrix(as.vector(g.sigmas),gnj,k.w);
    l.sel.w    <- matrix(as.vector(l.sel)   ,lnj,k.w);
    l.z1.w     <- matrix(as.vector(l.z1)    ,lnj,k.w);
    l.z2.w     <- matrix(as.vector(l.z2)    ,lnj,k.w);
    l.p.w      <- matrix(as.vector(l.p)     ,lnj,k.w);
    l.ce.w     <- matrix(as.vector(l.ce)    ,lnj,k.w);
    l.sigmas.w <- matrix(as.vector(l.sigmas),lnj,k.w);
    g.weights          <- matrix(weight,g.n,j,byrow=T);
    l.weights	       <- matrix(weight,l.n,j,byrow=T);
    g.weights.u        <- matrix(as.vector(g.weights) ,gnj,k.u);
    l.weights.u        <- matrix(as.vector(l.weights) ,lnj,k.u);
    g.weights.v        <- matrix(as.vector(g.weights) ,gnj,k.v);
    l.weights.v        <- matrix(as.vector(l.weights) ,lnj,k.v);
    g.weights.w        <- matrix(as.vector(g.weights) ,gnj,k.w);
    l.weights.w        <- matrix(as.vector(l.weights) ,lnj,k.w);
    der.a   <- g.un;
    der.b   <- l.un;
    der.g.g <- g.vn;
    der.l.g <- l.vn;
    der.g.d <- g.wn;
    der.l.d <- l.wn;
    a.z.grad     <- g.sel.u*g.weights.u*(-(((((1 - g.p.u)^g.g.u*g.z2.u^a.u + g.p.u^g.g.u*g.z1.u^a.u*g.d.u)/((1 - g.p.u)^g.g.u + g.p.u^g.g.u*g.d.u))^(-1 + a.u^(-1))*(-g.ce.u + (((1 - g.p.u)^g.g.u*g.z2.u^a.u + g.p.u^g.g.u*g.z1.u^a.u*g.d.u)/((1 - g.p.u)^g.g.u + g.p.u^g.g.u*g.d.u))^a.u^(-1))*(a.u*(g.p.u^g.g.u*g.z1.u^a.u*g.d.u*log(g.z1.u) + (1 - g.p.u)^g.g.u*g.z2.u^a.u*log(g.z2.u)) - ((1 - g.p.u)^g.g.u*g.z2.u^a.u + g.p.u^g.g.u*g.z1.u^a.u*g.d.u)*log(((1 - g.p.u)^g.g.u*g.z2.u^a.u + g.p.u^g.g.u*g.z1.u^a.u*g.d.u)/((1 - g.p.u)^g.g.u + g.p.u^g.g.u*g.d.u)))*der.a)/(g.sigmas.u^2*a.u^2*((1 - g.p.u)^g.g.u + g.p.u^g.g.u*g.d.u))));
    g.g.z.grad   <- g.sel.v*g.weights.v*(((-((-1 + g.p.v)*g.p.v))^g.g.v*(g.z1.v^a.v - g.z2.v^a.v)*g.d.v*(((1 - g.p.v)^g.g.v*g.z2.v^a.v + g.p.v^g.g.v*g.z1.v^a.v*g.d.v)/((1 - g.p.v)^g.g.v + g.p.v^g.g.v*g.d.v))^(-1 + a.v^(-1))*(-g.ce.v + (((1 - g.p.v)^g.g.v*g.z2.v^a.v + g.p.v^g.g.v*g.z1.v^a.v*g.d.v)/((1 - g.p.v)^g.g.v + g.p.v^g.g.v*g.d.v))^a.v^(-1))*(log(1 - g.p.v) - log(g.p.v))*der.g.g)/(g.sigmas.v^2*a.v*((1 - g.p.v)^g.g.v + g.p.v^g.g.v*g.d.v)^2));
    g.d.z.grad   <- g.sel.w*g.weights.w*(-(((-((-1 + g.p.w)*g.p.w))^g.g.w*(g.z1.w^a.w - g.z2.w^a.w)*(((1 - g.p.w)^g.g.w*g.z2.w^a.w + g.p.w^g.g.w*g.z1.w^a.w*g.d.w)/((1 - g.p.w)^g.g.w + g.p.w^g.g.w*g.d.w))^(-1 + a.w^(-1))*(-g.ce.w + (((1 - g.p.w)^g.g.w*g.z2.w^a.w + g.p.w^g.g.w*g.z1.w^a.w*g.d.w)/((1 - g.p.w)^g.g.w + g.p.w^g.g.w*g.d.w))^a.w^(-1))*der.g.d)/(g.sigmas.w^2*a.w*((1 - g.p.w)^g.g.w + g.p.w^g.g.w*g.d.w)^2)));
    g.grad.sigma <- g.sel*g.weights*(-1/g.sigma+(g.ce-g.cehat)^2/(g.sigma^3*(g.z1-g.z2)^2));
    b.z.grad     <- l.sel.u*l.weights.u*(-((((l.p.u^l.g.u*(-l.z1.u)^b.u + (1 - l.p.u)^l.g.u*(-l.z2.u)^b.u*l.d.u)/(l.p.u^l.g.u + (1 - l.p.u)^l.g.u*l.d.u))^(-1 + b.u^(-1))*(l.ce.u + ((l.p.u^l.g.u*(-l.z1.u)^b.u + (1 - l.p.u)^l.g.u*(-l.z2.u)^b.u*l.d.u)/(l.p.u^l.g.u + (1 - l.p.u)^l.g.u*l.d.u))^b.u^(-1))*(b.u*(l.p.u^l.g.u*(-l.z1.u)^b.u*log(-l.z1.u) + (1 - l.p.u)^l.g.u*(-l.z2.u)^b.u*l.d.u*log(-l.z2.u)) - (l.p.u^l.g.u*(-l.z1.u)^b.u + (1 - l.p.u)^l.g.u*(-l.z2.u)^b.u*l.d.u)*log((l.p.u^l.g.u*(-l.z1.u)^b.u + (1 - l.p.u)^l.g.u*(-l.z2.u)^b.u*l.d.u)/(l.p.u^l.g.u + (1 - l.p.u)^l.g.u*l.d.u)))*der.b)/(l.sigmas.u^2*b.u^2*(l.p.u^l.g.u + (1 - l.p.u)^l.g.u*l.d.u))));
    l.g.z.grad   <- l.sel.v*l.weights.v*(((-((-1 + l.p.v)*l.p.v))^l.g.v*((-l.z1.v)^b.v - (-l.z2.v)^b.v)*l.d.v*((l.p.v^l.g.v*(-l.z1.v)^b.v + (1 - l.p.v)^l.g.v*(-l.z2.v)^b.v*l.d.v)/(l.p.v^l.g.v + (1 - l.p.v)^l.g.v*l.d.v))^(-1 + b.v^(-1))*(l.ce.v + ((l.p.v^l.g.v*(-l.z1.v)^b.v + (1 - l.p.v)^l.g.v*(-l.z2.v)^b.v*l.d.v)/(l.p.v^l.g.v + (1 - l.p.v)^l.g.v*l.d.v))^b.v^(-1))*(log(1 - l.p.v) - log(l.p.v))*der.l.g)/(l.sigmas.v^2*b.v*(l.p.v^l.g.v + (1 - l.p.v)^l.g.v*l.d.v)^2));
    l.d.z.grad   <- l.sel.w*l.weights.w*(((-((-1 + l.p.w)*l.p.w))^l.g.w*((-l.z1.w)^b.w - (-l.z2.w)^b.w)*((l.p.w^l.g.w*(-l.z1.w)^b.w + (1 - l.p.w)^l.g.w*(-l.z2.w)^b.w*l.d.w)/(l.p.w^l.g.w + (1 - l.p.w)^l.g.w*l.d.w))^(-1 + b.w^(-1))*(l.ce.w + ((l.p.w^l.g.w*(-l.z1.w)^b.w + (1 - l.p.w)^l.g.w*(-l.z2.w)^b.w*l.d.w)/(l.p.w^l.g.w + (1 - l.p.w)^l.g.w*l.d.w))^b.w^(-1))*der.l.d)/(l.sigmas.w^2*b.w*(l.p.w^l.g.w + (1 - l.p.w)^l.g.w*l.d.w)^2));
    l.grad.sigma <- l.sel*l.weights*(-1/l.sigma+(l.ce-l.cehat)^2/(l.sigma^3*(l.z1-l.z2)^2));
    ret<-list(a.grad=colSums(a.z.grad),b.grad=colSums(b.z.grad),g.g.grad=colSums(g.g.z.grad),l.g.grad=colSums(l.g.z.grad),g.d.grad=colSums(g.d.z.grad),l.d.grad=colSums(l.d.z.grad),sigma.grad=c(colSums(g.grad.sigma),colSums(l.grad.sigma)));

    ret;
};

fmg <- function(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau){
	a.zeta   <- matrix(v[1:(nc*k.u)],k.u,nc);
	b.zeta   <- matrix(v[(nc*k.u+1):(2*nc*k.u)],k.u,nc);
	g.g.zeta <- matrix(v[(2*nc*k.u+1):(2*nc*k.u+nc)],k.v,nc);
	l.g.zeta <- matrix(v[(2*nc*k.u+nc*k.v+1):(2*nc*(1+k.v))],k.v,nc);
	g.d.zeta <- matrix(v[(2*nc*(k.u+k.v)+1):(2*nc*(1+k.v)+nc*k.w)],k.w,nc);
	l.d.zeta <- matrix(v[(2*nc*(k.u+k.v)+nc*k.w+1):(2*nc*(1+k.v+k.w))],k.w,nc);
	sigma    <- v[(2*nc*(k.u+k.v+k.w)+1):(2*(nc*(1+k.v+k.w)+j))];
	pre.ret.a      <- rep(0,nc*k.u);
	pre.ret.b      <- rep(0,nc);
	pre.ret.g.g    <- rep(0,nc*k.v);
	pre.ret.l.g    <- rep(0,nc);
	pre.ret.g.d    <- rep(0,nc*k.w);
	pre.ret.l.d    <- rep(0,nc);
	pre.ret.sigma  <- rep(0,2*j);
	for (c in 1:nc){
		vc <- c(a.zeta[,c],b.zeta[,c],g.g.zeta[,c],l.g.zeta[,c],g.d.zeta[,c],l.d.zeta[,c],sigma);
		pre.ret <- fmp(vc,weight=tau[,c],g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn);
		pre.ret.a[((c-1)*k.u+1):(c*k.u)]   <- pre.ret$a.grad;
		pre.ret.b[((c-1)*k.u+1):(c*k.u)]   <- pre.ret$b.grad;
		pre.ret.g.g[((c-1)*k.v+1):(c*k.v)] <- pre.ret$g.g.grad;
		pre.ret.l.g[((c-1)*k.v+1):(c*k.v)] <- pre.ret$l.g.grad;
		pre.ret.g.d[((c-1)*k.w+1):(c*k.w)] <- pre.ret$g.d.grad;
		pre.ret.l.d[((c-1)*k.w+1):(c*k.w)] <- pre.ret$l.d.grad;
		pre.ret.sigma                      <- pre.ret.sigma + pre.ret$sigma.grad
	};
	ret <- c(pre.ret.a, pre.ret.b, pre.ret.g.g, pre.ret.l.g, pre.ret.g.d, pre.ret.l.d, pre.ret.sigma);
	ret;
};

fes <- function(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,mix){
	m <- matrix(mix,j,nc, byrow=T);
	numerator   <- m*fcd(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,1,1,1,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc);
	denominator <- matrix(rowSums(numerator),j,nc);
	numerator/denominator;
};

fss <- function(em.tau,j,nc){
	sem.tau <- matrix(0,j,nc);
	for (i in 1:j){
		sem.tau[i,] <- t(rmultinom(1,1,em.tau[i,]));
	};
	sem.tau;
};

fms <- function(v0,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau){
	mix <- colMeans(tau);
	npar      <- 2*nc*3;
	nsigma    <- 2*j;
	ui1	  <- matrix(0,nsigma,npar);
	ui2	  <- matrix(0,nsigma,nsigma);
	diag(ui2) <- rep(1,nsigma);
	ui        <- cbind(ui1,ui2);
	ci	  <- rep(0,nsigma);
	fstart <- fml(v0,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau);
	v1 <- v0;
	while (is.na(fstart) | fstart==Inf | fstart == -Inf) {
		v1 <- c(runif(npar,.9,1.1), runif(nsigma,0.6,1.4))*v0[1:(npar+nsigma)];
		vsigma <- v1[(npar+1):(npar+nsigma)];
		vsigma[vsigma>1] <- 1-1e-4;
		v1[(npar+1):(npar+nsigma)] <- vsigma;
		v1 <- c(v1,(runif(nc-1.2,1.8)*v0[(npar+nsigma+1):length(v0)]));
		fstart <- fml(v1,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,1,1,1,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau);
	};
	v   <- constrOptim(v1,fml,fmg,method="BFGS",control=list(maxit=100,trace=0,REPORT=1, fnscale=-1), ui=ui, ci=ci, g.p=g.p,l.p=l.p,g.z1=g.z1,l.z1=l.z1,g.z2=g.z2,l.z2=l.z2,g.ce=g.ce,l.ce=l.ce,g.sel=g.sel,l.sel=l.sel,g.n=g.n,l.n=l.n,j=j,k.u=k.u,k.v=k.v,k.w=k.w,g.un=g.un,g.vn=g.vn,g.wn=g.wn,l.un=l.un,l.vn=l.vn,l.wn=l.wn,nc=nc,tau=tau);
	list(mix=mix, v=v$par);
};

fdm <- function(in.dat,unames,vnames,wnames,n,j,up,dn,e,epsilon) {
    tdat <- table(in.dat[,"id"]);
    j    <- length(tdat);
    pid  <- as.numeric(names(tdat));
    tdat.l <- table(in.dat[,"lottery"]);
    ac.n   <- length(tdat.l);
    lid    <- as.numeric(names(tdat.l));
    id   <- matrix(0,n,j);
    z1   <- matrix(0,n,j);
    z2   <- matrix(0,n,j);
    p    <- matrix(0,n,j);
    ce   <- matrix(0,n,j);
    sel  <- matrix(0,n,j);
    un      <- matrix(1,n*j,1);
    vn      <- matrix(1,n*j,1);
    wn      <- matrix(1,n*j,1);
    for (i in 1:j){
        sel[,i] <- c(rep(1,tdat[i]),rep(0,n-tdat[i])); 
        idx      <- in.dat[,"id"]==pid[i];
	id[,i]   <- rep(pid[i],n);
        z1[,i]   <- c(in.dat[idx,"z1"],rep(up,n-tdat[i]))+epsilon;
        z2[,i]   <- c(in.dat[idx,"z2"],rep(dn,n-tdat[i]))+epsilon;
        p[,i]    <- c(in.dat[idx,"p1"],rep(0.5,n-tdat[i]));
        ce[,i]   <- c(in.dat[idx,"ce"],rep(e,n-tdat[i]))+epsilon;
    };
    list(id=id, z1=z1, z2=z2, p=p, un=un, vn=vn, wn=wn, ce=ce, sel=sel, j=j, n=n);
};

fdp <- function(dat,jorig=0) {
    map <- max(c(abs(dat[,"z1"]),abs(dat[,"z2"])))
    dat[,"z1"]<-dat[,"z1"]/map;
    dat[,"z2"]<-dat[,"z2"]/map;
    dat[,"ce"]<-dat[,"ce"]/map;
    g.idx <- dat$z1>0 & dat$z2>=0 & dat$ce>=0;
    l.idx <- dat$z1<=0 & dat$z2<0 & dat$ce<=0;
    g.dat <- dat[g.idx, ];
    l.dat <- dat[l.idx, ];
    g.n <- max(table(g.dat[,"id"]));
    l.n <- max(table(l.dat[,"id"]));
    n   <- max(g.n,l.n);
    k.u   <- 1;
    k.v   <- 1;
    k.w   <- 1;
    g.ret <- fdm(g.dat,unames,vnames,wnames,n,j,1,0.5,0.75,1e-18);
    l.ret <- fdm(l.dat,unames,vnames,wnames,n,j,-0.5,-1,-0.75,-1e-18);
    jdiff <- F;
    {if (g.ret$j != l.ret$j) {jdiff <- T;}};
    {if (g.ret$j < jorig | l.ret$j < jorig) {jdiff <- T;}};
    list(id=g.ret$id, j=g.ret$j, k.u=k.u, k.v=k.v, k.w=k.w, g.p=g.ret$p, g.z1=g.ret$z1, g.z2=g.ret$z2, g.ce=g.ret$ce, g.sel=g.ret$sel, g.un=g.ret$un, g.vn=g.ret$vn, g.wn=g.ret$wn, g.n=g.ret$n, l.p=l.ret$p, l.z1=l.ret$z1, l.z2=l.ret$z2, l.ce=l.ret$ce, l.sel=l.ret$sel, l.un=l.ret$un, l.vn=l.ret$vn, l.wn=l.ret$wn, l.n=l.ret$n, jdiff=jdiff);
};

fan <- function(em.mix, em.vm, sem.mix, sem.vm, iter, b=20){
	{if (iter <= b)
		{q <- cos((iter/b)*acos(b/100));}
	 else
		{q <- (b/100)*sqrt(b/iter);}
	};
	mix <- q*sem.mix+(1-q)*em.mix;
	vm  <- q*sem.vm+(1-q)*em.vm;
	list(mix=mix, vm=vm);
};
#-------------------------------------------------------------------------------
fest <- function(datastring, v0file, nc){
	dat       <- read.table(datastring, header=T, sep=",");

	#Save the original Data-Dimensions
	tdat  <- table(dat[,c("id")]);
	obsid <- as.numeric(names(tdat));
	j     <- length(tdat); 
	n     <- max(tdat);   
	prep  <- fdp(dat);
	id    <- prep$id;
	j     <- prep$j;
	k.u   <- prep$k.u;
	k.v   <- prep$k.v;
	k.w   <- prep$k.w;
	g.z1  <- prep$g.z1;
	g.z2  <- prep$g.z2;
	g.ce  <- prep$g.ce;
	g.p   <- prep$g.p;
	g.sel <- prep$g.sel;
	g.un  <- prep$g.un;
	g.vn  <- prep$g.vn;
	g.wn  <- prep$g.wn;
	g.n   <- prep$g.n;
	l.z1  <- prep$l.z1;
	l.z2  <- prep$l.z2;
	l.ce  <- prep$l.ce;
	l.p   <- prep$l.p;
	l.sel <- prep$l.sel;
	l.un  <- prep$l.un;
	l.vn  <- prep$l.vn;
	l.wn  <- prep$l.wn;
	l.n   <- prep$l.n;
	rm(prep);
	parv <- read.table(v0file, sep=",", header=F)[,1];
	ab     <- c(rep(parv[1:k.u],nc),rep(parv[(k.u+1):(2*k.u)],nc));
	g.g    <- rep(parv[(2*k.u+1):(2*k.u+k.v)],nc)
	l.g    <- rep(parv[(2*k.u+k.v+1):(2*(k.u+k.v))],nc)
	g.d    <- rep(parv[(2*(k.u+k.v)+1):(2*(k.u+k.v)+k.w)],nc)
	l.d    <- rep(parv[(2*(k.u+k.v)+k.w+1):(2*(k.u+k.v+k.w))],nc)
	sigmav <- parv[(2*(k.u+k.v+k.w)+1):length(parv)];
	v    <- c(ab,g.g*runif(length(g.g),.9,1.1),l.g*runif(length(l.g),.9,1.1),g.d*runif(length(g.d),.9,1.1),l.d*runif(length(l.d),.9,1.1),sigmav);
	aa   <- runif(nc,0.5,2);
	mix  <- aa/sum(aa);
	rm(aa);
	np   <- length(v); 
  	iter  <- 1;
	diff  <- 1;
	llold <- 0;
	llnew <- 1;
	print("saem est",quote=F);
  	while (iter < 21 & abs(diff)>1e-3){
  		tau.em  <- fes(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,mix);
  		mstep   <- fms(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau.em);
  		v.em    <- mstep$v;
  		mix.em  <- mstep$mix;
  		rm(mstep);
		tau.sem <- fss(tau.em,j,nc);
		mstep   <- fms(v,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc,tau.sem);
		v.sem   <- mstep$v;
		mix.sem <- mstep$mix;
		rm(mstep);
		ann <- fan(mix.em, v.em, mix.sem, v.sem, iter, b=50);
		mix <- ann$mix;
		v   <- ann$vm;
		rm(ann);
		print(paste("iter",iter), quote=F)
		iter <- iter+1;
		llold <- llnew;
		v0    <- c(v,mix[1:(length(mix)-1)]);
		llnew <- ftl(v0,g.p,l.p,g.z1,l.z1,g.z2,l.z2,g.ce,l.ce,g.sel,l.sel,g.n,l.n,j,k.u,k.v,k.w,g.un,g.vn,g.wn,l.un,l.vn,l.wn,nc);
		diff  <- llnew-llold;
	};
	print("ml est",quote=F)
	ret <- optim(v0,ftl, ftg, method="BFGS", control=list(maxit=100000,trace=0,REPORT=1, fnscale=-1), hessian=F, g.p=g.p,l.p=l.p,g.z1=g.z1,l.z1=l.z1,g.z2=g.z2,l.z2=l.z2,g.ce=g.ce,l.ce=l.ce,g.sel=g.sel,l.sel=l.sel,g.n=g.n,l.n=l.n,j=j,k.u=k.u,k.v=k.v,k.w=k.w,g.un=g.un,g.vn=g.vn,g.wn=g.wn,l.un=l.un,l.vn=l.vn,l.wn=l.wn,nc=nc);
	v        <- ret$par;
	mix            <- v[(length(v)-nc+2):length(v)];
	mix            <- c(mix,1-sum(mix));
	om2      <- matrix(ret$par[1:6], 6, nc, byrow=T);
	dimnames(om2) <- list(c("alpha:","beta:","g.gamma:","l.gamma:","g.delta:","l.delta:"),paste("type",1:nc));
	print("mixture:")
	print(round(mix,3))
	print(round(om2,3))
};

china05<-fest("Desktop/cn05_data.csv","Desktop/cn05_startval_model1cl.csv", 2);
#swiss06<-fest("Desktop/ch06_data.csv","Desktop/ch06_startval_model1cl.csv", 2);
#swiss03<-fest("Desktop/ch03_data.csv","Desktop/ch03_startval_model1cl.csv", 2);
