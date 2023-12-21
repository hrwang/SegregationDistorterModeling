## 1. functions -----

## the order of haplotypes
# AB
# Ab
# aB
# ab

## The order of haplotype combinations
# AB|AB
# AB|Ab
# AB|aB
# AB|ab
# Ab|Ab
# Ab|aB
# Ab|ab
# aB|aB
# aB|ab
# ab|ab



## calculate haplotype frequency from the 10 possible haplotype combinations (genotypes)
genotype_to_allele_freq <- function(g){
  AB = g[1] + 0.5*g[2] + 0.5*g[3] + 0.5*g[4]
  Ab = 0.5*g[2] + g[5] + 0.5*g[6] + 0.5*g[7]
  aB = 0.5*g[3] + 0.5*g[6] + g[8] + 0.5*g[9]
  ab = 0.5*g[4] + 0.5*g[7] + 0.5*g[9] + g[10]
  return(c(AB, Ab, aB, ab))
}

## calculate the pool pollen number given pollen number for each flower (N), inbreeding coefficient (F) and killing efficiency (k)
pool_pollen_number <- function (g, N=4000, indF=0.5, k=1){
  M = N*(1-indF)/indF
  #M - g[4]*M*k/2 - g[6]*M*k/2 - g[7]*M*k/2 - g[9]*M*k - g[10]*M*k
}


## calculate pollen haplotype frequency given the genotype (g), recombination (r) and killing efficiency (k)
pool_allele_freq <- function (g, r=0.2, k=1) {
  
  G2 <- matrix( c(1, 0, 0, 0,
                  0.5, 0.5, 0, 0,
                  0.5, 0, 0.5, 0,
                  0.5*(1-r), 0.5*r, 0.5*r*(1-k), 0.5*(1-r)*(1-k),
                  0, 1, 0, 0,
                  0.5*r, 0.5*(1-r), 0.5*(1-r)*(1-k), 0.5*r*(1-k),
                  0, 0.5, 0, 0.5*(1-k),
                  0, 0, 1, 0,
                  0, 0, 0.5*(1-k), 0.5*(1-k),
                  0, 0, 0, (1-k)
  ), ncol=4, byrow = TRUE )
  
  x = g %*% G2
  x = x/sum(x)
  return(x)
  
}


## the function between pollen number and fertility.
fertility_grain_number <- function (n=1000, N=4000, indF = 0.5, R=10) {
  M = N*(1-indF)/indF
  x = n/(M+N) * R  ## change 1 to other number > 1 to set the redundancy levels
  x[x>1]=1
  return(x)
}


## the iteration function, calculate genotype of next generation given the genotype of current generation (g),
## pool haplotype frequency (f), pollen number from pool (m), a vector of pollen supply for each genotype (n) and pollen number generated for each flower (N).


nextgen_genotype = function (g, f, m, n, N=4000, indF = 0.5, R=10, r, k) {
  
  ## g: genotype matrix of current generation
  ## f: pool haplotype frequency
  ## m: pollen number from pool
  ## n: a vector of pollen supply for each genotype
  ## N: the number of pollen one flower generates
  ## indF: inbreeding coefficient
  ## R: pollen redundancy level
  ## r: recombination between the two loci.
  ## k: killing efficiency.
  
  fAB = f[1]
  fAb = f[2]
  faB = f[3]
  fab = f[4]
  
  P <- matrix( c((N+m*fAB)/n[1], m*fAb/n[1], m*faB/n[1], m*fab/n[1], 0, 0, 0, 0, 0, 0,
                 1/2*(N/2+m*fAB)/n[2], ## 2nd row.
                 1/2*(N/2+m*fAb)/n[2] + 1/2*(N/2+m*fAB)/n[2],
                 1/2*m*faB/n[2],
                 1/2*m*fab/n[2],
                 1/2*(N/2+m*fAb)/n[2],
                 1/2*m*faB/n[2],
                 1/2*m*fab/n[2],
                 0, 0, 0,
                 1/2*(N/2+m*fAB)/n[3], ## 3rd row.
                 1/2*m*fAb/n[3],
                 1/2*(N/2+m*faB)/n[3] + 1/2*(N/2+m*fAB)/n[3],
                 1/2*m*fab/n[3],
                 0,
                 1/2*m*fAb/n[3],
                 0,
                 1/2*(N/2+m*faB)/n[3],
                 1/2*m*fab/n[3], 0,
                 (1-r)*1/2*(N/2*(1-r)+m*fAB)/n[4],  ## 4th row
                 (1-r)*1/2*(N/2*r+m*fAb)/n[4] + r*1/2*(N/2*(1-r)+m*fAB)/n[4],
                 (1-r)*1/2*(N/2*r*(1-k)+m*faB)/n[4] + r*1/2*(N/2*(1-r)+m*fAB)/n[4],
                 (1-r)*1/2*(N/2*(1-r)*(1-k)+m*fab)/n[4] + (1-r)*1/2*(N/2*(1-r)+m*fAB)/n[4],
                 r*1/2*(N/2*r+m*fAb)/n[4],
                 r*1/2*(N/2*r*(1-k)+m*faB)/n[4] + r*1/2*(N/2*r+m*fAb)/n[4],
                 r*1/2*(N/2*(1-r)*(1-k)+m*fab)/n[4] + (1-r)*1/2*(N/2*r+m*fAb)/n[4],
                 r*1/2*(N/2*r*(1-k)+m*faB)/n[4],
                 r*1/2*(N/2*(1-r)*(1-k)+m*fab)/n[4] + (1-r)*1/2*(N/2*r*(1-k)+m*faB)/n[4],
                 (1-r)*1/2*(N/2*(1-r)*(1-k)+m*fab)/n[4],
                 0, ## 5th row
                 m*fAB/n[5], 0, 0, (N+m*fAb)/n[5], m*faB/n[5], m*fab/n[5], 0, 0, 0,
                 r*1/2*(N/2*r+m*fAB)/n[6], ## 6th row
                 r*1/2*(N/2*(1-r)+m*fAb)/n[6] + (1-r)*1/2*(N/2*r+m*fAB)/n[6],
                 r*1/2*(N/2*(1-r)*(1-k)+m*faB)/n[6] + (1-r)*1/2*(N/2*r+m*fAB)/n[6],
                 r*1/2*(N/2*r*(1-k)+m*fab)/n[6] + r*1/2*(N/2*r+m*fAB)/n[6],
                 (1-r)*1/2*(N/2*(1-r)+m*fAb)/n[6],
                 (1-r)*1/2*(N/2*(1-r)*(1-k)+m*faB)/n[6] + (1-r)*1/2*(N/2*(1-r)+m*fAb)/n[6],
                 (1-r)*1/2*(N/2*r*(1-k)+m*fab)/n[6] + r*1/2*(N/2*(1-r)+m*fAb)/n[6],
                 (1-r)*1/2*(N/2*(1-r)*(1-k)+m*faB)/n[6],
                 (1-r)*1/2*(N/2*r*(1-k)+m*fab)/n[6] + r*1/2*(N/2*(1-r)*(1-k)+m*faB)/n[6],
                 r*1/2*(N/2*r*(1-k)+m*fab)/n[6],
                 0, ## 7th row
                 1/2*m*fAB/n[7],
                 0,
                 1/2*m*fAB/n[7],
                 1/2*(N/2+m*fAb)/n[7],
                 1/2*m*faB/n[7],
                 1/2*(N/2*(1-k)+m*fab)/n[7] + 1/2*(N/2+m*fAb)/n[7],
                 0,
                 1/2*m*faB/n[7],
                 1/2*(N/2*(1-k)+m*fab)/n[7],
                 0, ## 8th row.
                 0, m*fAB/n[8], 0, 0, m*fAb/n[8], 0, (N+m*faB)/n[8], m*fab/n[8], 0,
                 0, ## 9th row
                 0,
                 1/2*m*fAB/n[9],
                 1/2*m*fAB/n[9],
                 0,
                 1/2*m*fAb/n[9],
                 1/2*m*fAb/n[9],
                 1/2*(N/2*(1-k)+m*faB)/n[9],
                 1/2*(N/2*(1-k)+m*fab)/n[9] + 1/2*(N/2*(1-k)+m*faB)/n[9],
                 1/2*(N/2*(1-k)+m*fab)/n[9],
                 0, 0, 0, m*fAB/n[10], 0, 0, m*fAb/n[10], 0, m*faB/n[10], (N*(1-k)+m*fab)/n[10]
  ), ncol=10, byrow = TRUE )
  
  P[is.na(P)] = 0  ## NaN caused by 0 pollen supply as divisor.
  
  fertility = fertility_grain_number(n = n, N = N, indF = indF, R=R)
  P = fertility * P  ## P is the genotype-iterating matrix. Each row has the same fertility, so can just multiple fertility.
  gout = g %*% P  ## this will give genotype of next generation.
  return(gout)
}


## iterate the nextgen_genotype function for given generations to find the state of SD. Terminate at certain gen or when change is small.
statefun2 <- function(pollen.N, inbreed.f, recomb.r, generation.g, f.Ab, killing.k, R=R){
  N = pollen.N
  inbreed.f = inbreed.f
  r = recomb.r
  k = killing.k
  
  G = c(0, 0, 0, 0, f.Ab, 0, 0, 1-f.Ab, 0, 0)  ## Ab invade aB.aB population
  #G = c((1-f.Ab)/2, 0, 0, 0, f.Ab, 0, 0, (1-f.Ab)/2, 0, 0)  ## Ab arise as a new mutation in AB.AB + aB.aB population
  AF = genotype_to_allele_freq(G)
  
  diff = 1
  gen = 0
  
  while (gen < generation.g && AF[2] > f.Ab/100 && AF[3] > (1-f.Ab)/100) {
    oldG= G
    mpool = pool_pollen_number(g = G, N = N, indF = inbreed.f, k = k)
    fpool = pool_allele_freq(g = G, r = r, k = k)
    genotype_pollen_supply = c(mpool+N, mpool+N, mpool+N, mpool+N-N*k/2, mpool+N, mpool+N-N*k/2, mpool+N-N*k/2, mpool+N, mpool+N-N*k, mpool+N-N*k)
    G = nextgen_genotype(g=G, f=fpool, m=mpool, n = genotype_pollen_supply, N=N, R=R, indF = inbreed.f, r=r, k=k)
    G = G/sum(G)
    AF = genotype_to_allele_freq(G)
    
    #diff = sum( abs(abs(oldG)-abs(G)) )
    gen = gen+1 
  }
  
  df = as.data.frame(t(genotype_to_allele_freq(G)))
  names(df) =  c('AB','Ab','aB','ab')
  
  
  state='NotAvail'
  
  delta = min (f.Ab/200, (1-f.Ab)/200)  ## arbitrary cutoff.
  
  if (df[1,]$Ab < f.Ab/100){  
    state = "removed"
  } else if ( df[1,]$Ab > 1- (1-f.Ab)/100 ){
    state = "fixed"
  } else if ( df[1,]$Ab>f.Ab + delta){
    state = "increase"
  } else if ( df[1,]$Ab<f.Ab - delta  ){
    state = "decrease"
  } else {
    state = "stable"
  }
  
  out = list('AB'= df[1,]$AB,'Ab'= df[1,]$Ab,'aB'= df[1,]$aB,'ab'= df[1,]$ab, 'state' = state, 'gen' = gen)
  return(out)
}



## eqmfun to find out the state of population


eqmfun <- function (x, f.Ab) {
  
  x = x[3:5]
  x = as.numeric(x)
  
  if (x[2]<f.Ab/100) {x[2]=0}
  if (x[3]<(1-f.Ab)/100) {x[3]=0}  ## either drive or sensitive drops to 1/100 initial freq, consider them as 0
  x = x/sum(x)
  
  delta = min (f.Ab/200, (1-f.Ab)/200) ## arbitrary cutoff.
  
  state = 'ThreeHap'
  index = which(x>1-delta)
  if (length(index)==1){
    if (index==1) {
      state='NeutralOnly'
    } else if (index==2) {
      state='DriveOnly'
    } else if (index==3) {
      state='SensitiveOnly'
    } 
  } else {
    index = which(x<delta)
    if (length(index)==1){
      if (index==1) {
        state='DriveSensitive'
      } else if (index==2) {
        state='NeutralSensitive'
      } else if (index==3) {
        state='NeutralDrive'
      } 
    }
  }
  
  state
  
}


trace.allele.freq.df <- function (indF = 0.00001, r = 0, generation.g = 50, f.Ab = 0.1, redundancy = 10){
  
  N = 4000
  k = 1
  
  G = c(0, 0, 0, 0, f.Ab, 0, 0, 1-f.Ab, 0, 0)  ## AbAb invade aB.aB population
  df = as.data.frame(t(genotype_to_allele_freq(G)))
  #dg = as.data.frame(t(G))
  
  for (gen in 1:generation.g) {
    mpool = pool_pollen_number(g = G, N = N, indF = indF, k = k)
    fpool = pool_allele_freq(g = G, r = r, k = k)
    genotype_pollen_supply = c(mpool+N, mpool+N, mpool+N, mpool+N-N*k/2, mpool+N, mpool+N-N*k/2, mpool+N-N*k/2, mpool+N, mpool+N-N*k, mpool+N-N*k)
    G = nextgen_genotype(g=G, f=fpool, m=mpool, n = genotype_pollen_supply, N=N, indF = indF, R=redundancy, r=r, k=k)
    G = G/sum(G)
    #dg = rbind(dg, as.data.frame(G))
    x = as.data.frame(t(genotype_to_allele_freq(G)))
    df = rbind(df, x)
  }
  
  names(df) =  c('AB','Ab','aB','ab')
  #names(dg) = c('AB.AB', 'AB.Ab','AB.aB','AB.ab', 'Ab.Ab','Ab.aB','Ab.ab','aB.aB','aB.ab','ab.ab')
  
  df = df[-c(1),]  ## first generation is like burn-in.
  #dg = dg[-c(1),]
  df
  #dg
}


plot.allele.freq.df <- function(df2) {
  par(mar=c(3,3,1,2),mgp=c(1.5, 0.5, 0))
  plot(1:nrow(df2), df2$Ab, type='l',ylim=c(0,1), col='firebrick4', lwd=4, yaxt='n', 
       xlab='Generations',
       ylab='Allele frequency')
  
  lines(1:nrow(df2), df2$aB, lwd=4, col='steelblue')
  lines(1:nrow(df2), df2$AB, lwd=4)
  axis(2, at=c(0,0.5,1))
  legend('topright', legend=c('Drive', 'Sensitive','Neutral'), 
         col=c('firebrick4','steelblue','black'), pch=20, cex=0.7 )
}





## 2. run statefun and plot ------
#if (TRUE) {
if (FALSE) {
df <- read.table(text ="", col.names = c('F', 'r', 'AB','Ab','aB','ab','state'))
F.seq = seq(0.01, 0.99, by=0.01)
r.seq = seq(0, 0.5, by=0.01)

for (i in F.seq){
  for (j in r.seq){
    out=statefun2(pollen.N = 4000, inbreed.f = i, recomb.r = j, generation.g = 10000, f.Ab = 0.5, killing.k = 1, R = 10)
    x = c(i, j, out$AB, out$Ab, out$aB, out$ab, out$state, out$gen)
    x = as.data.frame(t(x))
    names(x) = c('F', 'r', 'AB','Ab','aB','ab','state','gen')
    df = rbind(df, x)
  }
}



order = c('increase', 'decrease', 'stable', 'fixed', 'removed')
statecol = c('deeppink4', 'royalblue', 'forestgreen','black','grey90')
df$state = factor(df$state, levels = order)
df$F = as.numeric(as.character(df$F))
df$r = as.numeric(as.character(df$r))
df$AB = as.numeric(as.character(df$AB))
df$Ab = as.numeric(as.character(df$Ab))
df$aB = as.numeric(as.character(df$aB))
df$ab = as.numeric(as.character(df$ab))
df$gen = as.numeric(as.character(df$gen))

plot(df$F, df$r, pch=15, col=statecol[df$state], xlab = 'Selfing rate', ylab='Recombination rate')


}



## 3. trace the allele frequency to make plots ----
#if (TRUE) {
if (FALSE) {


par(mfrow=c(2,2))

  
df = trace.allele.freq.df(indF = 0.00001, r = 0, generation.g = 100, f.Ab = 0.1, redundancy = 10)
plot.allele.freq.df(df)

df = trace.allele.freq.df(indF = 0.9, r = 0, generation.g = 100)
plot.allele.freq.df(df)


df = trace.allele.freq.df(indF = 0.00001, r = 0.1, generation.g = 100)
plot.allele.freq.df(df)

df = trace.allele.freq.df(indF = 0.9, r = 0.1, generation.g = 100)
plot.allele.freq.df(df)


}












