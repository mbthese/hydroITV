data {
  int<lower=1>  N ; // # observations
  int<lower=1> S ; // # species
  vector[N] Y ; // trait value
  vector[N] DBH ; // Diameter at Breast Height
  vector[S] TWIs ; // Topographic Wetness Index of species
  vector[N] TWIis ; // Topographic Wetness Index of individual knowing species TWI_i|s
  int<lower=1, upper=S> species[N] ; // species index
}
parameters {
  vector [S] alpha_s ; // species intercepts
  vector<lower=0> [S]  betaDBH ; // DBH half-load
  real betaTWI ; // individual TWI slope
  vector[S] gammaTWI ; // species TWI slopee
  real<lower=0> sigma ; // variance
}
model {
  Y ~ normal((alpha_s[species] + betaTWI*TWIs[species] + gammaTWI[species] .* TWIis) .* 
                  (DBH ./ (betaDBH[species] + DBH)), sigma) ; // Likelihood
  alpha_s ~ normal(0, 1) ;
  betaDBH ~ normal(0,1) ;
  betaTWI ~ normal(0,1) ;
  gammaTWI ~ normal(0,1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[N] Yp = (alpha_s[species] + betaTWI*TWIs[species] + gammaTWI[species] .* TWIis) .* 
                 (DBH ./ (betaDBH[species] + DBH)) ;
  vector[N] Ydbh = (alpha_s[species] + betaTWI*mean(TWIs[species]) + gammaTWI[species]*mean(TWIis)) .* 
                 (DBH ./ (betaDBH[species] + DBH)) ;
  vector[N] Ytwii = (alpha_s[species] + betaTWI*mean(TWIs[species]) + gammaTWI[species] .* TWIis) .* 
                 (mean(DBH) ./ (betaDBH[species] + mean(DBH))) ;
}
