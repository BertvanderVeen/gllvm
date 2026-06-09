// family_va_nll.h — VA and EVA family-specific NLL contributions
//
// Include inside if((method<1)||(method>1)){ ... } in a TMB objective function.
// Requires the same outer-scope variables as family_nll.h (see that file for contract).
// Handles method<1 (VA) and method>1 (EVA) paths; does NOT close the enclosing if block.
//
// Variables declared internally: idx, has12
    int idx = 0; // initialize indexing for zeta
    
    bool has12 = false;
    for (int j = 0; j < family.size(); ++j) {
      if (family(j) == 12) {
        has12 = true;
        break;
      }
    }
    
    // Distributions (family)
    // truep =p-p_betaH
  for (int j=0; j<(truep);j++){
    
    switch (family(j)) {
    
    case POISSON: {//poisson family 0
      for (int i=0; i<n; i++) {
        // for (int j=0; j<p;j++){
        if(!gllvmutils::isNA(y(i,j)))nll -= dpois(y(i,j), exp(eta(i,j)+cQ(i,j)), true)-y(i,j)*cQ(i,j);
        // }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0r(i)/sigma,2))*random(0);
      }
      break;
    }
      
    case NEG_BINOMIAL: {//NB family 1
      if(method<1){//NB VA
        if(extra(j) == 0){
          //nb2
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p;j++){
              // nll -= Type(gllvm::dnegbinva(y(i,j), eta(i,j), iphi(j), cQ(i,j)));
              // if(!gllvmutils::isNA(y(i,j)))nll -= y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*cQ(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
              if(!gllvmutils::isNA(y(i,j))){
                nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
                Type log_term_1 = gllvmutils::log1plus(exp(eta(i,j)-lg_phi(j)));
                Type log_term_2 = gllvmutils::log1plus(exp(eta(i,j)-cQ(i,j)-lg_phi(j)));
                nll -= (y(i,j)+iphi(j))*(log_term_1-log_term_2)-(y(i,j)+iphi(j))*cQ(i,j);
              }
            // }
          }
        }else if(extra(j) == 1){
          //nb1
          const double gamma = 0.57721566490153286060651209008240243;
          
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p;j++){
              if(!gllvmutils::isNA(y(i,j))){
                nll -= -(y(i,j) + exp(eta(i,j) + cQ(i,j))*iphi(j))*gllvmutils::log1plus(iphi(j)) - lfactorial(y(i,j)) + iphi(j)*exp(eta(i,j) + cQ(i,j))*(gamma + lg_phi(j)) + eta(i,j) +lg_phi(j) - gamma*exp(eta(i,j)+2*cQ(i,j))*iphi(j) - lgamma(exp(eta(i,j)+2*cQ(i,j))*iphi(j)+1.0) + lgamma(y(i,j) + exp(eta(i,j)+cQ(i,j))*iphi(j));
              }
            // }
          }
        }
      } else if (method>1) { // NB EVA
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p;j++){
            if(!gllvmutils::isNA(y(i,j))){
              nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
              nll += (((iphi(j)+y(i,j)) / (iphi(j)+exp(eta(i,j)))) * exp(eta(i,j)) - ((iphi(j)+y(i,j))*pow(iphi(j)+exp(eta(i,j)),-2))*pow(exp(eta(i,j)),2)) * cQ(i,j);
            }
            // nll += gllvm::nb_Hess(y(i,j), eta(i,j), iphi(j)) * cQ(i,j);
            // nll -= lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) - lgamma(y(i,j)+1) + y(i,j)*eta(i,j) + iphi(j)*log(iphi(j))-(y(i,j)+iphi(j))*log(exp(eta(i,j))+iphi(j));
            // nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
            // nll += (((iphi(j)+y(i,j)) / (iphi(j)+exp(eta(i,j)))) * exp(eta(i,j)) - ((iphi(j)+y(i,j))*pow(iphi(j)+exp(eta(i,j)),-2))*pow(exp(eta(i,j)),2)) * cQ(i,j);
          // }
        }
      }
      break;
    }
      
    case BINOMIAL: {//binomial family 2
      if(method<1) {//binomial VA
        if (extra(j) == 0) { //logit
          // for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              // Type a = 0.5*sqrt(squeeze(eta(i,j)*eta(i,j) + 2*cQ(i,j)));//bound it because derivative logcosh = tanh(10) = 1 flattens
              // Type a = sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
              // Type softplus_neg_a = CppAD::CondExpGt(a, Type(15), exp(-a), gllvmutils::log1plus(exp(-a)));
              
              //Type b = CppAD::CondExpGt(a, 10, a/8-log(2.0), gllvmutils::logcosh(0.5*sqrt(a)));
              // Type b = CppAD::CondExpGt(a, 10, 10, gllvmutils::logcosh(0.5*sqrt(squeeze(eta(i,j)*eta(i,j) + 2*cQ(i,j)))));
              // nll -= (y(i,j)-Ntrials(i,j)/2)*eta(i,j) - Ntrials(i,j)*(0.5*a+softplus_neg_a);//logspace_add(Type(0),-a));//gllvmutils::log1plus(exp(-a)));//log(invlogit(a)));//Ntrials(i,j)*gllvmutils::logcosh(a);//-0.5*tanh(0.5)*(eta(i,j)*eta(i,j)+2*cQ(i,j))+0.5*tanh(a)*(eta(i,j)*eta(i,j)+2*cQ(i,j));
              Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
              // Type wij = 0.5*gllvmutils::hypo(eta(i,j), sqrt(2*cQ(i,j)));
              nll -= (y(i,j)-Ntrials(i, j)*0.5)*eta(i,j) - Ntrials(i, j)*logspace_add(wij, -wij);
               // nll -= (y(i,j)-Ntrials(i, j)*0.5)*eta(i,j) - Ntrials(i, j)*gllvmutils::log1plus(exp(-2*wij));
              if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
              }
            }
            // nll += n*Ntrials(i,j)*log(2.0);
          // }
        }else if(extra(j)==1){//probit
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p;j++){
            mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            if(!gllvmutils::isNA(y(i,j))){
              nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(i,j)-y(i,j));
              nll += cQ(i,j)*Ntrials(i,j);
              if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
              }
            }
          // }
        }
        }else if(extra(j)==2){//cloglog
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p;j++){
              mu(i,j) = exp(eta(i,j)+cQ(i,j));
              if(!gllvmutils::isNA(y(i,j))){
                nll -= y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }
            // }
          }
        }
      } else if (method>1) { // Binomial EVA
        if (extra(j) == 0) { // logit
          //Type mu_prime;
          //CppAD::vector<Type> z(4);
          
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p; j++) {
              if (!gllvmutils::isNA(y(i,j))) {
                Type log_1mp = -CppAD::CondExpLe(eta(i,j), Type(18.), gllvmutils::log1plus(exp(eta(i,j))), eta(i,j));
                Type log_p = -CppAD::CondExpLe(-eta(i,j), Type(18.), gllvmutils::log1plus(exp(-eta(i,j))), -eta(i,j));
                nll -= y(i,j)*log_p + (Type(1.)-y(i,j))*log_1mp;
                nll += gllvmutils::mfexp(log_1mp + log_p)*cQ(i,j);
              }
          // nll -= gllvm::dbinom_logit_eva(y(i,j), eta(i,j), cQ(i,j));
              
          //    mu(i,j) = 0.0;
          //    mu_prime = 0.0;
              
          //    z[0] = eta(i,j);
          //    z[1] = 0;
          //    z[2] = 1/(1+exp(-z[0]));
          //    z[3] = exp(z[0])/(exp(z[0])+1);
              
          //    mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
          //    mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
          //    mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
              
          //    mu_prime = mu(i,j) * (1-mu(i,j));
          //    if(!gllvmutils::isNA(y(i,j))){
          //      nll -= y(i,j) * eta(i,j) + log(1-mu(i,j));
          //      nll += mu_prime*cQ(i,j);
            // }
          }
        } else if (extra(j) == 1) { // probit
          // Type etaP;
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p; j++) {
              if(!gllvmutils::isNA(y(i,j))){
                Type etaD =  dnorm(eta(i,j), Type(0), Type(1), 1);   // normal density evaluated at eta(i,j)
                Type logit_p = gllvmutils::logit_pnorm(eta(i,j));
  
                Type log_p = -CppAD::CondExpLe(-logit_p, Type(18.0), gllvmutils::log1plus(exp(-logit_p)), -logit_p); 
                Type log_1mp = -CppAD::CondExpLe(logit_p, Type(18.0), gllvmutils::log1plus(exp(logit_p)), logit_p); 
                
                nll -= y(i,j)*log_p + (Type(1.0)-y(i,j))*log_1mp;
                //Type tmp = CppAD::CondExpLt(logit_p, Type(0.0), -logit_p + 2.*log_1mp, logit_p + 2.*log_p);
                
                //nll -= -pow(y(i,j)-exp(log_p),2)* exp(2*tmp+2*etaD)*cQ(i,j) - (y(i,j)-exp(log_p))*eta(i,j)*exp(tmp+etaD)*cQ(i,j);
                nll -= ((y(i,j)*(gllvmutils::mfexp(log_p + etaD)*(-eta(i,j))-gllvmutils::mfexp(2.*etaD))*gllvmutils::mfexp(2.*log_1mp) + (1.-y(i,j))*(gllvmutils::mfexp(log_1mp+etaD)*eta(i,j)-gllvmutils::mfexp(2.*etaD))*gllvmutils::mfexp(2.*log_p) )/(gllvmutils::mfexp(2*log_p)*(gllvmutils::mfexp(2*log_p)-2*gllvmutils::mfexp(log_p)+1)))*cQ(i,j);
                //etaP = pnorm_approx(Type(eta(i,j)));
                
                //etaP = Type(CppAD::CondExpEq(etaP, Type(1), etaP-Type(1e-12), etaP));//check if on the boundary
                //etaP = Type(CppAD::CondExpEq(etaP, Type(0), etaP+Type(1e-12), etaP));//check if on the boundary
                
                //nll -= y(i,j)*log(etaP) + (1-y(i,j))*log(1-etaP); //
                //Type etaD =  dnorm(Type(eta(i,j)), Type(0), Type(1), true);   // log normal density evaluated at eta(i,j)
                //nll -= ((y(i,j)*(etaP*exp(etaD)*(-eta(i,j))-pow(exp(etaD),2))*pow(1-etaP,2) + (1-y(i,j))*((1-etaP)*exp(etaD)*eta(i,j)-pow(exp(etaD),2))*pow(etaP,2) )/(etaP*etaP*(etaP*etaP-2*etaP+1)))*cQ(i,j);
              }
            // }
          }
        }
      }
      break;
    }
    
    case GAUSSIAN: {//gaussian family 3
      for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j))) - log(M_PI)/2;
      }
      break;
    }
    
    case GAMMA: {//gamma family 4
      for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));
      }
      break;
    } 
      
    case TWEEDIE: {// Tweedie family 5
      if(method >1){ // Tweedie EVA
        Type ePower1 = invlogit(ePower) + Type(1);
        for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j))){
              // Tweedie log-likelihood:
              nll -= dtweedie(y(i,j), exp(eta(i,j)), iphi(j), ePower1, true);
              if (y(i,j) == 0) {
                // Hessian-trace part:
                nll += (1/iphi(j)) * (2-ePower1)*exp(2*eta(i,j))*exp(-ePower1*eta(i,j)) * cQ(i,j);
              } else if (y(i,j) > 0) {
                nll -= (1/iphi(j)) * (y(i,j)*(1-ePower1)*exp((1-ePower1)*eta(i,j)) - (2-ePower1)*exp((2-ePower1)*eta(i,j))) * cQ(i,j);
              }
            }
        }
      }else if(method <1){ // Tweedie VA
        Type ePower1 = invlogit(ePower) + Type(1);
        for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j))){
              Type p1 = ePower1 - 1.0, p2 = 2.0 - ePower1;
              Type ans = -pow(exp(eta(i,j)+p2*cQ(i,j)), p2)/(iphi(j)*p2);
               if(y(i,j)>0){
                CppAD::vector<Type> tx(4);
                tx[0] = y(i,j);
                tx[1] = iphi(j);
                tx[2] = ePower1;
                tx[3] = 0;
                ans += atomic::tweedie_logW(tx)[0];
                ans += -y(i,j) / (iphi(j) * p1 * pow(exp(eta(i,j)-p1*cQ(i,j)), p1)) - log(y(i,j));
              }
               nll -= ans;
            }
        }
      }
      break;
    }
    
    case ZIP: { //ZIP family 6
      Type iphij = iphi(j)/(1+iphi(j));
      Type pVA;
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            if(y(i,j)>0){
              nll -= log(1-iphij)+y(i,j)*eta(i,j)-exp(eta(i,j)+cQ(i,j))-lfactorial(y(i,j));
            }else{
              pVA = exp(log(-iphij+1)-exp(eta(i,j)+cQ(i,j))-log((1-iphij)*exp(-exp(eta(i,j)+cQ(i,j)))+iphij));
              pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
              pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
              nll -= log(iphij)-log(1-pVA);
            }
          }
        }
      break;
    }
      
    case ORDINAL: {//ordinal family 7
      if(zetastruc == 1){//ordinal with species specific cutoffs
        int ymax =  CppAD::Integer(y.maxCoeff());
        int K = ymax - 1;
        
        // matrix <Type> zetanew(p,K);
        vector <Type> zetanew(K);
        zetanew.setZero();
        
        // int idx = 0; // indexing for zeta moved before for j
          int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
          int Kj = ymaxj - 1;
          if(Kj>1){
            for(int k=0; k<(Kj-1); k++){
              zetanew(k+1) = zeta.segment(idx,k+1).array().exp().sum();
            }
            idx += Kj-1; 
          }
        
        if (method<1) { // VA
          if(extra(j) == 0){ //va logit
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                if(!gllvmutils::isNA(y(i,j))){
                  int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                  //yik = 1 if yi >=k and 0 otherwise
                  // p(yik = 0) for k<y(i,j)
                  for (int l=0; l<CppAD::Integer(y(i,j)-1); l++) {
                    Type wij = 0.5*sqrt((zetanew(l)-eta(i,j))*(zetanew(l)-eta(i,j)) + 2*cQ(i,j));
                    nll -= -0.5*(zetanew(l)-eta(i,j)) - logspace_add(wij, -wij);
                    // Type wij = 0.5*sqrt((zetanew(j,l)-eta(i,j))*(zetanew(j,l)-eta(i,j)) + 2*cQ(i,j));
                    // nll -= -0.5*(zetanew(j,l)-eta(i,j)) - logspace_add(wij, -wij);
                  }
                  // p(yik = 1)  for k>= y(i,j)
                  for (int l=CppAD::Integer(y(i,j)-1); l< (ymaxj -1); l++) {
                    Type wij = 0.5*sqrt((zetanew(l)-eta(i,j))*(zetanew(l)-eta(i,j)) + 2*cQ(i,j));
                    nll -= 0.5*(zetanew(l)-eta(i,j)) - logspace_add(wij, -wij);
                    // Type wij = 0.5*sqrt((zetanew(j,l)-eta(i,j))*(zetanew(j,l)-eta(i,j)) + 2*cQ(i,j));
                    // nll -= 0.5*(zetanew(j,l)-eta(i,j)) - logspace_add(wij, -wij);
                  }
                }
              // }
            }
          }else{//va probit
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                if(!gllvmutils::isNA(y(i,j))){
                  int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                  //minimum category
                  if(y(i,j)==1){
                    mu(i,j) = pnorm(zetanew(0) - eta(i,j), Type(0), Type(1));
                    // mu(i,j) = pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1));
                    mu(i,j) = Type(CppAD::CondExpLt(mu(i,j), Type(1e-12), mu(i,j)+Type(1e-12), mu(i,j)));
                    nll -= log(mu(i,j));
                  }else if(y(i,j)==ymaxj){
                    //maximum category
                    int idxj = ymaxj-2;
                    mu(i,j) = pnorm(zetanew(idxj) - eta(i,j), Type(0), Type(1));
                    // mu(i,j) = pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1));
                    mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));
                    nll -= log(1 - mu(i,j));
                  }else if(ymaxj>2){
                    for (int l=2; l<ymaxj; l++) {
                      if((y(i,j)==l) && (l != ymaxj)){
                        mu(i,j) = pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1));
                        // mu(i,j) = pnorm(zetanew(j,l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,l-2)-eta(i,j), Type(0), Type(1));
                        mu(i,j) = Type(CppAD::CondExpLt(mu(i,j), Type(1e-12), mu(i,j)+Type(1e-12), mu(i,j)));
                        nll -= log(mu(i,j));
                      }
                    }
                  }
                  
                  nll += cQ(i,j);
                }
                //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));//
              // }
            }
          }
        } else if (method>1) { // EVA ordinal
          if (extra(j)==0) { // logit
            for (int i=0; i<n; i++) {
                if(!gllvmutils::isNA(y(i,j))){
                  int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                  //minimum category
                  if(y(i,j)==1){
                    nll -= -gllvmutils::log1plus(exp(eta(i,j))-zetanew(0));
                    nll -= -dlogis(zetanew(0), eta(i,j), Type(1), 0)*cQ(i,j);
                    // //nll -= -logspace_add(Type(0), eta(i,j)-zetanew(j,0));
                    // nll -= -dlogis(zetanew(j,0), eta(i,j), Type(1), 0)*cQ(i,j);
                  }else if(y(i,j)==ymaxj){
                    //maximum category
                    int idxj = ymaxj-2;
                    nll -= -gllvmutils::log1plus(exp(zetanew(idxj)-eta(i,j)));
                    nll -= -dlogis(zetanew(idxj), eta(i,j), Type(1), 0)*cQ(i,j);
                  }else if(ymaxj>2){
                    for (int l=2; l<ymaxj; l++) {
                      if((y(i,j)==l) && (l != ymaxj)){
                        nll -= logspace_sub(-logspace_add(Type(0), eta(i,j)-zetanew(l-1)), -logspace_add(Type(0), eta(i,j)-zetanew(l-2)));
                        nll -= -dlogis(zetanew(l-2),eta(i,j),Type(1),0)*cQ(i,j);
                        nll -= -dlogis(zetanew(l-1),eta(i,j),Type(1),0)*cQ(i,j);
                      }
                    }
                  }
                }
              
            }
          }
        }
        
      } else if(zetastruc==0){//ordinal with common cutoffs
        // int ymax =  CppAD::Integer(y.maxCoeff());
        // int K = ymax - 1;
        
        int ymax = CppAD::Integer(y.col(j).maxCoeff());
        int K = ymax - 1;
        
        vector <Type> zetanew(K);
        zetanew.setZero();

        if(has12) idx = 2; // start from 2 if there are orderedBeta columns in the model
        for(int k=0; k<(K-1); k++){
          zetanew(k+1) = zeta.segment(idx, k+1).array().exp().sum();//second cutoffs must be positive
        }
        
        if (method<1) {
          if(extra(j) == 0){ // va logit
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                if(!gllvmutils::isNA(y(i,j))){
                  int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                  //yik = 1 if yi >=k and 0 otherwise
                  // p(yik = 0) for k<y(i,j)
                  for (int l=0; l<CppAD::Integer(y(i,j)-1); l++) {
                    Type wij = 0.5*sqrt((zetanew(l)-eta(i,j))*(zetanew(l)-eta(i,j)) + 2*cQ(i,j)); 
                    nll -= -0.5*(zetanew(l)-eta(i,j)) - logspace_add(wij, -wij);
                  }
                  // p(yik = 1)  for k>= y(i,j)
                  for (int l=CppAD::Integer(y(i,j)-1); l< (ymaxj -1); l++) {
                    Type wij = 0.5*sqrt((zetanew(l)-eta(i,j))*(zetanew(l)-eta(i,j)) + 2*cQ(i,j)); 
                    nll -= 0.5*(zetanew(l)-eta(i,j)) - logspace_add(wij, -wij);
                  }
                }
              // }
            }
          }else{ // va probit
            for (int i=0; i<n; i++) {
                if(!gllvmutils::isNA(y(i,j))){
                  //minimum category
                  if(y(i,j)==1){
                    mu(i,j) = pnorm(zetanew(0) - eta(i,j), Type(0), Type(1));
                    mu(i,j) = Type(CppAD::CondExpLt(mu(i,j), Type(1e-12), mu(i,j)+Type(1e-12), mu(i,j)));
                    nll -= log(mu(i,j));
                  }else if(y(i,j)==ymax){
                    //maximum category
                    int idxj = ymax-2;
                    mu(i,j) = pnorm(zetanew(idxj) - eta(i,j), Type(0), Type(1));
                    mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));
                    nll -= log(1 - mu(i,j));
                  }else if(ymax>2){
                    for (int l=2; l<ymax; l++) {
                      if((y(i,j)==l) && (l != ymax)){
                        mu(i,j) = pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1));
                        mu(i,j) = Type(CppAD::CondExpLt(mu(i,j), Type(1e-12), mu(i,j)+Type(1e-12), mu(i,j)));
                        nll -= log(mu(i,j));
                      }
                    }
                  }
                  nll += cQ(i,j);
                }
              // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0r(i)/sigma,2))*random(0);
            }
          }
        } else if (method>1) {
          if (extra(j)==0) {
            for (int i=0; i<n; i++) {
              // for (int j=0; j<p; j++) {
                if (y(i,j)==1) { // min category
                  nll -= -logspace_add(Type(0),eta(i,j)-zetanew(0));
                  nll -= -dlogis(zetanew(0,0),eta(i,j),Type(1),0)*cQ(i,j);
                } else if(y(i,j)==ymax) { // max category
                  int idxj = ymax-2;
                  nll -= -logspace_add(Type(0),zetanew(idxj)-eta(i,j));
                  nll -= -dlogis(zetanew(idxj),eta(i,j),Type(1),0)*cQ(i,j);
                } else if(ymax>2) {
                  for (int l=2; l<ymax; l++) {
                    if ((y(i,j)==l) && (l != ymax)) {
                      //nll(i,j) -= logspace_sub(-log1plus(exp(eta(i,j)-zetanew(0,l-1))),-log1plus(exp(eta(i,j)-zetanew(0,l-2))));
                      nll -= logspace_sub(-logspace_add(Type(0),eta(i,j)-zetanew(l-1)),-logspace_add(Type(0),eta(i,j)-zetanew(l-2)));
                      nll -= (-dlogis(zetanew(l-2),eta(i,j),Type(1),0) - dlogis(zetanew(l-1),eta(i,j),Type(1),0))*cQ(i,j);
                    }
                  }
                }
              // }
            }
          }
        }
      }
      break;
    }
    
    case EXPONENTIAL: {// exp family 8
      for (int i=0; i<n; i++) {
        // for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j))) nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) );
        // }
      }
      break;
    } 
    
    case BETA: { // Beta family 9  (EVA only)
      Type mu_prime;
      Type mu_prime2;
      CppAD::vector<Type> z;
      if(extra(j)==0){
        z = CppAD::vector<Type> (4);
      }
      CppAD::vector<Type> a(2);
      CppAD::vector<Type> b(2);
      CppAD::vector<Type> aa;
      CppAD::vector<Type> bb;
      Type dig_a;
      Type dig_b;
      Type trig_a;
      Type trig_b;
      for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            // define mu, mu' and mu''
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            mu_prime2 = 0.0;
            if (extra(j) == 0) { // logit
              
              z[0] = eta(i,j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
              
              mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              mu_prime = mu(i,j) * (1-mu(i,j));
              mu_prime2 = mu_prime * (1-2*mu(i,j));
              
            } else if (extra(j) == 1) { // probit
              mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
              mu_prime = dnorm(eta(i,j), Type(0), Type(1));
              mu_prime2 = (-eta(i,j))*mu_prime;
            }
            a[0] = mu(i,j)*iphi(j);
            a[1] = 1;
            b[0] = (1-mu(i,j))*iphi(j);
            b[1] = 1;
            aa = a;
            bb = b;
            aa[1] = 2;
            bb[1] = 2;
            dig_a = Type(atomic::D_lgamma(a)[0]);
            dig_b = Type(atomic::D_lgamma(b)[0]);
            trig_a = Type(atomic::D_lgamma(aa)[0]);
            trig_b = Type(atomic::D_lgamma(bb)[0]);
            
            nll -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
            nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
            nll -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
            
          }
      }
      break;
    }
    
    case BETA_HURDLE: {// hurdle Beta family 10
      if(method<1)  { // hurdle Beta VA-EVA hybrid
        Type mu_prime;
        Type mu_prime2;
        CppAD::vector<Type> z;
        if(extra(j)==0){
          z = CppAD::vector<Type> (4);
        }
        CppAD::vector<Type> a(2);
        CppAD::vector<Type> b(2);
        CppAD::vector<Type> aa;
        CppAD::vector<Type> bb;
        Type dig_a;
        Type dig_b;
        Type trig_a;
        Type trig_b;
        for (int i=0; i<n; i++) {
          // for (int j=0; j<truep; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              // define mu, mu' and mu''
              mu(i,j) = 0.0;
              mu_prime = 0.0;
              mu_prime2 = 0.0;
              if (extra(j) == 0) { // logit
                // mu(i,truep+j) = Type(CppAD::CondExpGe(eta(i,truep+j), type(0), 1/(1+exp(-eta(i,truep+j)) ), exp(eta(i,truep+j))/(exp(eta(i,truep+j))+1) ));
                z[0] = eta(i,truep+j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
                
                mu(i,truep+j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                
                z[0] = eta(i,j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
                
                mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                mu_prime = mu(i,j) * (1-mu(i,j));
                mu_prime2 = mu_prime * (1-2*mu(i,j));
                
              } else if (extra(j) == 1) { // probit
                mu(i,truep+j) = pnorm(eta(i,truep+j), Type(0), Type(1));
                mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
                mu_prime = dnorm(eta(i,j), Type(0), Type(1));
                mu_prime2 = (-eta(i,j))*mu_prime;
              }
              
              if(y(i,j)==0){
                nll -= log( 1.0 - mu(i,truep+j) ) - cQ(i,truep+j);
              } else{
                nll -= log( mu(i,truep+j) ) - cQ(i,truep+j);
                
                a[0] = mu(i,j)*iphi(j);
                a[1] = 1;
                b[0] = (1-mu(i,j))*iphi(j);
                b[1] = 1;
                aa = a;
                bb = b;
                aa[1] = 2;
                bb[1] = 2;
                dig_a = Type(atomic::D_lgamma(a)[0]);
                dig_b = Type(atomic::D_lgamma(b)[0]);
                trig_a = Type(atomic::D_lgamma(aa)[0]);
                trig_b = Type(atomic::D_lgamma(bb)[0]);
                
                nll -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
                nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
                nll -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
              }
              
            }
            
          // }
        }
        
      } else if (method>1) { // hurdle beta EVA
        
        Type mu_prime;
        Type mu_prime2;
        Type mu0_prime;
        Type mu0_prime2;
        
        CppAD::vector<Type> z;
        if(extra(j)==0){
          z = CppAD::vector<Type> (4);
        }
        CppAD::vector<Type> a(2);
        CppAD::vector<Type> b(2);
        CppAD::vector<Type> aa;
        CppAD::vector<Type> bb;
        Type dig_a;
        Type dig_b;
        Type trig_a;
        Type trig_b;
        
        for (int i=0; i<n; i++) {
          // for (int j=0; j<truep; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              // define mu, mu' and mu''
              mu(i,j) = 0.0;
              mu_prime = 0.0;
              mu_prime2 = 0.0;
              mu0_prime = 0.0;
              mu0_prime2 = 0.0;
              if (extra(j) == 0) { // logit
                // mu(i,truep+j) = Type(CppAD::CondExpGe(eta(i,truep+j), type(0), 1/(1+exp(-eta(i,truep+j)) ), exp(eta(i,truep+j))/(exp(eta(i,truep+j))+1) ));
                z[0] = eta(i,truep+j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
                
                mu(i,truep+j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                
                z[0] = eta(i,j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
                
                mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                mu_prime = mu(i,j) * (1-mu(i,j));
                mu_prime2 = mu_prime * (1-2*mu(i,j));
                
                mu0_prime = mu(i,truep+j) * (1-mu(i,truep+j));
                mu0_prime2 = mu0_prime * (1-2*mu(i,truep+j));
                
              } else if (extra(j) == 1) { // probit
                mu(i,truep+j) = pnorm(eta(i,truep+j), Type(0), Type(1));
                mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
                mu_prime = dnorm(eta(i,j), Type(0), Type(1));
                mu_prime2 = (-eta(i,j))*mu_prime;
                
                mu0_prime = dnorm(eta(i,truep+j), Type(0), Type(1));
                mu0_prime2 = (-eta(i,truep+j))*mu0_prime;
              }
              
              if(y(i,j)==0){
                nll -= log( 1.0 - mu(i,truep+j) );
                //nll -= -dlogis(Type(0), eta(i,truep+j), Type(1), 0)*cQ(i,truep+j);
                nll -= -(mu0_prime2 * (1-mu(i,truep+j)) + pow(mu0_prime,2))/pow(1-mu(i,truep+j),2) * cQ(i,truep+j);            
              } else{
                nll -= log( mu(i,truep+j) );
                //nll -= -dlogis(eta(i,truep+j), Type(0.0), Type(1), 0)*cQ(i,truep+j);
                nll -= (mu(i,truep+j)*mu0_prime2-pow(mu0_prime,2))/pow(mu(i,truep+j),2) * cQ(i,truep+j);
                
                a[0] = mu(i,j)*iphi(j);
                a[1] = 1;
                b[0] = (1-mu(i,j))*iphi(j);
                b[1] = 1;
                aa = a;
                bb = b;
                aa[1] = 2;
                bb[1] = 2;
                dig_a = Type(atomic::D_lgamma(a)[0]);
                dig_b = Type(atomic::D_lgamma(b)[0]);
                trig_a = Type(atomic::D_lgamma(aa)[0]);
                trig_b = Type(atomic::D_lgamma(bb)[0]);
                
                nll -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
                nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
                nll -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
              }
            }
          // }
        }
      }
      break;
    }
    
    
    case ZINB: { // ZINB family 11
      Type iphij = iphi(j)/(1+iphi(j));
      Type iphiZINB = exp(lg_phiZINB(j));
      Type pVA;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            if(y(i,j)>0){
              nll -= log(1-iphij)+y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphiZINB)*log(iphiZINB+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphiZINB) - iphiZINB*cQ(i,j) + iphiZINB*log(iphiZINB) - lgamma(iphiZINB) -lfactorial(y(i,j));
            }else{
              pVA = exp(log(1-iphij)- iphiZINB*log(iphiZINB+exp(eta(i,j)-cQ(i,j))) + lgamma(iphiZINB) - iphiZINB*cQ(i,j) + iphiZINB*log(iphiZINB) - lgamma(iphiZINB)-log((1-iphij)*exp(- iphiZINB*log(iphiZINB+exp(eta(i,j)-cQ(i,j))) + lgamma(iphiZINB) - iphiZINB*cQ(i,j) + iphiZINB*log(iphiZINB) - lgamma(iphiZINB))+iphij));
              pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
              pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
              nll -= log(iphij)-log(1-pVA);
            }
          }
        }
      }
      break;
    } 
    
    case ORDERED_BETA: {// ordered Beta 12
      vector <Type> zetacutoffnew(2);
      zetacutoffnew.setZero();
      
      if(zetastruc==0){ // common cutoffs
        zetacutoffnew(0)= zeta(0);
        zetacutoffnew(1)= exp(zeta(1));
      } else { // species specific cutoffs
        zetacutoffnew(0)= zeta(idx);
        zetacutoffnew(1)= exp(zeta(idx+1));
        idx += 2;
      }
      if(method<1) { // ordered Beta VA-EVA hybrid
        if(extra(j)==1){
          //probit
        Type mu_prime;
        Type mu_prime2;
        CppAD::vector<Type> a(2);
        CppAD::vector<Type> b(2);
        CppAD::vector<Type> aa;
        CppAD::vector<Type> bb;
        Type dig_a;
        Type dig_b;
        Type trig_a;
        Type trig_b;
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              // define mu, mu' and mu''
              mu(i,j) = 0.0;
              mu_prime = 0.0;
              mu_prime2 = 0.0;
              // probit link
              if((y(i,j)==0)){
                // mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
                // nll -= log(pow(1.0 - pnorm(zetacutoffnew(j,1) - eta(i,j), Type(0), Type(1)), y(i,j)) * pow(pnorm(zetacutoffnew(j,0) - eta(i,j), Type(0), Type(1)),(1-y(i,j)))) - cQ(i,j);
                mu(i,j) = pnorm(zetacutoffnew(0) - eta(i,j), Type(0), Type(1));
                mu(i,j) = CppAD::CondExpLt(mu(i,j), Type(1e-12), mu(i,j)+1e-12, mu(i,j));
                nll -= (1-y(i,j))*log(mu(i,j)) - cQ(i,j); //
              } else if((y(i,j)==1)){
                mu(i,j) = pnorm(zetacutoffnew(1) - eta(i,j), Type(0), Type(1));
                mu(i,j) = CppAD::CondExpLt(mu(i,j), Type(1.0), mu(i,j), mu(i,j)-1e-12);
                nll -= y(i,j)*log(1.0 - mu(i,j)) - cQ(i,j); //
              } else{
                // if (extra(j) == 1) { // probit
                // if(zetacutoff.size()>p) {
                mu(i,j) = pnorm(zetacutoffnew(1) - eta(i,j), Type(0), Type(1)) - pnorm(zetacutoffnew(0) - eta(i,j), Type(0), Type(1));
                mu(i,j) = CppAD::CondExpGt(mu(i,j), Type(1e-12), mu(i,j), mu(i,j)+1e-12);  
                nll -= log(mu(i,j)) - cQ(i,j); //
                  // Type a1 = pnorm(zetacutoffnew(1) - eta(i,j), Type(0), Type(1)) - pnorm(zetacutoffnew(0) - eta(i,j), Type(0), Type(1));
                  // a1 = CppAD::CondExpLe(a1, Type(1.0), a1, a1-1e-12);  
                  // nll -= log(a1) - cQ(i,j); //
                // } else { // Case where there is no upperbound, atm not used 
                //   mu(i,j) = pnorm(zetacutoffnew(0) - eta(i,j), Type(0), Type(1));
                //   mu(i,j) = CppAD::CondExpLe(mu(i,j), Type(1.0), mu(i,j), mu(i,j)-1e-12);
                //   nll -= log(1 - mu(i,j)) - cQ(i,j); //
                // }
                mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
                mu_prime = dnorm(eta(i,j), Type(0), Type(1));
                mu_prime2 = (-eta(i,j))*mu_prime;
                // }
                a[0] = mu(i,j)*iphi(j);
                a[1] = 1;
                b[0] = (1-mu(i,j))*iphi(j);
                b[1] = 1;
                aa = a;
                bb = b;
                aa[1] = 2;
                bb[1] = 2;
                dig_a = Type(atomic::D_lgamma(a)[0]);
                dig_b = Type(atomic::D_lgamma(b)[0]);
                trig_a = Type(atomic::D_lgamma(aa)[0]);
                trig_b = Type(atomic::D_lgamma(bb)[0]);
                
                nll -= dbeta(y(i,j), Type(a[0]), Type(b[0]), 1);
                nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
                nll -= iphi(j) * mu_prime2 * (log(y(i,j)) - log(1-y(i,j))) * cQ(i,j) ;
              }
              
            }
          // }
        }
        }else if(extra(j)==0){
          //logit
          Type mu_prime;
          Type mu_prime2;
          CppAD::vector<Type> z;
          z = CppAD::vector<Type> (4);
          CppAD::vector<Type> a(2);
          CppAD::vector<Type> b(2);
          CppAD::vector<Type> aa;
          CppAD::vector<Type> bb;
          Type dig_a;
          Type dig_b;
          Type trig_a;
          Type trig_b;
          
          for (int i=0; i<n; i++) {
            // for (int j=0; j<p; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              // logit link
              if((y(i,j)==0)){
                  Type wij = 0.5*sqrt((zetacutoffnew(0)-eta(i,j))*(zetacutoffnew(0)-eta(i,j)) + 2*cQ(i,j));
                  nll -= 0.5*(zetacutoffnew(0)-eta(i,j)) - logspace_add(wij, -wij);
              } else if((y(i,j)==1)){
                Type wij = 0.5*sqrt((eta(i,j)-zetacutoffnew(1))*(eta(i,j)-zetacutoffnew(1)) + 2*cQ(i,j));
                nll -= 0.5*(eta(i,j)-zetacutoffnew(1)) - logspace_add(wij, -wij);
              } else{
                nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(0), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(0))), eta(i,j)-zetacutoffnew(0));
                nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(1), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(1))), eta(i,j)-zetacutoffnew(1));
                nll -= eta(i,j) - zetacutoffnew(0);
                nll -= CppAD::CondExpLe(zetacutoffnew(1)-zetacutoffnew(0), log(Type(2.)), log(-gllvmutils::expminus1(zetacutoffnew(0)-zetacutoffnew(1))),  gllvmutils::log1plus(-exp(zetacutoffnew(0)-zetacutoffnew(1))));
                nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(0), eta(i,j), Type(1), 1))*cQ(i,j); 
                nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(1), eta(i,j), Type(1), 1))*cQ(i,j);
              
                CppAD::vector<Type> z(4);
                z[0] = eta(i,j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
                
                mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                mu_prime = mu(i,j) * (1-mu(i,j));
                mu_prime2 = mu_prime * (1-2*mu(i,j));
                
                a[0] = mu(i,j)*iphi(j);
                a[1] = 1;
                b[0] = (1-mu(i,j))*iphi(j);
                b[1] = 1;
                aa = a;
                bb = b;
                aa[1] = 2;
                bb[1] = 2;
                dig_a = Type(atomic::D_lgamma(a)[0]);
                dig_b = Type(atomic::D_lgamma(b)[0]);
                trig_a = Type(atomic::D_lgamma(aa)[0]);
                trig_b = Type(atomic::D_lgamma(bb)[0]);
                
                nll -= dbeta(y(i,j), Type(a[0]), Type(b[0]), 1);
                nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
                nll -= iphi(j) * mu_prime2 * logit(y(i,j)) * cQ(i,j);
              }
              
            }
            // }
          }
        }
        
      } else if (method>1) {  // Ordered beta EVA

        Type mu_prime;
        Type mu_prime2;
        CppAD::vector<Type> z;
        if(extra(j)==0){
          z = CppAD::vector<Type> (4);
        }
        CppAD::vector<Type> a(2);
        CppAD::vector<Type> b(2);
        CppAD::vector<Type> aa;
        CppAD::vector<Type> bb;
        Type dig_a;
        Type dig_b;
        Type trig_a;
        Type trig_b;
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              // define mu, mu' and mu''
              mu(i,j) = 0.0;
              mu_prime = 0.0;
              mu_prime2 = 0.0;
              if((y(i,j)==0)){
                  //nll -= -logspace_add(Type(0),eta(i,j)-zetacutoffnew(j,0));
                  nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(0), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(0))), eta(i,j)-zetacutoffnew(0));
                  nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(0), eta(i,j), Type(1), 1))*cQ(i,j);
              } else if((y(i,j)==1)){
                //nll -= -logspace_add(Type(0),zetacutoffnew(j,1)-eta(i,j));
                nll -= -CppAD::CondExpLe(zetacutoffnew(1)-eta(i,j), Type(18.), gllvmutils::log1plus(exp(zetacutoffnew(1)-eta(i,j))), zetacutoffnew(1)-eta(i,j));
                nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(1), eta(i,j), Type(1), 1))*cQ(i,j);
              } else{
                // if(zeta.size()>p) {
                  //nll -= log(pnorm(zetacutoffnew(j,1) - eta(i,j), Type(0), Type(1)) - pnorm(zetacutoffnew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
                  //nll -= logspace_sub(-logspace_add(Type(0),zetacutoffnew(j,0)-eta(i,j)), -logspace_add(Type(0),zetacutoffnew(j,1)-eta(i,j)));
                  nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(0), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(0))), eta(i,j)-zetacutoffnew(0));
                  nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(1), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(1))), eta(i,j)-zetacutoffnew(1));
                  nll -= eta(i,j) - zetacutoffnew(0);
                  nll -= CppAD::CondExpLe(zetacutoffnew(1)-zetacutoffnew(0), log(Type(2.)), log(-gllvmutils::expminus1(zetacutoffnew(0)-zetacutoffnew(1))),  gllvmutils::log1plus(-exp(zetacutoffnew(0)-zetacutoffnew(1))));
                  nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(0), eta(i,j), Type(1), 1))*cQ(i,j); 
                  nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(1), eta(i,j), Type(1), 1))*cQ(i,j);
                // } else { //Model without upper bound, not implemented in R side
                //   nll -= eta(i,j) - zetacutoffnew(0); 
                //   nll -= -CppAD::CondExpLe(eta(i,j)-zetacutoffnew(0), Type(18.), gllvmutils::log1plus(exp(eta(i,j)-zetacutoffnew(0))), eta(i,j)-zetacutoffnew(0));
                //   nll -= -gllvmutils::mfexp(dlogis(zetacutoffnew(0), eta(i,j), Type(1), 1))*cQ(i,j);
                // }
                CppAD::vector<Type> z(4);
                z[0] = eta(i,j);
                z[1] = 0;
                z[2] = 1/(1+exp(-z[0]));
                z[3] = exp(z[0])/(exp(z[0])+1);
            
                mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
                mu_prime = mu(i,j) * (1-mu(i,j));
                mu_prime2 = mu_prime * (1-2*mu(i,j));
                
                a[0] = mu(i,j)*iphi(j);
                a[1] = 1;
                b[0] = (1-mu(i,j))*iphi(j);
                b[1] = 1;
                aa = a;
                bb = b;
                aa[1] = 2;
                bb[1] = 2;
                dig_a = Type(atomic::D_lgamma(a)[0]);
                dig_b = Type(atomic::D_lgamma(b)[0]);
                trig_a = Type(atomic::D_lgamma(aa)[0]);
                trig_b = Type(atomic::D_lgamma(bb)[0]);
                
                nll -= dbeta(y(i,j), Type(a[0]), Type(b[0]), 1);
                nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
                nll -= iphi(j) * mu_prime2 * logit(y(i,j)) * cQ(i,j);
              }
              
            }
          // }
        }
      }
      break;
    }
    
    case ZIB: { // ZIB family 13 VA
      Type iphij = iphi(j)/(1+iphi(j));
      Type pVA;
      if(method == 0 && extra(j)<1){
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0){
                nll -= log(1-iphij);
                Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
                nll -= (y(i,j)-Ntrials(i,j)*0.5)*eta(i,j) - Ntrials(i,j)*logspace_add(wij, -wij);

                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else{
                Type LL = 0;
                Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
                LL += (-Ntrials(i,j)*0.5)*eta(i,j) - Ntrials(i,j)*logspace_add(wij, -wij);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  LL += lgamma(Ntrials(i,j)+1.) - lgamma(Ntrials(i,j)+1.);//norm.const.
                }
                
                pVA = exp(log(-iphij+1)+LL-log((1-iphij)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }
            }
          }
        // }
      }else if(method == 0 && extra(j)==1){
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0){
                nll -= log(1-iphij);
                nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(i,j)-y(i,j));
                nll += cQ(i,j)*Ntrials(i,j);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else{
                Type LL = 0;
                LL += log(1-mu(i,j))*Ntrials(i,j);
                LL -= cQ(i,j)*Ntrials(i,j);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  LL += lgamma(Ntrials(i,j)+1.) - lgamma(Ntrials(i,j)+1.);//norm.const.
                }
                
                pVA = exp(log(-iphij+1)+LL-log((1-iphij)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }
            }
          }
        // }
      }else if(method == 0 && extra(j)==2){
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            mu(i,j) = exp(eta(i,j) + cQ(i,j));
  
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0){
                nll -= log(1-iphij);
                nll -= y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else{
                Type LL = 0;
                LL += y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  LL += lgamma(Ntrials(i,j)+1.) - lgamma(Ntrials(i,j)+1.);//norm.const.
                }
                
                pVA = exp(log(-iphij+1)+LL-log((1-iphij)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }
            }
          }
        // }
      }
      break;
    }
    
    case ZNIB: { // ZNIB family 14 (VA)
      Type iphij = exp(lg_phi(j))/(1+exp(lg_phi(j)) + exp(lg_phiZINB(j)));
      // vector<Type> iphi2 = exp(lg_phiZINB)/(1+exp(lg_phi) + exp(lg_phiZINB));
      // vector<Type> iphi3 = iphi+iphi2;
      Type iphi2 = exp(lg_phiZINB(j))/(1+exp(lg_phi(j)) + exp(lg_phiZINB(j)));
      Type iphi3 = iphij+iphi2;
      Type pVA;
      Type pVA2;
      if(method == 0 && extra(j)<1){
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0 && y(i,j)< Ntrials(i,j)){
                nll -= log(1-iphi3);
                Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
                nll -= (y(i,j)-Ntrials(i,j)*0.5)*eta(i,j) - Ntrials(i,j)*logspace_add(wij, -wij);

                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else if(y(i,j)==0){
                Type LL = 0;
                Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
                LL += (-Ntrials(i,j)*0.5)*eta(i,j) - Ntrials(i,j)*logspace_add(wij, -wij);

                pVA = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }else if(y(i,j) == Ntrials(i,j)){
                Type LL = 0;
                Type wij = 0.5*sqrt(eta(i,j)*eta(i,j) + 2*cQ(i,j));
                LL += (y(i,j)-Ntrials(i,j)*0.5)*eta(i,j) - Ntrials(i,j)*logspace_add(wij, -wij);
                
                pVA2 = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphi2));
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(1), pVA2-Type(1e-12), pVA2));//check if pVA is on the boundary
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(0), pVA2+Type(1e-12), pVA2));//check if pVA is on the boundary
                nll -= log(iphi2)-log(1-pVA2);
              }
            }
          }
        // }
      }else if(method == 0 && extra(j)==1){
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0 && y(i,j)< Ntrials(i,j)){
                nll -= log(1-iphi3);
                nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(i,j)-y(i,j));
                nll += cQ(i,j)*Ntrials(i,j);

                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else if(y(i,j)==0){
                Type LL = 0;
                LL += log(1-mu(i,j))*Ntrials(i,j);
                LL -= cQ(i,j)*Ntrials(i,j);
                
                
                pVA = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }else if(y(i,j) == Ntrials(i,j)){
                Type LL = 0;
                LL += y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(i,j)-y(i,j));
                LL -= cQ(i,j)*Ntrials(i,j);
                
                pVA2 = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphi2));
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(1), pVA2-Type(1e-12), pVA2));//check if pVA is on the boundary
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(0), pVA2+Type(1e-12), pVA2));//check if pVA is on the boundary
                nll -= log(iphi2)-log(1-pVA2);
              }
            }
          }
        // }
      }else if(method == 0 && extra(j)==2){
          for (int i=0; i<n; i++) {
            mu(i,j) = exp(eta(i,j) + cQ(i,j));

            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0 && y(i,j)< Ntrials(i,j)){
                nll -= log(1-iphi3);
                nll -= y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);
                
                if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                  nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
                }
              }else if(y(i,j)==0){
                Type LL = 0;
                LL += y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);

                
                pVA = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphij));
                pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
                pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
                nll -= log(iphij)-log(1-pVA);
              }else if(y(i,j) == Ntrials(i,j)){
                Type LL = 0;
                LL += y(i,j)*gllvmutils::log1plus(-exp(-mu(i,j)*exp(-cQ(i,j))))-(Ntrials(i,j)-y(i,j))*mu(i,j) + mu(i,j)*(exp(-cQ(i,j))-1);
                
                pVA2 = exp(log(1-iphi3)+LL-log((1-iphi3)*exp(LL)+iphi2));
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(1), pVA2-Type(1e-12), pVA2));//check if pVA is on the boundary
                pVA2 = Type(CppAD::CondExpEq(pVA2, Type(0), pVA2+Type(1e-12), pVA2));//check if pVA is on the boundary
                nll -= log(iphi2)-log(1-pVA2);
              }
            }
          }
      }
      break;
    }
    
      default: {
        // Error message for non-available family
        error("%s", ("Unsupported family at column " + std::to_string(j) +
          std::string(": ") + std::to_string(static_cast<int>(family(j)))).c_str());
      }
    } // switch
  } // for j
