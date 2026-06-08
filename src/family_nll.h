// family_nll.h — family-specific NLL contributions for gllvm VA/EVA
//
// Include this file inside the body of a TMB objective function, after the
// following outer-scope variables have been declared and populated:
//
//   Scalars / integers:
//     int n          — number of sites
//     int p          — total number of response columns (incl. betaH)
//     int truep      — p - p_betaH (columns in the family switch)
//     int method     — 0=VA, 2=EVA
//     int zetastruc  — 0=common cutoffs, 1=species-specific cutoffs
//
//   Vectors (length p or p-related):
//     vector<Type> family      — family code per column
//     vector<Type> extra       — link/variant flag per column
//     vector<Type> iphi        — exp(lg_phi), dispersion
//     vector<Type> lg_phi      — log(dispersion)
//     vector<Type> lg_phiZINB  — log(dispersion) for ZINB
//     vector<Type> zeta        — ordinal cutoff raw parameters
//     Type          ePower     — Tweedie power (logit-scale)
//
//   Matrices (n x p):
//     matrix<Type> y     — response
//     matrix<Type> eta   — linear predictor
//     matrix<Type> cQ    — VA correction term (½ * tr(Hess * Var))
//     matrix<Type> mu    — scratch/output (written for probit etc.)
//     matrix<Type> Ntrials — binomial trials
//
//   Scalar accumulator (mutated here):
//     parallel_accumulator<Type> nll  (or Type nll)
//
// Variables declared internally (do NOT pre-declare in caller):
//     int  idx, idxj, K, Kj, ymax, ymaxj
//     bool has12
    #include "family_va_nll.h"

    
  }else{
    // method = "LA"
    
    using namespace density;
    if(random(2)>0){
      // REPORT(Sigmab_lv); //!!!!
      if((randomB>0) && (csb_lv.cols()<2)){
        //randomB == "lv" without correlation
        MVNORM_t<Type> mvnorm(Sigmab_lv(0));
        for (int klv=0; klv<Klv; klv++) {
          nll += mvnorm(b_lv.row(klv));
        }
      }else if((randomB>0) && (csb_lv.cols()==2)){
        //randomB == "lv" with correlation
        matrix<Type> SigmaB_lvC = Sigmab_lv(0)*Sigmab_lv(0).transpose();//correlation matrix
        for (int q=0; q<(num_lv_c+num_RR); q++) {
          matrix<Type>SigmaB_lv = SigmaB_lvC*exp(sigmab_lv(q))*exp(sigmab_lv(q));
          nll += MVNORM(SigmaB_lv)(b_lv.col(q));
        }
      }else if((randomB<1) && (csb_lv.cols()<2)){
        for (int q=0; q<(num_lv_c+num_RR); q++) {
          nll += MVNORM(Sigmab_lv(q))(b_lv.col(q));
        }
      }else if((randomB<1) && (csb_lv.cols()>1)){
        //randomB = "P" with correlation
        vector<Type>sigma2(num_lv_c+num_RR);
        sigma2.fill(1.0);
        sigma2.tail(num_lv_c+num_RR-1) = pow(exp(sigmab_lv.segment(x_lv.cols(), num_lv_c+num_RR-1)), 2);
        
        matrix<Type>SigmaB_lv = Sigmab_lv(0)*Sigmab_lv(0).transpose();
        for (int q=0; q<(num_lv_c+num_RR); q++) {
          SigmaB_lv *= sigma2(q);
          nll += MVNORM(SigmaB_lv)(b_lv.col(q));
          SigmaB_lv /= sigma2(q);
        }
      }
    }
    
    // REPORT(ucopy);
    // REPORT(num_corlv);
    // REPORT(nu);
    // REPORT(dr0);
    // REPORT(cstruc);
    
    // matrix<Type> etaH(n,p); 
    // etaH.setZero();
    
    //For fixed-effects RRR with and without quadratic term
    if(num_RR>0){
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      if(quadratic>0){
        matrix<Type> D_RR(num_RR,num_RR);
        D_RR.setZero();
        if(lambda2.cols()==1){
          for (int d=0; d<num_RR;d++){
            D_RR.diagonal()(d) = fabs(lambda2(d,0));
          }
          for (int i=0; i<n; i++) {
            eta.row(i).array() -=  (x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose()).value();
          }
          
        }else{
          for (int j=0; j<p;j++){
            D_RR.setZero();
            for (int d=0; d<num_RR;d++){
              D_RR.diagonal()(d) = fabs(lambda2(d,j));
            }
            
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
        }
        
      }
    }
    
    // Laplace approximation
    
    // add offset to lin. predictor 
    if(offset.rows()==n){
      eta += offset;
    }
    // if(r0f.size() == n && (random(0)==0) && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){
    //   eta += r0f.replicate(1,p);
    if(xr.rows()==n){
      eta += (xr*r0f).replicate(1,p);
    }
    
    if((random(1)>0) || (random(3)>0)){
      if(random(1)>0){
        // random slopes in TMBtrait.R
        eta += xb*Br;  
      }else if(random(3)>0){
        //random slopes in gllvm.TMB
        eta += xb*(Br.colwise()+B.col(0));
      }
      
      matrix <Type> Spr(xb.cols(),xb.cols());
      matrix <Type> SprI(xb.cols(),xb.cols());
      Type logdetSpr = 0;
      
      int l = xb.cols();
      if(random(1)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB);
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL(l,l);
        SprL.fill(0.0);
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaij(i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        matrix <Type> SprIL(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      if(random(3)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB.segment(0,xb.cols()));
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL;
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaB(xb.cols()+i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        matrix <Type> SprIL(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      if(colMatBlocksI(0)(0,0)==p){
        Type logdetColCorMat = 0;
        vector <Type> rhoSP(1);
        
        if(random(1)>0 && sigmaB.size()>xb.cols()){
          rhoSP.resize(sigmaB.size()-xb.cols());
          
          rhoSP.fill(1.0);
            if(sigmaB.size()>xb.cols()){//traitTMB has correlation parameters in sigmaij. Ultimately, these should also go via sigmaij in gllvm.TMB.
              rhoSP = exp(-exp(sigmaB.segment(xb.cols(),sigmaB.size()-xb.cols())));
              if(nncolMat.rows()<p){
                //need to cap this on the lower end for numerical stability
                for(int re=0; re<rhoSP.size(); re++){
                  rhoSP(re) = CppAD::CondExpLt(rhoSP(re), Type(1e-12), Type(1e-12), rhoSP(re));
                }
              }
            }

        }else if(random(3)>0){
          rhoSP.resize(sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1));
          
            rhoSP.fill(1.0);
            if(sigmaB.size()>(xb.cols()+cs.rows()*(cs.cols()>1))){//unLike traitTMB, gllvm.TMB has correlation parameters also in sigmaB. Ultimately, these should also go via sigmaij in gllvm.TMB.
              rhoSP = exp(-exp(sigmaB.segment(xb.cols()+cs.rows()*(cs.cols()>1),sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1))));
              if(nncolMat.rows()<p){
                //need to cap this on the lower end for numerical stability
                for(int re=0; re<rhoSP.size(); re++){
                  rhoSP(re) = CppAD::CondExpLt(rhoSP(re), Type(1e-12), Type(1e-12), rhoSP(re));
                }
              }
            }
        }
        
        
        int sp = 0;
        
        if(nncolMat.rows()<p && rhoSP.size()==1){
          //only go here if rhoSP.size()==1. Other case we need a cholesky.
          logdetColCorMat = colMatBlocksI(0).col(1).segment(1,colMatBlocksI.size()-1).sum();
          for(int cb=1; cb<colMatBlocksI.size(); cb++){
            //efficiently update inverse and determinant using rank 1 updates
            matrix<Type> colCorMatI(colMatBlocksI(cb).cols(), colMatBlocksI(cb).cols());
            gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb), logdetColCorMat, rhoSP(0));
            //formulate MVNORM manually because we have all the components
            nll += 0.5*(Br.middleCols(sp, colCorMatI.cols())*colCorMatI*Br.middleCols(sp, colCorMatI.cols()).transpose()*SprI).trace();
            sp += colCorMatI.cols();
          }
          //add cheap normalizing constant
          nll -= 0.5*(p*xb.cols()-p*logdetSpr-xb.cols()*logdetColCorMat);
        }else if(nncolMat.rows()==p){
          if(rhoSP.size()==1){//p(block) sized matrix
            int sp = 0;
            for(int cb=1; cb<colMatBlocksI.size(); cb++){
              Eigen::SparseMatrix<Type> AL(colMatBlocksI(cb).cols(),colMatBlocksI(cb).cols());
              gllvmutils::nngp(AL, colMatBlocksI(cb), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb).cols()));
              nll += 0.5*(Br.middleCols(sp, AL.cols())*AL*AL.transpose()*Br.middleCols(sp, AL.cols()).transpose()*SprI).trace();
              sp += colMatBlocksI(cb).cols();
            }
            //add cheap normalizing constant
            nll -= 0.5*(p*xb.cols()-p*logdetSpr-xb.cols()*logdetColCorMat);
            
          }else{//now colCorMatI is updated for every covariate
            int sp = 0;
            for(int cb=1; cb<colMatBlocksI.size(); cb++){
              Eigen::SparseMatrix<Type> AL(colMatBlocksI(cb).cols(), colMatBlocksI(cb).cols());
              for(int re=0; re<rhoSP.size(); re++){
                gllvmutils::nngp(AL, colMatBlocksI(cb), logdetColCorMat, rhoSP(re), nncolMat.middleCols(sp, colMatBlocksI(cb).cols()));
                Br.row(re).middleCols(sp,colMatBlocksI(cb).cols()) *= AL;
              }
              sp += colMatBlocksI(cb).cols();
            }
            nll += 0.5*(Br*Br.transpose()*SprI).trace();
            nll -= 0.5*(p*xb.cols()-p*logdetSpr-logdetColCorMat);
          }
          
        }
      }else{
        //independence across species
        MVNORM_t<Type> mvnorm(Spr);
        for (int j=0; j<p;j++){
          nll += mvnorm(Br.col(j));
        }
      }
    }
    
    //latent variables
    if(nlvr>0){
      if(num_corlv==0){
        for (int i=0; i<n; i++) {
          for(int q=0; q<u.cols(); q++){
            nll -= dnorm(u(i,q), Type(0), Type(1), true);
          }
        }
      }
      //variances of LVs
      u *= Delta;
      if(num_lv_c>0){
        matrix<Type> b_lv2(x_lv.cols(),nlvr);
        
        b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        eta += x_lv*b_lv2*newlam;
      }
      // add LV term to lin. predictor 
      lam += u*newlam;
      eta += lam;
      // if(family(j)==10){
      //   // etaH += lam;
      //   etaH += ucopy*thetaH;
      // }
    }
    
    
    // Row/Site effects
    if((random(0)>0)){
      vector<Type> sigma = exp(log_sigma);
      eta += (dr0*r0r).replicate(1,p);//matrix<Type>(Eigen::MatrixXd::Ones(1,p));
      
      int dccounter = 0; // tracking used dc entries
      int sigmacounter = 0; // tracking used sigma entries
      int ucount = 0;
      int propcount = 0;
      for(int re=0; re<trmsize.cols();re++){
        
        if(cstruc(re)<0 || cstruc(re)>5){
          matrix<Type> Sr(trmsize(0,re), trmsize(0,re));

          matrix<Type> sds = Eigen::MatrixXd::Zero(trmsize(0,re),trmsize(0,re));
          sds.diagonal() =  sigma.segment(sigmacounter, trmsize(0,re));
          sigmacounter += trmsize(0,re);
          
          vector<Type>sigmaRij((trmsize(0,re)*trmsize(0,re)-trmsize(0,re))/2);
          sigmaRij.fill(0.0);
          //covariances of random effects
          matrix<Type> SrL(trmsize(0,re),trmsize(0,re));
          SrL.fill(0.0);
          if(csR.cols()>1){
            //need a vector with covariances and zeros in the right places
            for(int i=0; i<sigmaRij.size(); i++){
              sigmaRij((csR(ucount,0) - 1) * (csR(ucount,0) - 2) / 2 + csR(ucount,1)-1) = sigmaijr(ucount);
              ucount++;
            }
            SrL = sds*gllvmutils::constructL(sigmaRij);
          }else{
            SrL = sds;
          }
          Sr = SrL*SrL.transpose();
          
      if(cstruc(re)<0){  
        MVNORM_t<Type> MVNSr(Sr);
        for (int q=0; q<trmsize(1,re); q++){//loop over blocks
          if(re==0){
            vector<Type> r0s = r0r.col(0).segment(trmsize(0,re)*q,trmsize(0,re));
            nll += MVNSr(r0s);
          }else{
            vector<Type> r0s = r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum()+trmsize(0,re)*q,trmsize(0,re));
            nll += MVNSr(r0s);
          }
        }
      }else if(cstruc(re) > 5){
        matrix <Type> invSr(trmsize(0,re),trmsize(0,re));invSr.setZero();
        matrix <Type> Ir = Eigen::MatrixXd::Identity(SrL.cols(),SrL.cols());
        matrix <Type> SrIL(SrL.cols(),SrL.cols());
        SrIL = SrL.template triangularView<Eigen::Lower>().solve(Ir);
        SrIL = SrIL.transpose()*SrIL;
        invSr=SrIL*SrIL.transpose();
        
        Type logdetSr = 2*SrL.diagonal().array().log().sum();
        matrix<Type>invMat(trmsize(1,re), trmsize(1,re));
        
       if(cstruc(re)>6){
         // here we need to calculate the inverse of our second covariance matrix
         // as we have a kronecker product, and variances are in SrL, the matrices below are correlation matrices.
         // this keeps the number of constraints similar to the proptoustruc case
         matrix<Type>Sr(trmsize(1,re), trmsize(1,re));
         Sr.setZero();
         
         if(cstruc(re) == 7){ // corAR1
           Sr = gllvm::corAR1(Type(1), log_sigma(sigmacounter), trmsize(1,re));
           sigmacounter+= 1;
         }else if(cstruc(re) == 9){ // corCS
           Sr = gllvm::corCS(Type(1), log_sigma(sigmacounter), trmsize(1,re));
           sigmacounter += 1;
         }else if((cstruc(re) == 8) || (cstruc(re) == 10)){ // corMatern, corExp
           // Distance matrix calculated from the coordinates for rows
           matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
           matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
           DiSc.setZero();
           DiSc.diagonal().array() += 1/sigma(sigmacounter);
           sigmacounter++;
           dc_scaled = dc(dccounter)*DiSc;
           if(cstruc(re) == 8){ // corExp
             Sr = gllvm::corExp(Type(1), Type(0), trmsize(1,re), dc_scaled);
           } else if(cstruc(re) == 10) { // corMatern
             Sr = gllvm::corMatern(Type(1), Type(1), sigma(sigmacounter), trmsize(1,re), dc_scaled);
             sigmacounter += 1;
           }
           dccounter++;
         }
         
         //TMB's matinvpd function: inverse of matrix with logdet for free
         CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
         logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*res[0];
         invMat = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
       }else if(cstruc(re)==6){
         // here we have a known inverse
         invMat = proptoMats(propcount)(0);
         logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*proptoMats(propcount)(1)(0); //logdet kronecker
         
         propcount ++;
       }
        
        if(re==0){
          Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(0, trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
          nll -=  -0.5*(bm*invMat*bm.transpose()*invSr).trace();
        }else{
          Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
          nll -=  -0.5*(bm*invMat*bm.transpose()*invSr).trace();
        }
        
        // determinants of each block of the covariance matrix
        nll -= -0.5*trmsize(0,re)*trmsize(1,re)*log(2*M_PI)-0.5*logdetSr;
        
      }
        }else{
          matrix<Type> Sr(trmsize(1,re),trmsize(1,re));Sr.setZero();
          
        if(cstruc(re) < 5){
        if(cstruc(re) == 0){
          Sr.diagonal().array() = pow(sigma(sigmacounter), 2);
          sigmacounter++;
        }else if(cstruc(re) == 1){ // corAR1
          Sr = gllvm::corAR1(sigma(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
          sigmacounter+=2;
        }else if(cstruc(re) == 3){ // corCS
          Sr = gllvm::corCS(sigma(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
          sigmacounter += 2;
        }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
          // Distance matrix calculated from the coordinates for rows
          matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
          matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
          DiSc.setZero();
          DiSc.diagonal().array() += 1/sigma(sigmacounter);
          sigmacounter++;
          dc_scaled = dc(dccounter)*DiSc;
          if(cstruc(re)==2){ // corExp
            Sr = gllvm::corExp(sigma(sigmacounter), Type(0), trmsize(1,re), dc_scaled);
            sigmacounter++;
          } else if(cstruc(re)==4) { // corMatern
            Sr = gllvm::corMatern(sigma(sigmacounter), Type(1), sigma(sigmacounter+1), trmsize(1,re), dc_scaled);
            sigmacounter += 2;
          }
          dccounter++;
        }
        
        if(cstruc(re)==0){
          //independence of REs
          if(re==0){
          vector<Type> r0s = r0r.col(0).segment(0,trmsize(1,re));
            
          for(int ir=0; ir<r0s.size(); ir++){
            nll -= dnorm(r0s(ir), Type(0), sigma(sigmacounter-1), true);
          }
          }else{
            vector<Type> r0s = r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re));

            for(int ir=0; ir<r0s.size(); ir++){
              nll -= dnorm(r0s(ir), Type(0), sigma(sigmacounter-1), true);
            }
          }
        }else{
          if(re==0){
            vector<Type> r0s = r0r.col(0).segment(0,trmsize(1,re));
            nll += MVNORM(Sr)(r0s);
          }else{
            vector<Type> r0s = r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re));
            nll += MVNORM(Sr)(r0s);
          }
        }
        }else{
          matrix<Type> invSr(trmsize(1,re), trmsize(1,re));
          invSr = pow(sigma(sigmacounter), -2)*proptoMats(propcount)(0);
          Type logdetSr = proptoMats(propcount)(1)(0) + 2*proptoMats(propcount)(0).cols()*log_sigma(sigmacounter);
          sigmacounter++;
          propcount++;
          if(re==0){
          nll -= -Type(trmsize(1,re))/2*log(2*M_PI) - 0.5*logdetSr -0.5*r0r.col(0).segment(0,trmsize(1,re)).transpose()*invSr*r0r.col(0).segment(0,trmsize(1,re));
          }else{
          nll -= -Type(trmsize(1,re))/2*log(2*M_PI) - 0.5*logdetSr -0.5*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)).transpose()*invSr*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re));  
          }
        }
        }
        
      }
    }
    
    // Correlated LVs
    if(num_corlv>0) {
      int i;
      // if(ucopy.rows() == nu){
      if(cw == 0){
          // eta += (dLV*ucopy)*newlamCor;
        // if(family(j)==10){ // betaH
        //   etaH += (dLV*ucopy)*thetaH;
        // }
        
        // group specific lvs
        if(cstruclv(0)==0){// no covariance
          matrix<Type> Slv(num_corlv,num_corlv);
          Slv.setZero();
          Slv.diagonal().fill(1.0);
          MVNORM_t<Type> mvnorm(Slv);
          for (int i=0; i<nu; i++) {
            nll += mvnorm(ucopy.row(i));
          }
          // REPORT(Slv);
        } else {
          
          matrix<Type> Slv(nu,nu);
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated between groups
            Slv.setZero();
            
            if(cstruclv(0)==1){// AR1 covariance
              Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
            } else if(cstruclv(0)==3) {// Compound Symm  if(cstruclv==3)
              Slv = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
            } else {
              DiSc_lv.setZero();
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,0));
              }
              dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv(0)==2){// exp decaying
                Slv = gllvm::corExp(Type(1), Type(0), nu, dc_scaled_lv);
                // Slv = gllvm::corExp(Type(1), (rho_lvc(q,0)), nu, DistM);
              } else if(cstruclv(0)==4) {// matern
                Slv = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,rho_lvc.cols()-1)), nu, dc_scaled_lv);
              }
            }
            
            MVNORM_t<Type> mvnormS1(Slv);
            nll += mvnormS1(ucopy.col(q));
          }
          // REPORT(Slv);
        }
      } else {
        int it_ind = 0;
        matrix<Type> Slv;
        for (i=0; i<times.row(0).size(); i++) {
          Slv.resize(times(0,i),times(0,i));
          
          // eta += ucopy*newlamCor;
          // if(family(j)==10){// betaH
          //   etaH += ucopy*thetaH;
          // }
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            Slv.setZero();
            // Define covariance matrix
            int ics =0;
            if(cstruclv.size() >= nu) ics =i;
            if(cstruclv(ics)==1){// AR1 covariance
              Slv = gllvm::corAR1(Type(1), rho_lvc(q,i), times(0,i));
            } else if(cstruclv(ics)==3) {// Compound Symm  if(cstruclv==3)
              Slv = gllvm::corCS(Type(1), rho_lvc(q,i), times(0,i));
            } else {
              DiSc_lv.setZero();
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,i));
                // DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled_lv = dc_lv.block(it_ind,0,times(0,i),dc_lv.cols())*DiSc_lv;
              // dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv(ics)==2){// exp decaying
                Slv = gllvm::corExp(Type(1), Type(0), times(0,i), dc_scaled_lv);
              } else if(cstruclv(ics)==4) {// matern
                Slv = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,rho_lvc.cols()-1)), times(0,i), dc_scaled_lv);
              }
            }
            
            MVNORM_t<Type> mvnormS2(Slv);
            
            nll += mvnormS2(ucopy.block(it_ind,q,times(0,i),1));
            // for (i=0; i<nu; i++) {
            //   nll += mvnormS2(ucopy.block(i*times,q,times,1));
            // }
            
          }
          it_ind += times(0,i);
        }
        // REPORT(Slv);
      }
      // REPORT(nu);
    }
    
    if(model<1){
      // gllvm.TMB.R
      // if(family(j)==10){
      //   etaH += x*bH;
      // }
      eta += x*b;
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          mu(i,j) = exp(eta(i,j));
        }
      }
      
    } else {
      // Fourth corner model, TMBtrait.R
      // if(family(j)==10){
      //   matrix<Type> eta1h=x*bH;
      //   eta1h.resize(n, p);
      //   etaH += eta1h;
      // }
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)*extra(p)+eta1(m,0);
          m++;
          mu(i,j) = exp(eta(i,j));
        }
      }
    }
    
    
    int idx = 0; // initialize indexing for zeta
    
    //likelihood model with the log link function
    for (int j=0; j<truep; j++){
      
      switch (family(j)) {

      case POISSON: { //poisson family 0
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= dpois(y(i,j), exp(eta(i,j)), true);
        }
        break;
      }
      
      case NEG_BINOMIAL: {//negative.binomial family 1
        if(extra(j)==0){
          //nb2
          if((num_RR>0) && (nlvr == 0) && (random(2)<1)){
            //use dnbinom_robust in this case - below code does not function well
            //for constrained ordination without any random-effects
              for (int i=0; i<n; i++) {
                if(!gllvmutils::isNA(y(i,j)))nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
              }
          }else{
              for (int i=0; i<n; i++) {
                if(!gllvmutils::isNA(y(i,j)))nll -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
              }
          }
        }else if(extra(j)==1){
          //nb1
            for (int i=0; i<n; i++) {
              if(!gllvmutils::isNA(y(i,j)))nll -= dnbinom_robust(y(i,j), eta(i,j), eta(i,j) - lg_phi(j), 1);
            }
        }
        break;
      } 
      
      case BINOMIAL: {//binomial family 2
          for (int i=0; i<n; i++) {
            if(extra(j)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else if(extra(j)==1){mu(i,j) = pnorm(eta(i,j));
            }else if(extra(j)==2)mu(i,j) = 1-exp(-exp(eta(i,j)));
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            if(!gllvmutils::isNA(y(i,j))){
              nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(i,j)-y(i,j));
              if(Ntrials(i,j)>1 && (Ntrials(i,j)>y(i,j))){
                nll -= lgamma(Ntrials(i,j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(i,j)-y(i,j)+1.);//norm.const.
              }
            }
          }
        break;
      } 
      
      case GAUSSIAN: {//gaussian family 3
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))) nll -= dnorm(y(i,j), eta(i,j), iphi(j), true); 
        }
        break;
      } 
      
      case GAMMA: {//gamma family 4
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= dgamma(y(i,j), iphi(j), exp(eta(i,j))/iphi(j), true); 
        }
        break;
      } 
      
      case TWEEDIE: {//tweedie family 5
        Type ePower1 = invlogit(ePower) + Type(1);
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))) nll -= dtweedie(y(i,j), exp(eta(i,j)),iphi(j),ePower1, true); 
        }
        break;
      } 
      
      case ZIP: {//zero-infl-poisson 6
        Type iphij=iphi(j)/(1+iphi(j));
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j)))nll -= dzipois(y(i,j), exp(eta(i,j)),iphij, true); 
          }
          break;
      } 
      
      case ORDINAL: { //ordinal family 7
        if(zetastruc == 1){//ordinal, only here for models without random-effects
          int ymax =  CppAD::Integer(y.maxCoeff());
          int K = ymax - 1;
          
          // matrix <Type> zetanew(p,K);
          vector <Type> zetanew(K);
          zetanew.setZero();
          
          // int idx = 0; // indexing moved before for j
          // for(int j=0; j<p; j++){
            int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
            int Kj = ymaxj - 1;
            if(Kj>1){
              for(int k=0; k<(Kj-1); k++){
                zetanew(k+1) = zeta.segment(idx,k+1).array().exp().sum();
              }
            }
            idx += Kj-1;
          // }
          
          if(extra(j)==0){
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                //minimum category
                if(y(i,j)==1){
                  nll -= log(invlogit(zetanew(0) - eta(i,j)));
                }else if(y(i,j)==ymaxj){
                  //maximum category
                  int idxj = ymaxj-2;
                  nll -= log(1 - invlogit(zetanew(idxj) - eta(i,j)));
                }else if(ymaxj>2){
                  for (int l=2; l<ymaxj; l++) {
                    if((y(i,j)==l) && (l != ymaxj)){
                      nll -= log(invlogit(zetanew(l-1)-eta(i,j))-invlogit(zetanew(l-2)-eta(i,j)));
                    }
                  }
                }
              // }
            }
          }else if(extra(j)==1){
          for (int i=0; i<n; i++) {
            // for(int j=0; j<p; j++){
              int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
              //minimum category
              if(y(i,j)==1){
                nll -= log(pnorm(zetanew(0) - eta(i,j), Type(0), Type(1)));
              }else if(y(i,j)==ymaxj){
                //maximum category
                int idxj = ymaxj-2;
                nll -= log(1 - pnorm(zetanew(idxj) - eta(i,j), Type(0), Type(1)));
              }else if(ymaxj>2){
                for (int l=2; l<ymaxj; l++) {
                  if((y(i,j)==l) && (l != ymaxj)){
                    nll -= log(pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1)));
                  }
                }
              }
            // }
          }
          }
        } else if(zetastruc==0){
          int ymax =  CppAD::Integer(y.col(j).maxCoeff());
          int K = ymax - 1;
          
          vector <Type> zetanew(K);
          zetanew.setZero();
          for(int k=0; k<(K-1); k++){
            zetanew(k+1) = zeta.head(k+1).array().exp().sum();
          }

          if(extra(j)==0){
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                if(!gllvmutils::isNA(y(i,j))){
                  //minimum category
                  if(y(i,j)==1){
                    nll -= log(invlogit(zetanew(0) - eta(i,j)));
                  }else if(y(i,j)==ymax){
                    //maximum category
                    int idxj = ymax-2;
                    nll -= log(1 - invlogit(zetanew(idxj) - eta(i,j)));
                  }else if(ymax>2){
                    for (int l=2; l<ymax; l++) {
                      if((y(i,j)==l) && (l != ymax)){
                        nll -= log(invlogit(zetanew(l-1)-eta(i,j))-invlogit(zetanew(l-2)-eta(i,j)));
                      }
                    }
                  }
                }
              // }
            }
          }else if(extra(j)==1){
            for (int i=0; i<n; i++) {
              // for(int j=0; j<p; j++){
                if(!gllvmutils::isNA(y(i,j))){
                  //minimum category
                  if(y(i,j)==1){
                    nll -= log(pnorm(zetanew(0) - eta(i,j), Type(0), Type(1)));
                  }else if(y(i,j)==ymax){
                    //maximum category
                    int idxj = ymax-2;
                    nll -= log(1 - pnorm(zetanew(idxj) - eta(i,j), Type(0), Type(1)));
                  }else if(ymax>2){
                    for (int l=2; l<ymax; l++) {
                      if((y(i,j)==l) && (l != ymax)){
                        nll -= log(pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1)));
                      }
                    }
                  }
                }
              // }
            }
          }
        }
        break;
      }
      
      case EXPONENTIAL: {// exponential family 8
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p;j++){
            if(!gllvmutils::isNA(y(i,j)))nll -= dexp(y(i,j), exp(-eta(i,j)), true);  // (-eta(i,j) - exp(-eta(i,j))*y(i,j) );
          // }
        }
        break;
      } 
      
      case BETA: {// beta family 9
        for (int i=0; i<n; i++) {
          // for (int j=0; j<p;j++){
            if(extra(j)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            if(!gllvmutils::isNA(y(i,j)))nll -= dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
          // }
        }
        break;
      } 
      
      case BETA_HURDLE: {// beta hurdle family 10
        for (int i=0; i<n; i++) {
          // for (int j=0; j<truep; j++){
            if(extra(j)<1) {
              // etaH(i,j) = exp(etaH(i,j))/(exp(etaH(i,j))+1);
              mu(i,j) = mu(i,j)/(mu(i,j)+1);
              mu(i,truep+j) = mu(i,truep+j)/(mu(i,truep+j)+1);
            } else {
              // etaH(i,j) = pnorm(etaH(i,j));
              mu(i,j) = pnorm(eta(i,j));
              mu(i,truep+j) = pnorm(eta(i,truep+j));
            }
            if(!gllvmutils::isNA(y(i,j))){
              if (y(i,j) == 0) {
                // nll -= log(1-mu(i,j));
                nll -= log(1-mu(i,truep+j));
              } else{
                // nll -= log(mu(i,j)) + dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
                nll -= log(mu(i,truep +j)) + dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
              }
            }
          // }
        }
        // REPORT(mu);
        // REPORT(etaH);
        break;
      } 
      
      case ZINB: {//zero-infl-NB 11
        Type iphij=iphi(j)/(1+iphi(j));
        // vector<Type> iphiZINB = exp(lg_phiZINB);
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0){
                nll -= log(1-iphij) + dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phiZINB(j), 1);
              }else{
                nll -= log(iphij + (Type(1)-iphij)*dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phiZINB(j), 0)); 
              }
            }
          }
        // }
        break;
      }
      
      case ZIB: { // Zero-Inflated-Binomial, ZIB 13
        Type iphij=iphi(j)/(1+iphi(j));
          for (int i=0; i<n; i++) {
            if(extra(j)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0){
                nll -= log(1-iphij) + dbinom(y(i,j), Type(Ntrials(i,j)), mu(i,j), 1);
              }else{
                nll -= log(iphij + (Type(1)-iphij)*dbinom(y(i,j), Type(Ntrials(i,j)), mu(i,j), 0)); 
              }
            }
          }
        break;
      }
      
      case ZNIB: { // ZNIB 14
        Type iphij = exp(lg_phi(j))/(1+exp(lg_phi(j)) + exp(lg_phiZINB(j)));
        // vector<Type> iphi2 = exp(lg_phiZINB)/(1+exp(lg_phi) + exp(lg_phiZINB));
        Type iphi2 = exp(lg_phiZINB(j))/(1+exp(lg_phi(j)) + exp(lg_phiZINB(j)));
        Type iphi3 = iphij+iphi2;
        // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(extra(j)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            if(!gllvmutils::isNA(y(i,j))){
              if(y(i,j)>0 && y(i,j) < Ntrials(i,j)){
                nll -= log(1-iphi3) + dbinom(y(i,j), Type(Ntrials(i,j)), mu(i,j), 1);
              }else if(y(i,j)==0){
                nll -= log(iphij + (Type(1)-iphi3)*dbinom(y(i,j), Type(Ntrials(i,j)), mu(i,j), 0)); 
              }else if(y(i,j) == Ntrials(i,j)){
                nll -= log(iphi2 + (Type(1)-iphi3)*dbinom(y(i,j), Type(Ntrials(i,j)), mu(i,j), 0)); 
                
              }
            }
          }
        // }
      }
      
      case BETA_BINOMIAL: { // beta-binomial family 15
        for (int i=0; i<n; i++) {
          if(extra(j)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
          } else if(extra(j)==1){mu(i,j) = pnorm(eta(i,j));
          }else if(extra(j)==2) mu(i,j) = 1-exp(-exp(eta(i,j)));
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));
          if(!gllvmutils::isNA(y(i,j))){
            Type alpha = mu(i,j) * iphi(j);
            Type beta_shape = (1 - mu(i,j)) * iphi(j);
            nll -= lgamma(alpha + beta_shape) + lgamma(alpha + y(i,j)) + lgamma(beta_shape + Ntrials(i,j) - y(i,j))
                   - lgamma(alpha) - lgamma(beta_shape) - lgamma(alpha + beta_shape + Ntrials(i,j))
                   + lgamma(Ntrials(i,j) + 1.) - lgamma(y(i,j) + 1.) - lgamma(Ntrials(i,j) - y(i,j) + 1.);
          }
        }
        break;
      }

      default: {
        // Error message for non-available family
        error("%s", ("Unsupported family at column " + std::to_string(j) +
          std::string(": ") + std::to_string(static_cast<int>(family(j)))).c_str());
        break;
      }

      } // switch
    } // for j end
