// row_effects.h
// Shared row/site random effects NLL for gllvm.cpp and gllvm_HO.cpp.
// Requires in scope: random, dr0, r0r, lg_Ar, log_sigma, sigmaijr, trmsize,
//   cstruc, csR, proptoMats, dc, n, p, eta, nll.
// NOTE: The local 'sigma' variable in gllvm.cpp is renamed to 'sigma_re' here
//   to avoid conflicts with the HO model's LV-scale sigma vector.
#pragma once
    // Row/Site effects
    if((random(0)>0)){
      vector<Type> sigma_re = exp(log_sigma);
      eta += (dr0*r0r).replicate(1,p);//*matrix<Type>(Eigen::MatrixXd::Ones(1,p));
      int dccounter = 0; // tracking used dc entries
      int sigmacounter = 0; // tracking used sigma_re entries
      
      // One: build Arm, variational covariance matrix for all random effects as a list
      int sdcounter = 0;
      int covscounter = 0; //starts at # of VA scale pars
      
      //not ideal nor necessary, should eventually be replaced
      for(int i=0; i<trmsize.cols(); i++){
        if(cstruc(i)<6){
          covscounter += trmsize(0,i)*trmsize(1,i);  
        }else if(cstruc(i) >5){ //kronecker product VA
          covscounter += trmsize(0,i) + trmsize(1,i) -1;
        }
        
      }
      int VAcovs = covscounter; //to see if we are doing unstructured VA or diagonal
      int ucount = 0;
      int propcount = 0;
      for(int re=0; re<trmsize.cols();re++){
        
        //unstructured row cov
        // still need to get the parameter count here right
        if((lg_Ar.size()>VAcovs && cstruc(re)>0) || (lg_Ar.size()>VAcovs && cstruc(re)<0)){
          Type logdetSr;
          if(cstruc(re)<0 || cstruc(re) > 5){
            ///////////////////////////////////////////////////////////////////////////////////////////////
            // we go here if we have an unstructured row covariance matrix (i.e., between random effects)//
            ///////////////////////////////////////////////////////////////////////////////////////////////
            
            matrix <Type> invSr(trmsize(0,re),trmsize(0,re));invSr.setZero();
            
              matrix<Type> sds = Eigen::MatrixXd::Zero(trmsize(0,re),trmsize(0,re));
              sds.diagonal() =  sigma_re.segment(sigmacounter, trmsize(0,re));
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
              matrix <Type> Ir = Eigen::MatrixXd::Identity(SrL.cols(),SrL.cols());
              matrix <Type> SrIL(SrL.cols(),SrL.cols());
              SrIL = SrL.template triangularView<Eigen::Lower>().solve(Ir);
              SrIL = SrIL.transpose()*SrIL;
              invSr=SrIL*SrIL.transpose();
              logdetSr = 2*SrL.diagonal().array().log().sum();
              
              matrix<Type> Arm(trmsize(0,re),trmsize(0,re));
              
              if(cstruc(re)<0){
              // we go here if we have no second covariance matrix
              for (int q=0; q<trmsize(1,re); q++){//loop over blocks
                Arm.setZero();  
                for (int d=0; d<(trmsize(0,re)); d++){ // diagonals of varcov
                  Arm(d,d)=exp(lg_Ar(sdcounter));
                  sdcounter++;
                }
                
                // off-diagonals
                for (int c=0; c<(trmsize(0,re)); c++){
                  for (int r=c+1; r<(trmsize(0,re)); r++){
                    Arm(r,c)=lg_Ar(covscounter);
                    covscounter++;
                  }}

                matrix<Type> ArmMat = Arm*Arm.transpose();
              
              cQ += (0.5*(dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re))*ArmMat*dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re)).transpose()).diagonal()).replicate(1,p);
              
              if(re==0){
                nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(trmsize(0,re)*q,trmsize(0,re)).transpose()*(invSr*r0r.col(0).segment(trmsize(0,re)*q,trmsize(0,re)))).sum());  
              }else{
                nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum()+trmsize(0,re)*q,trmsize(0,re)).transpose()*(invSr*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum()+trmsize(0,re)*q,trmsize(0,re)))).sum());
              }
              
              // determinants of each block of the covariance matrix
              nll -= 0.5*(trmsize(0,re)-logdetSr);
              }
              }else if(cstruc(re) > 5){
                // we go here if we have a second covariance matrix; our RE covariance is a kronecker matrix
                matrix<Type> invMat(trmsize(1,re), trmsize(1,re));
                invMat.setZero();
                
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
                  DiSc.diagonal().array() += 1/sigma_re(sigmacounter);
                  sigmacounter++;
                  dc_scaled = dc(dccounter)*DiSc;
                  if(cstruc(re) == 8){ // corExp
                    Sr = gllvm::corExp(Type(1), Type(0), trmsize(1,re), dc_scaled);
                  } else if(cstruc(re) == 10) { // corMatern
                    Sr = gllvm::corMatern(Type(1), Type(1), sigma_re(sigmacounter), trmsize(1,re), dc_scaled);
                    sigmacounter += 1;
                  }
                  dccounter++;
                }
                
                //TMB's matinvpd function: inverse of matrix with logdet for free
                CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
                logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*res[0];
                invMat = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
                REPORT(Sr);
                }else if(cstruc(re)==6){
                // here we have a known inverse
                invMat = proptoMats(propcount)(0);
                logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*proptoMats(propcount)(1)(0); //logdet kronecker
                
                propcount ++;
                }
                
                  Arm.setZero();  
                  for (int d=0; d<(trmsize(0,re)); d++){ // diagonals of varcov
                    Arm(d,d)=exp(lg_Ar(sdcounter));
                    sdcounter++;
                  }
                  
                  // off-diagonals
                  for (int c=0; c<(trmsize(0,re)); c++){
                    for (int r=c+1; r<(trmsize(0,re)); r++){
                      Arm(r,c)=lg_Ar(covscounter);
                      covscounter++;
                    }}
                  
                  matrix<Type> ArmMat = Arm*Arm.transpose();
                  
                  matrix<Type> ArmP(trmsize(1,re), trmsize(1,re));
                  ArmP.setZero();  
                  ArmP(0,0) = 1; // identifiability
                  for (int d=1; d<(trmsize(1,re)); d++){ // diagonals of varcov
                    ArmP(d,d)=exp(lg_Ar(sdcounter));
                    sdcounter++;
                  }
                  
                  // off-diagonals
                  for (int c=0; c<(trmsize(1,re)); c++){
                    for (int r=c+1; r<(trmsize(1,re)); r++){
                      ArmP(r,c)=lg_Ar(covscounter);
                      covscounter++;
                    }}
                  
                  matrix<Type> ArmMatP = ArmP*ArmP.transpose();
                  
                  for (int q=0; q<trmsize(1,re); q++){//loop over blocks
                  cQ += ((0.5*(dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re))*ArmMat*dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re)).transpose()).diagonal())*ArmMatP.diagonal()(q)).replicate(1,p);
                  }
                  
                  if(re==0){
                    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(0, trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
                    nll -= ArmP.cols()*Arm.diagonal().array().log().sum() + Arm.cols()*ArmP.diagonal().array().log().sum() - 0.5*((invMat*ArmMatP).trace()*(invSr*ArmMat).trace()+(bm*invMat*bm.transpose()*invSr).trace());                                                   
                  }else{
                    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
                    nll -= ArmP.cols()*Arm.diagonal().array().log().sum() + Arm.cols()*ArmP.diagonal().array().log().sum() - 0.5*((invMat*ArmMatP).trace()*(invSr*ArmMat).trace()+(bm*invMat*bm.transpose()*invSr).trace());                                                   
                  }
                  
          
                  // determinants of each block of the covariance matrix
                  nll -= 0.5*(trmsize(0,re)*trmsize(1,re)-logdetSr);
                }
          }else{
            ///////////////////////////////////////////////////////////////////////////////////////////////
            /// we go here if we have an diagonal row covariance matrix (i.e., no between random effects)//
            ///////////////////////////////////////////////////////////////////////////////////////////////
            
            // here we have no 0 or 1 that represent diagonal and propto. In those cases VA covariance is always unstructured
            matrix <Type> invSr(trmsize(1,re),trmsize(1,re));invSr.setZero();
            
            // unstructured Var.cov for cstruc<5 except block diagonal for cstruc = -1, and kronecker >5
            matrix<Type> Arm(trmsize(1,re),trmsize(1,re));
            matrix<Type> Sr(trmsize(1,re), trmsize(1,re));
            Arm.setZero();Sr.setZero();
            
            for (int d=0; d<(trmsize(1,re)); d++){ // diagonals of varcov
              Arm(d,d)=exp(lg_Ar(sdcounter));
              sdcounter++;
            }
            
            for (int d=0; d<(trmsize(1,re)); d++){
              for (int r=d+1; r<(trmsize(1,re)); r++){
                Arm(r,d)=lg_Ar(covscounter);
                covscounter++;
              }}
            
            // add terms to cQ
            matrix<Type> ArmMat = Arm*Arm.transpose();
            cQ += (0.5*(dr0.middleCols(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(1,re))*ArmMat*dr0.middleCols(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(1,re)).transpose()).diagonal()).replicate(1,p);
            
          if(cstruc(re)<5){
          if(cstruc(re) == 1){ // corAR1
            Sr = gllvm::corAR1(sigma_re(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
            sigmacounter+= 2;
          }else if(cstruc(re) == 3){ // corCS
            Sr = gllvm::corCS(sigma_re(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
            sigmacounter += 2;
          }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
            // Distance matrix calculated from the coordinates for rows
            matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
            matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
            DiSc.setZero();
            DiSc.diagonal().array() += 1/sigma_re(sigmacounter);
            sigmacounter++;
            dc_scaled = dc(dccounter)*DiSc;
            if(cstruc(re)==2){ // corExp
              Sr = gllvm::corExp(sigma_re(sigmacounter), Type(0), trmsize(1,re), dc_scaled);
              sigmacounter++;
            } else if(cstruc(re)==4) { // corMatern
              Sr = gllvm::corMatern(sigma_re(sigmacounter), Type(1), sigma_re(sigmacounter+1), trmsize(1,re), dc_scaled);
              sigmacounter += 2;
            }
            dccounter++;
          }
          
          //TMB's matinvpd function: inverse of matrix with logdet for free
          CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
          logdetSr = res[0];
          invSr = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
          }else{
            invSr = pow(sigma_re(sigmacounter), -2)*proptoMats(propcount)(0);
            logdetSr = proptoMats(propcount)(1)(0) + 2*proptoMats(propcount)(0).cols()*log_sigma(sigmacounter);
            sigmacounter++;
            propcount++;
          }
          
          if(re==0){
          //diagonal RE
            nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(0,trmsize(1,re)).transpose()*(invSr*r0r.col(0).segment(0,trmsize(1,re)))).sum());
          }else{
          //struc RE
            nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)).transpose()*(invSr*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)))).sum());
          }
          // determinants of each block of the covariance matrix
          nll -= 0.5*(trmsize(1,re)-logdetSr);
          
          }
        }else{
          Type logdetSr = 0;
          
          if(cstruc(re)<0 || cstruc(re) > 5){
            matrix <Type> invSr(trmsize(0,re),trmsize(0,re));invSr.setZero();
            
            matrix<Type> sds = Eigen::MatrixXd::Zero(trmsize(0,re),trmsize(0,re));
            sds.diagonal() =  sigma_re.segment(sigmacounter, trmsize(0,re));
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
            matrix <Type> Ir = Eigen::MatrixXd::Identity(SrL.cols(),SrL.cols());
            matrix <Type> SrIL(SrL.cols(),SrL.cols());
            SrIL = SrL.template triangularView<Eigen::Lower>().solve(Ir);
            SrIL = SrIL.transpose()*SrIL;
            invSr=SrIL*SrIL.transpose();
            logdetSr = 2*SrL.diagonal().array().log().sum();
            
            matrix<Type> Arm(trmsize(0,re),trmsize(0,re));
            if(cstruc(re)<0){
              for (int q=0; q<trmsize(1,re); q++){//loop over blocks
                Arm.setZero();  
                for (int d=0; d<(trmsize(0,re)); d++){ // diagonals of varcov
                  Arm(d,d)=exp(lg_Ar(sdcounter));
                  sdcounter++;
                }

                matrix<Type> ArmMat = Arm*Arm.transpose();
                
                cQ += (0.5*(dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re))*ArmMat*dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re)).transpose()).diagonal()).replicate(1,p);
                
                if(re==0){
                  nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(trmsize(0,re)*q,trmsize(0,re)).transpose()*(invSr*r0r.col(0).segment(trmsize(0,re)*q,trmsize(0,re)))).sum());  
                }else{
                  nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum()+trmsize(0,re)*q,trmsize(0,re)).transpose()*(invSr*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum()+trmsize(0,re)*q,trmsize(0,re)))).sum());
                }
                
                // determinants of each block of the covariance matrix
                nll -= 0.5*(trmsize(0,re)-logdetSr);
              }
            }else if(cstruc(re) > 5){
              matrix<Type> invMat(trmsize(1,re), trmsize(1,re));
              invMat.setZero();
              
                if(cstruc(re)>6){
                // here we need to calculate the inverse of our second covariance matrix
                // as we have a kronecker product, and variances are in SrL, the matrices below are correlation matrices.
                // this keeps the number of constraints similar to the proptoustruc case
                matrix<Type> Sr(trmsize(1,re), trmsize(1,re));
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
                  DiSc.diagonal().array() += 1/sigma_re(sigmacounter);
                  sigmacounter++;
                  dc_scaled = dc(dccounter)*DiSc;
                  if(cstruc(re) == 8){ // corExp
                    Sr = gllvm::corExp(Type(1), Type(0), trmsize(1,re), dc_scaled);
                  } else if(cstruc(re) == 10) { // corMatern
                    Sr = gllvm::corMatern(Type(1), Type(1), sigma_re(sigmacounter), trmsize(1,re), dc_scaled);
                    sigmacounter += 1;
                  }
                  dccounter++;
                }
                
                //TMB's matinvpd function: inverse of matrix with logdet for free
                CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
                logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*res[0];
                invMat = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
                REPORT(Sr);
                }else if(cstruc(re)==6){
                  // we have a known inverse here
                  logdetSr = logdetSr*trmsize(1,re) + trmsize(0,re)*proptoMats(propcount)(1)(0); //logdet kronecker
                  invMat = proptoMats(propcount)(0);
                  propcount ++;
                }
              
              Arm.setZero();  
              for (int d=0; d<(trmsize(0,re)); d++){ // diagonals of varcov
                Arm(d,d)=exp(lg_Ar(sdcounter));
                sdcounter++;
              }

              matrix<Type> ArmMat = Arm*Arm.transpose();
              
              vector<Type> ArmP(trmsize(1,re));
              ArmP(0) = 1; //identifiability
              for (int d=1; d<(trmsize(1,re)); d++){ // diagonals of varcov
                ArmP(d)=exp(lg_Ar(sdcounter));
                sdcounter++;
              }

              for (int q=0; q<trmsize(1,re); q++){//loop over blocks
                cQ += ((0.5*(dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re))*ArmMat*dr0.middleCols(trmsize.row(1).head(re).sum()+trmsize(0,re)*q, trmsize(0,re)).transpose()).diagonal())*ArmP(q)*ArmP(q)).replicate(1,p);
              }
              
              if(re==0){
                Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(0, trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
                nll -= ArmP.size()*Arm.diagonal().array().log().sum() + Arm.cols()*ArmP.array().log().sum() - 0.5*((invMat*(ArmP.array()*ArmP.array()).matrix().asDiagonal()).trace()*(invSr*ArmMat).trace()+(bm*invMat*bm.transpose()*invSr).trace());
              }else{
                Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bm = Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(0,re)*trmsize(1,re)).data(), trmsize(0,re), trmsize(1,re));
                nll -= ArmP.size()*Arm.diagonal().array().log().sum() + Arm.cols()*ArmP.array().log().sum() - 0.5*((invMat*(ArmP.array()*ArmP.array()).matrix().asDiagonal()).trace()*(invSr*ArmMat).trace()+(bm*invMat*bm.transpose()*invSr).trace());                                                   
              }
              
              // determinants of each block of the covariance matrix
              nll -= 0.5*(trmsize(0,re)*trmsize(1,re)-logdetSr);
            }
            
          }else if(cstruc(re) == 0 && trmsize(0,re) > 1){
            // diag() with nc > 1 covariates: nc independent per-covariate prior variances.
            // Va covariance is nc x nc diagonal per group, looped over nl groups.
            matrix<Type> sds = Eigen::MatrixXd::Zero(trmsize(0,re), trmsize(0,re));
            sds.diagonal() = sigma_re.segment(sigmacounter, trmsize(0,re));
            sigmacounter += trmsize(0,re);
            matrix<Type> invSrDiag(trmsize(0,re), trmsize(0,re));
            invSrDiag.setZero();
            invSrDiag.diagonal() = sds.diagonal().array().pow(-2);
            Type logdetSr = 2 * sds.diagonal().array().log().sum(); // log-det per group

            matrix<Type> Arm(trmsize(0,re), trmsize(0,re));
            for(int q = 0; q < trmsize(1,re); q++){
              Arm.setZero();
              for(int d = 0; d < trmsize(0,re); d++){
                Arm(d,d) = exp(lg_Ar(sdcounter));
                sdcounter++;
              }
              matrix<Type> ArmMat = Arm * Arm.transpose();
              int col_offset = trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum() + trmsize(0,re)*q;
              cQ += (0.5*(dr0.middleCols(col_offset, trmsize(0,re)) * ArmMat * dr0.middleCols(col_offset, trmsize(0,re)).transpose()).diagonal()).replicate(1,p);
              int r0r_offset = trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum() + trmsize(0,re)*q;
              nll -= Arm.diagonal().array().log().sum()
                     - 0.5*((invSrDiag * ArmMat).trace() + (r0r.col(0).segment(r0r_offset, trmsize(0,re)).transpose() * (invSrDiag * r0r.col(0).segment(r0r_offset, trmsize(0,re)))).sum());
              nll -= 0.5*(trmsize(0,re) - logdetSr);
            }
          }else{
          Eigen::DiagonalMatrix<Type, Eigen::Dynamic> Arm(trmsize(1,re));
          matrix<Type> Sr(trmsize(1,re), trmsize(1,re));Sr.setZero();

          for (int d=0; d<(trmsize(1,re)); d++){ // diagonals of varcov
            Arm.diagonal()(d)=exp(lg_Ar(sdcounter));
            sdcounter++;
          }
          // add terms to cQ
          cQ += (0.5*(dr0.middleCols(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(1,re))*Arm*Arm*dr0.middleCols(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(), trmsize(1,re)).transpose()).eval().diagonal()).replicate(1,p);

          // We build the actual covariance matrix
          // This can straightforwardly be extended to estimate correlation between effects
          matrix <Type> invSr(trmsize(1,re), trmsize(1,re));invSr.setZero();
          // diagonal row effect
          if(cstruc(re)<5){
          if(cstruc(re) == 0){
            // inverse and log determinant are straighforwardly available here
            logdetSr = 2*trmsize(1,re)*log(sigma_re(sigmacounter));
            sigmacounter++;
          }else if(cstruc(re) == 1){ // corAR1
            Sr = gllvm::corAR1(sigma_re(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
            sigmacounter+= 2;
          }else if(cstruc(re) == 3){ // corCS
            Sr = gllvm::corCS(sigma_re(sigmacounter), log_sigma(sigmacounter+1), trmsize(1,re));
            sigmacounter += 2;
          }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
            // Distance matrix calculated from the coordinates for rows
            matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
            matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
            DiSc.setZero();
            DiSc.diagonal().array() += 1/sigma_re(sigmacounter);
            sigmacounter++;
            dc_scaled = dc(dccounter)*DiSc;
            if(cstruc(re)==2){ // corExp
              Sr = gllvm::corExp(sigma_re(sigmacounter), Type(0), trmsize(1,re), dc_scaled);
              sigmacounter++;
            } else if(cstruc(re)==4) { // corMatern
              Sr = gllvm::corMatern(sigma_re(sigmacounter), Type(1), sigma_re(sigmacounter+1), trmsize(1,re), dc_scaled);
              sigmacounter += 2;
            }
            dccounter++;
          }
          if(cstruc(re)>0){
            //TMB's matinvpd function: inverse of matrix with logdet for free
            CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
            logdetSr = res[0];
            invSr = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
          }
          }else{
            invSr = pow(sigma_re(sigmacounter), -2)*proptoMats(propcount)(0);
            logdetSr = proptoMats(propcount)(1)(0) + 2*proptoMats(propcount)(0).cols()*log_sigma(sigmacounter);
            sigmacounter++;
            propcount++;
          }
          
          if(re==0){
            if(cstruc(re)==0){
              nll -= Arm.diagonal().array().log().sum() - 0.5*pow(sigma_re(sigmacounter-1), -2)*(Arm.diagonal().array().pow(2).sum()+(r0r.col(0).segment(0,trmsize(1,re)).transpose()*r0r.col(0).segment(0,trmsize(1,re))).sum());              
            }else{
              nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*Arm*Arm).trace()+(r0r.col(0).segment(0,trmsize(1,re)).transpose()*(invSr*r0r.col(0).segment(0,trmsize(1,re)))).sum());
            }
          }else{
            if(cstruc(re)==0){
              nll -= Arm.diagonal().array().log().sum() - 0.5*pow(sigma_re(sigmacounter-1), -2)*(Arm.diagonal().array().pow(2).sum()+(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)).transpose()*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re))).sum());
            }else{
              nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*Arm*Arm).trace()+(r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)).transpose()*(invSr*r0r.col(0).segment(trmsize.row(1).cwiseProduct(trmsize.row(0)).head(re).sum(),trmsize(1,re)))).sum());
            }
          }
          // determinants of each block of the covariance matrix
          nll -= 0.5*(trmsize(1,re)-logdetSr);
        }
        }
      }
    }
