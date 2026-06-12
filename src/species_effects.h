// species_effects.h
// Shared species-specific random effects NLL for gllvm.cpp and gllvm_HO.cpp.
// Requires in scope: random, xb, Br, B, sigmaB, sigmaij, cs, colMatBlocksI,
//   nncolMat, Abranks, Abstruc, Abb, cQ, n, p, eta, nll.
#pragma once
    if((random(1)>0) || (random(3)>0)){
      if(random(1)>0){
        // random slopes in TMBtrait.R
        eta += xb*Br;  
      }else if(random(3)>0){
        //random slopes in gllvm.TMB
        if(B.rows()==xb.cols()){
          //RE means or community-level effects
          eta += (xb*B).replicate(1,p);
        }
        eta += xb*Br;
      }
      matrix <Type> Spr(xb.cols(),xb.cols());
      matrix <Type> SprI(xb.cols(),xb.cols());
      matrix <Type> SprIL(xb.cols(),xb.cols());
      
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
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      
      // vector<matrix<Type>> colCorMatIblocks(colMatBlocksI.size()-1);
      // Type logdetColCorMat =0;
      vector <Type> rhoSP(1);
      if(random(1)>0 && sigmaB.size()>xb.cols()){
        rhoSP.resize(sigmaB.size()-xb.cols());
        
        if(colMatBlocksI(0)(0,0) == p){ //extra check, just in case
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
        }
        
      }else if(random(3)>0){
        rhoSP.resize(sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1));
        
        if(colMatBlocksI(0)(0,0) == p){
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
      }
      
      //Variational components
      int ncov = xb.cols();
      
      if(Abstruc == 0){
        //((Abb.size() ==(p*ncov)) || (Abb.size() == (p*ncov+p*ncov*(ncov-1)/2))) && Abranks(0) == 0
        // Ab.struct == "diagonal" or "blockdiagonal", i.e., independence over species
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        // vector<Type> SprIAtrace(p);
        vector<matrix <Type>> SArm(p);
        
        for(int j=0; j<p; j++){//limit the scope of SArmLst(j) to this iteration for memory efficiency
          SArm(j).resize(ncov,ncov);
          SArm(j).setZero();
          for (int d=0; d<ncov; d++){ // diagonals of varcov
            SArm(j)(d,d)=exp(Abb(sdcounter));
            sdcounter++;
          }
          
          if((Abb.size()>(p*ncov))){ // unstructured block Var.cov
            for (int d=0; d<(ncov); d++){
              for (int r=d+1; r<(ncov); r++){
                SArm(j)(r,d)=Abb(covscounter);
                covscounter++;
              }}
          }
          
          nll -= SArm(j).diagonal().array().log().sum();
          cQ.col(j) += 0.5*((xb*SArm(j)*SArm(j).transpose()).cwiseProduct(xb)).rowwise().sum();
        }
        // tr(S⁻¹A)+Br'S⁻¹Br+logdet(S)
        if((colMatBlocksI(0)(0,0)==p) && (rhoSP.size()==1)){
          int sp = 0;
          if(nncolMat.rows()<p){
            for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
              matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
              colCorMatI.setZero();
              gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
              nll -= -0.5*(Br.middleCols(sp, colMatBlocksI(cb+1).cols())*colCorMatI*Br.middleCols(sp, colMatBlocksI(cb+1).cols()).transpose()*SprI).trace();
              
              for(int j=0; j<(colMatBlocksI(cb+1).cols()); j++){
                nll -= -0.5*colCorMatI(j,j)*(SprI*SArm(sp+j)*SArm(sp+j).transpose()).trace();  
              }
              sp += colMatBlocksI(cb+1).cols();
            }
          }else{
            for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
              Eigen::SparseMatrix<Type> colCorMatUI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
              gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb+1).cols()));
              Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
              nll -= -0.5*(Br.middleCols(sp, colMatBlocksI(cb+1).cols())*colCorMatI*Br.middleCols(sp, colMatBlocksI(cb+1).cols()).transpose()*SprI).trace();
              vector<Type>colCorMatIdiag = colCorMatI.diagonal();
              for(int j=0; j<(colMatBlocksI(cb+1).cols()); j++){
                nll -= -0.5*colCorMatIdiag(j)*(SprI*SArm(sp+j)*SArm(sp+j).transpose()).trace();  
              }
              sp += colMatBlocksI(cb+1).cols();
            }
          }
          
          //determinant
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
          
        }else if((colMatBlocksI(0)(0,0)==p) && (rhoSP.size()>1)){
          // calculate trace term without building the covariance matrix
          // requires a little bit of work with a temporary matrix..
          int sp = 0;
          for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
            int blocksize = colMatBlocksI(cb+1).cols();
            vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUI(ncov);
            
            for (int d=0; d<ncov; d++){
              colCorMatUI(d).resize(blocksize,blocksize);
              gllvmutils::nngp(colCorMatUI(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              Br.row(d).middleCols(sp, blocksize) *= colCorMatUI(d);
              colCorMatUI(d) = colCorMatUI(d).transpose(); //now lower triangular, needed to expose .col below
            }
            nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            
            matrix<Type>tempMat(ncov,ncov);
            for(int j=0; j<blocksize; j++){
              for (int d=0; d<ncov; d++){
                for (int d2=d+1; d2<ncov; d2++){
                  tempMat(d2,d) = (colCorMatUI(d).col(j).transpose()*colCorMatUI(d2).col(j)).sum()*SprI(d,d2);
                  tempMat(d,d2) = tempMat(d2,d);
                }
                tempMat(d,d) = (colCorMatUI(d).col(j).transpose()*colCorMatUI(d).col(j)).sum()*SprI(d,d);
              }
              nll -= -0.5*(tempMat*SArm(sp+j)*SArm(sp+j).transpose()).trace();
            }
            sp += colMatBlocksI(cb+1).cols();
          }
          
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }else{
          for(int j=0; j<p; j++){
            nll -= -0.5*(SprI*SArm(j)*SArm(j).transpose()).trace();
          }
          nll -= 0.5*(p*ncov-p*logdetSpr); 
          nll -= -0.5*(Br*Br.transpose()*SprI).trace();
        }

        
      }else if(Abstruc == 1){
        //(Abb.size()==(ncov+p-1+p*(p-1)/2)) || (Abb.size() == (p+ncov-1 + p*(p-1)/2 + ncov*(ncov-1)/2)) || (Abb.size() == (p+ncov-1+ncov*(ncov-1)/2+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum())) || (Abb.size() == (p+ncov-1+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "MNdiagonal" OR "MNunstructured" with blocks due to Phylogeny
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        matrix<Type>SArmR(ncov, ncov);//row covariance matrix
        SArmR.setZero();
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p+ncov-1;
        
        //row covariance matrix
        for (int d=0; d<(ncov); d++){
          SArmR.diagonal()(d)=exp(Abb(sdcounter));
          sdcounter++;
        }
        
        if((Abb.size()>(ncov+p-1+p*(p-1)/2)) || (Abb.size()>(p+ncov-1+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))){ // unstructured row covariance
          for (int d=0; d<(ncov); d++){
            for (int r=d+1; r<(ncov); r++){
              SArmR(r,d)=Abb(covscounter);
              covscounter++;
            }}
        }
        
        // determinant of kronecker product of two matrices based on their cholesky factors: first part
        nll -= p*SArmR.diagonal().array().log().sum();
        SArmR *= SArmR.transpose();
        
        //for later use
        matrix<Type> xbSArmxb = ((xb*SArmR).cwiseProduct(xb)).rowwise().sum();
        
        int sp = 0;//tracking how many species we have had so far
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          if(Abranks(cb)<blocksize){
            //use sparse matrices if we can to speed things up with many species
            // Eigen::SparseMatrix<Type>SArmC(blocksize,blocksize);
            // std::vector<T> tripletList;
            matrix<Type> SArmC(blocksize, CppAD::Integer(Abranks(cb)));
            SArmC.setZero();
            Eigen::Vector<Type, Eigen::Dynamic>SArmCD(blocksize);//remaining variances
            SArmCD.setZero();
            // set-up column (block) covariance matrix
            int jstart = 0;
            if(cb==0){
              //first entry fixed for identifiability
              SArmC(0,0) = 0.3;
              jstart = 1;
            }
            
            for (int j=jstart; j<blocksize; j++){
              if(j<Abranks(cb)){
                SArmC(j,j) = exp(Abb(sdcounter));
                SArmCD(j) = 0;
              }else{
                SArmCD(j) = exp(Abb(sdcounter));
              }
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<blocksize; r++){
                SArmC(r,j) =  Abb(covscounter);
                covscounter++;
              }
            }
            // determinant of kronecker product of two matrices based on their cholesky factors: second part
            nll -= SArmR.rows()*(SArmC.diagonal().array().log().sum() + SArmCD.array().tail(blocksize-SArmC.cols()).log().sum());
            // matrix<Type>SArmP = SArmC*SArmC.transpose();
            // SArmP.diagonal() += SArmCD.array().pow(2).matrix();
            
            //tr(S⁻¹A) and Br'S⁻¹Br
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                //not sparse
                matrix<Type>colCorMatI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(((SArmC.transpose()*colCorMatI*SArmC).trace()+(SArmCD.array().pow(2)*colCorMatI.diagonal().array()).sum())*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else{
                //NN sparse approximation
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Type traceUUtSArmP = ((SArmC.transpose()*colCorMatUI)*(colCorMatUI.transpose()*SArmC)).trace();
                //remaining diagonal entries
                // equivalent to diag(SArmCDs^2)UU', but via the cholesky instead
                // diag(SArmCDs)U, or the squared norm of the rows (),
                // so that I do not need to compute the whole product
                for (int j=0; j<blocksize; j++){
                  vector<Type> temp = colCorMatUI.col(j).cwiseProduct(SArmCD);
                  traceUUtSArmP += temp.pow(2).sum();
                }
                // Type traceold = (colCorMatUI*colCorMatUI.transpose()*SArmP).trace();
                // REPORT(traceold);
                // REPORT(traceUUtSArmP);
                nll -= -0.5*(traceUUtSArmP*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
                
                // nll -= -0.5*((colCorMatUI*colCorMatUI.transpose()*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              //tr(S^-1A)
              vector<Eigen::SparseMatrix<Type>> colCorMatUIBlocks(ncov);
              for (int d=0; d<ncov;d++){
                colCorMatUIBlocks(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIBlocks(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                //Br'S⁻¹Br
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIBlocks(d);
              }
              //this is written in a way that keeps dimensions of the
              //objects small thanks to the reduced rank
              //which is much more friendly in memory
              //based on the formulation A_m = LL'+D, where D is diagonal given by SArmCD^2
              // trace: tr((LL'+D)UU')*Ar(d,d2)*SprI(d,d2), split over the addition for efficiency
              
              matrix<Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
              vector<Type> temp(blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmC.transpose()*colCorMatUIBlocks(d);
                for (int d2=d+1; d2<ncov;d2++){
                  // nll -= -(SArmP*colCorMatUIBlocks(d)*colCorMatUIBlocks(d2).transpose()).trace()*SArmR(d,d2)*SprI(d,d2);
                  //reduced rank part: L
                  nll -= -((tempMat*(colCorMatUIBlocks(d2).transpose()*SArmC)).trace())*SArmR(d,d2)*SprI(d,d2);
                  // apparently faster than "just" computing the diagonal entries via an elementwise product
                  // I suppose Eigen can do some optimization with this
                  vector<Type> diags(blocksize);
                  diags.setZero();
                  for (int j=0; j<blocksize; j++){
                    diags += colCorMatUIBlocks(d).col(j).cwiseProduct(colCorMatUIBlocks(d2).col(j));
                  }
                  nll -= -(SArmCD.array().pow(2)*diags).sum()*SArmR(d,d2)*SprI(d,d2);
                }
                //reduced rank part
                nll -= -0.5*(tempMat).rowwise().squaredNorm().sum()*SArmR(d,d)*SprI(d,d);
                Type trace = 0;
                //remaining diagonal entries
                //need to do this column-wise so I do not need to compute the whole product
                //similar to above, trace via squared norm
                for (int j=0; j<blocksize; j++){
                  temp = colCorMatUIBlocks(d).col(j).cwiseProduct(SArmCD);
                  trace += temp.pow(2).sum();
                }
                nll -= -0.5*trace*SArmR(d,d)*SprI(d,d);
              }
              //nll -= -0.5*(SprIL.transpose()*Br.middleCols(sp, blocksize)).rowwise().squaredNorm().sum();//is slower it seems
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
            //retrieve diagonal outside loop for more efficient memory handling
            //calling it inside the loop constructs a temporary for each i
            matrix<Type>SArmPdiag = SArmC.rowwise().squaredNorm().array()+SArmCD.array().pow(2);
            cQ.middleCols(sp, blocksize).noalias() += 0.5*xbSArmxb*SArmPdiag.transpose();
            // for (int i=0; i<n;i++){
            //   cQ.row(i).middleCols(sp, blocksize) += 0.5*SArmPdiag*xbSArmxb(i);
            // } 
            sp += blocksize;
          }else{
            matrix<Type>SArmC;
            SArmC.resize(blocksize,blocksize);//column VA covariance
            SArmC.setZero();
            
            int jstart = 0;
            // set-up column (block) covariance matrix
            if(cb==0){
              SArmC(0,0) = 0.3;//first entry fixed for identifiability
              jstart = 1;
            }
            
            for (int j=jstart; j<SArmC.cols(); j++){
              SArmC(j,j) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                SArmC(r,j) = Abb(covscounter);
                covscounter++;
              }
            }
            
            // determinant of kronecker product of two matrices based on their cholesky factors: second part
            nll -= SArmR.rows()*SArmC.diagonal().array().log().sum();
            matrix<Type>SArmP = SArmC*SArmC.transpose();
            
            //tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                //not sparse
                matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*((colCorMatI*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else{
                //NN sparse approximation
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                nll -= -0.5*((colCorMatUI*colCorMatUI.transpose()*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              vector<Eigen::SparseMatrix<Type>> colCorMatUIBlocks(ncov);
              for (int d=0; d<ncov;d++){
                colCorMatUIBlocks(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIBlocks(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              }
              //tr(S⁻¹A)
              matrix<Type>tempMat(blocksize,blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmP*colCorMatUIBlocks(d);
                for (int d2=d+1; d2<ncov;d2++){
                  nll -= -(tempMat*colCorMatUIBlocks(d2).transpose()).diagonal().sum()*SArmR(d,d2)*SprI(d,d2);
                }
                nll -= -0.5*(tempMat*colCorMatUIBlocks(d).transpose()).diagonal().sum()*SArmR(d,d)*SprI(d,d);
              
              //Br'S⁻¹Br
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIBlocks(d);
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
            //retrieve diagonal outside loop for more efficient memory handling
            //calling it inside the loop constructs a temporary for each i
            matrix<Type> SArmPdiag = SArmP.diagonal();
            // for (int i=0; i<n;i++){
            //   cQ.row(i).middleCols(sp, blocksize) += 0.5*SArmPdiag*xbSArmxb(i);
            // } 
            cQ.middleCols(sp, blocksize).noalias() += 0.5*xbSArmxb*SArmPdiag.transpose();
            sp += blocksize;
          }
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 2){
        //(Abb.size() == ((ncov*p +ncov*p*(p-1)/2))) || (Abb.size()==(ncov*p+sum(nsp)*(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "diagonalCL2" with block structure due to phylogeny
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        
        int sp = 0;//keeping track of #species we've had due to blocks
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          matrix<Type>colCorMatI(blocksize,blocksize);
          vector<Eigen::SparseMatrix<Type>> colCorMatUIs(ncov);
          
          if((rhoSP.size()) == 1){
            gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
            if(nncolMat.rows()<p){
              nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }else if(nncolMat.rows()==p){
              colCorMatUIs(0).resize(blocksize, blocksize);
              gllvmutils::nngp(colCorMatUIs(0), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
              //sparse approximation to inverse
              nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatUIs(0)*colCorMatUIs(0).transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
          }else{
            //NN sparse approximation
            //Br'S⁻¹Br'
            for (int d=0; d<(ncov); d++){
              colCorMatUIs(d).resize(blocksize, blocksize);
              gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
            }
            nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
          }
          
          if(Abranks(cb)<colMatBlocksI(0)(cb+1,0)){
            //can use truncated matrices
            for (int d=0; d<(ncov); d++){
              matrix<Type> SArm(blocksize,CppAD::Integer(Abranks(cb)));
              vector<Type> SArmD(blocksize);
              for (int j=0; j<blocksize; j++){ // diagonals
                if(j<Abranks(cb)){
                  SArm(j,j) = exp(Abb(sdcounter));
                  SArmD(j) = 0;
                }else{
                  SArmD(j) = exp(Abb(sdcounter));
                }
                
                sdcounter++;
              }
              
              for (int j=0; j<Abranks(cb); j++){
                for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                  SArm(r,j) = Abb(covscounter);
                  covscounter++;
                }
              }
              //determinant with first term reduced rank and second term remaining diagonal entries
              nll -= SArm.diagonal().array().log().sum() + SArmD.tail(blocksize-SArm.cols()).log().sum();
              
              matrix<Type> SArmpd = SArm*SArm.transpose();
              SArmpd.diagonal() += SArmD.pow(2).matrix();
              //remaining likelihood terms
              //tr(S⁻¹A)
              if((rhoSP.size() == 1)){
                if(nncolMat.rows()<p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatI*SArmpd).trace());
                }else{
                  //sparse approximation to inverse
                  nll -= -0.5*(SprI(d,d)*(colCorMatUIs(0)*colCorMatUIs(0).transpose()*SArmpd).trace());
                }
              }else{
                //sparse approximation to inverse
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArmpd).trace()*SprI(d,d);
              }
              
              matrix<Type>xbxb = (xb.col(d).cwiseProduct(xb.col(d))).rowwise().sum();
              cQ.middleCols(sp, blocksize) += 0.5*xbxb*SArmpd.diagonal().transpose();
            }
          }else{
            for (int d=0; d<(ncov); d++){
              matrix<Type> SArm(blocksize,blocksize);//limit the scope of SArm(d) to this iteration for memory efficiency due to potentially large d
              SArm.setZero();
              for (int j=0; j<blocksize; j++){ // diagonals
                SArm(j,j)=exp(Abb(sdcounter));
                sdcounter++;
              }
              
              for (int j=0; j<blocksize; j++){
                for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                  SArm(r,j)=Abb(covscounter);
                  covscounter++;
                }
              }
              
              nll -= SArm.diagonal().array().log().sum();
              
              SArm *= SArm.transpose();
              //tr(S⁻¹A)
              if((rhoSP.size() == 1)){
                if(nncolMat.rows()<p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatI*SArm).trace());
                }else if(nncolMat.rows()==p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatUIs(0)*colCorMatUIs(0).transpose()*SArm).trace());
                }
              }else{
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArm).trace()*SprI(d,d);
              }
              
              matrix<Type>xbxb = (xb.col(d).cwiseProduct(xb.col(d))).rowwise().sum();
              cQ.middleCols(sp,blocksize) += 0.5*xbxb*SArm.diagonal().transpose();
            }
          }
          
          sp += blocksize;
          
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 3){
        //(Abb.size() == (p*ncov + p*(p-1)/2 + p-colCorMatIblocks.size())) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + p*ncov*(ncov-1)/2 + p*(p-1)/2)) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + p*ncov*(ncov-1)/2 + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum())) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "CL1" OR "diagonalCL1". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        // always has p*p matrix with fixed diagonals, in combination with diagonal/blockdiagonal matrix for sp*ncov
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov+p-(colMatBlocksI.size()-1);
        int sp = 0;
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          //construct blockdiagonal list
          vector<matrix<Type>> SArmb(blocksize);
          
          for (int j=0; j<blocksize; j++){
            SArmb(j).resize(ncov,ncov);
            SArmb(j).setZero();
            for (int d=0; d<(ncov); d++){ // diagonals of varcov
              SArmb(j)(d,d) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            if((Abb.size()>(p*ncov + p-(colMatBlocksI.size()-1) + p*(p-1)/2)) || (Abb.size()>(p*ncov + p-(colMatBlocksI.size()-1) + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))){ // unstructured block Var.cov
              //only for "blockdiagonalsp"
              for (int d=0; d<(ncov); d++){
                for (int r=d+1; r<(ncov); r++){
                  SArmb(j)(r,d) = Abb(covscounter);
                  covscounter++;
                }}
            }
          }
          
          //construct p by p covariance matrix
          if(Abranks(cb)<blocksize){
            //construct as sparse matrix
            //use sparse matrices if we can to speed things up with many species
            matrix<Type>SArmC(blocksize,CppAD::Integer(Abranks(cb)));
            SArmC.setZero();
            vector<Type> SArmCD(blocksize);
            SArmCD.setZero();
            
            // set-up column (block) covariance matrix
            SArmC(0,0) = 1;//first entry fixed for identifiability
            for (int j=1; j<blocksize; j++){
              if(j<Abranks(cb)){
                SArmC(j,j) = exp(Abb(sdcounter));
              }else{
                SArmCD(j) = exp(Abb(sdcounter));
              }
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<blocksize; r++){
                SArmC(r,j) =  Abb(covscounter);
                covscounter++;
              }
            }
            
            //p*p matrix in the middle: correlation for identifiability
            matrix<Type> SArmP = SArmC*SArmC.transpose();
            SArmP.diagonal() += SArmCD.pow(2).matrix();
            SArmCD.array() /= SArmP.diagonal().array().cwiseSqrt();
            SArmC.array().colwise() /= SArmP.diagonal().array().cwiseSqrt();
            
            //determinant of this matrix is nsp*sum(log(diag(SArmC)))+2*sum(log(diags(SArmb)))
            nll -=  ncov*SArmC.diagonal().array().log().sum() + ncov*SArmCD.tail(blocksize-SArmC.cols()).log().sum();
            for (int j=0; j<blocksize; j++){
              nll -= SArmb(j).diagonal().array().log().sum();
            }
            
            SArmP = gllvmutils::cov2cor(SArmP);
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              matrix<Type> colCorMatI(blocksize, blocksize);
              if(nncolMat.rows()<p){
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
                matrix<Type>tempMat(ncov,ncov);
                for (int j=0; j<blocksize;j++){
                  tempMat = SprI*SArmb(j);
                  for (int j2=j+1; j2<blocksize;j2++){
                    // temp.diagonal().fill(SArmP(j,j2));
                    nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                  }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
                }
                
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                colCorMatI = colCorMatUI*colCorMatUI.transpose();
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }
              //write trace separately so we don't need to compute the whole product
              //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
              matrix<Type>tempMat(ncov,ncov);
              for (int j=0; j<blocksize;j++){
                tempMat = SprI*SArmb(j);
                for (int j2=j+1; j2<blocksize;j2++){
                  // temp.diagonal().fill(SArmP(j,j2));
                  nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
              }
            }else{
              vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
                colCorMatUIs(d) = colCorMatUIs(d).transpose();//to expose .col below
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              matrix<Type>tempMat(ncov,ncov);
              //this trace cannot be separated further due to element-wise product..
              for (int j=0; j<blocksize;j++){
                for (int j2=j+1; j2<blocksize;j2++){
                  for (int d=0; d<(ncov); d++){
                    for (int d2=0; d2<(ncov); d2++){
                      tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j2)).sum()*SprI(d,d2);
                    }
                  }
                  nll-= -(tempMat.cwiseProduct(SArmb(j)*SArmb(j2).transpose())*SArmP(j,j2)).trace();
                }
                for (int d=0; d<(ncov); d++){
                  for (int d2=d+1; d2<(ncov); d2++){
                    tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j)).sum()*SprI(d,d2);
                    tempMat(d2,d) = tempMat(d,d2);
                  }
                  tempMat(d,d) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d).col(j)).sum()*SprI(d,d);
                }
                nll-= -0.5*(tempMat.cwiseProduct(SArmb(j)*SArmb(j).transpose())).trace();
              }
            }
            
            //remaining likelihood terms
            matrix<Type>SArmB(ncov, ncov);
            for (int j=0; j<blocksize;j++){
              SArmB = SArmb(j)*SArmb(j).transpose();
              for(int i=0; i<n;i++){
                cQ(i, j+sp) += 0.5*((xb.row(i)*SArmB).cwiseProduct(xb.row(i))).sum();
              }
            }
          }else{
            //use sparse matrices if we can to speed things up with many species
            matrix<Type>SArmC(blocksize,blocksize);
            SArmC.setZero();
            // set-up column (block) covariance matrix
            
            SArmC(0,0) = 1;//first entry fixed to 1 for identifiability
            for (int j=1; j<SArmC.cols(); j++){
              SArmC(j,j) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            for (int j=0; j<SArmC.cols(); j++){
              for (int r=j+1; r<SArmC.cols(); r++){
                SArmC(r,j) = Abb(covscounter);
                covscounter++;
              }
            }
            matrix<Type> SArmP = SArmC*SArmC.transpose();
            SArmC *= SArmP.diagonal().cwiseInverse().cwiseSqrt().asDiagonal();//ID constraint
            
            //determinant of this matrix is 2*nsp*sum(log(diag(SArmC)))+2*sum(log(diags(SArmb)))
            nll -=  ncov*SArmC.diagonal().array().log().sum();
            for (int j=0; j<blocksize; j++){
              nll -= SArmb(j).diagonal().array().log().sum();
            }
            SArmP = gllvmutils::cov2cor(SArmP);
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              matrix<Type> colCorMatI(blocksize, blocksize);
              if(nncolMat.rows()<p){
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }else if(nncolMat.cols()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approxiamtion to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }
              
              //write trace separately so we don't need to compute the whole product
              //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
              matrix<Type> tempMat(blocksize,blocksize);
              for (int j=0; j<blocksize;j++){
                tempMat = SprI*SArmb(j);
                for (int j2=j+1; j2<blocksize;j2++){
                    // temp.diagonal().fill(SArmP(j,j2));
                    nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
              }
            }else{
              vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
                colCorMatUIs(d) = colCorMatUIs(d).transpose();//needed to expose .col below
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              matrix<Type>tempMat(ncov,ncov);
              //this trace cannot be separated further due to element-wise product..
              for (int j=0; j<blocksize;j++){
                for (int j2=j+1; j2<blocksize;j2++){
                  for (int d=0; d<(ncov); d++){
                    for (int d2=0; d2<(ncov); d2++){
                      tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j2)).sum()*SprI(d,d2);
                    }
                  }
                  nll-= -(tempMat.cwiseProduct(SArmb(j)*SArmb(j2).transpose())*SArmP(j,j2)).trace();
                }
                for (int d=0; d<(ncov); d++){
                  for (int d2=d+1; d2<(ncov); d2++){
                    tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j)).sum()*SprI(d,d2);
                    tempMat(d2,d) = tempMat(d,d2);
                  }
                  tempMat(d,d) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d).col(j)).sum()*SprI(d,d);
                }
                nll-= -0.5*(tempMat.cwiseProduct(SArmb(j)*SArmb(j).transpose())).trace();
              }
            }

            //remaining likelihood terms
            matrix<Type>SArmB(ncov, ncov);
            for (int j=0; j<blocksize;j++){
              SArmB = SArmb(j)*SArmb(j).transpose();
              for(int i=0; i<n;i++){
                cQ(i, j+sp) += 0.5*((xb.row(i)*SArmB).cwiseProduct(xb.row(i))).sum();
              }
            }
            
          }
          sp += blocksize;
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 4){
        // Ab.struct == "CL2". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        // always has k*k matrix with fixed diagonals, in combination with k p*p matrices
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = ncov-1+p*ncov;//(p-(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()==1).count())*ncov;
        int sp = 0;
        
        matrix<Type> SArmR(ncov, ncov);
        SArmR.setZero();
        // row variance, fix first entry for identifiability (cov2cor)
        SArmR(0,0) = 1;
        for (int d=1; d<(ncov); d++){
          SArmR.diagonal()(d)=exp(Abb(sdcounter));
          sdcounter++;
        }
        
        for (int d=0; d<(ncov); d++){
          for (int r=d+1; r<(ncov); r++){
            SArmR(r,d)=Abb(covscounter);
            covscounter++;
          }}
        
        // to cholesky of correlation matrix
        SArmR.array().colwise() /= SArmR.array().rowwise().norm();
        //determinant first part
        nll -=  p*SArmR.diagonal().array().log().sum();
        //define correlation matrix for further down
        matrix<Type>SArmRc = SArmR*SArmR.transpose();
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          //construct blockdiagonal list
          vector<matrix<Type>> SArmPs(ncov);
          
          //construct p by p covariance matrix
          if(Abranks(cb)<blocksize){
            vector<vector<Type>> SArmCDs(ncov);
            
            for (int d=0; d<(ncov); d++){
              SArmPs(d).resize(blocksize,CppAD::Integer(Abranks(cb)));
              SArmPs(d).setZero();
              SArmCDs(d).resize(blocksize);
              SArmCDs(d).setZero();
              for (int j=0; j<blocksize; j++){
                if(j<Abranks(cb)){
                  SArmPs(d)(j,j) = exp(Abb(sdcounter));
                }else{
                  SArmCDs(d)(j) = exp(Abb(sdcounter));
                }
                sdcounter++;
              }
              for (int j=0; j<Abranks(cb); j++){
                for (int r=j+1; r<blocksize; r++){
                  SArmPs(d)(r,j) =  Abb(covscounter);
                  covscounter++;
                }
              }
            }
            
            //determinant of this matrix is 2*p*sum(log(diag(SArmR)))+2*sum(log(diags(SArmCs)))
            //determinant second part
            //expanding SArmPs so the remaining diagonal entries can be added to the reduced rank part
            for (int d=0; d<ncov; d++){
              nll -= SArmPs(d).diagonal().array().log().sum() + SArmCDs(d).tail(blocksize-SArmPs(d).cols()).log().sum();
            }
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                vector<Type> colCorMatIDiag = colCorMatI.diagonal();
                vector<Type> temp(blocksize);
                matrix <Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
                for (int d=0; d<ncov;d++){
                  tempMat = SArmPs(d).transpose()*colCorMatI;
                  temp = colCorMatIDiag*SArmCDs(d);
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(tempMat*SArmPs(d2)).trace()*SprI(d,d2)*SArmRc(d,d2);
                    nll -= -(temp*SArmCDs(d2)).sum()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(SArmPs(d).transpose()*colCorMatI*SArmPs(d)).trace()*SprI(d,d);
                  nll -= -0.5*(colCorMatIDiag*SArmCDs(d)*SArmCDs(d)).sum()*SprI(d,d)*SArmRc(d,d);
                }
              }else if(nncolMat.rows()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approximation to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                //this is the same as in matrix normal, except the m by m matrix is d-specific
                matrix<Type>tempMat(blocksize,blocksize);
                vector<Type> colCorMatIDiag = colCorMatI.diagonal();
                vector<Type>temp(blocksize);
                vector<Type> temp2(blocksize);
                for (int d=0; d<ncov;d++){
                  tempMat = SArmPs(d).transpose()*colCorMatUI;
                  temp2 = colCorMatIDiag*SArmCDs(d);
                  for (int d2=d+1; d2<ncov;d2++){
                    //reduced rank part: L
                    nll -= -(((tempMat)*(colCorMatUI.transpose()*SArmPs(d))).trace())*SArmRc(d,d2)*SprI(d,d2);
                    //remaining diagonal entries
                    nll -= -(temp2*SArmCDs(d2)).sum()*SArmRc(d,d2)*SprI(d,d2);
                  }
                  //reduced rank part
                  nll -= -0.5*(SArmPs(d).transpose()*colCorMatUI).rowwise().squaredNorm().sum()*SprI(d,d);
                  Type trace = 0;
                  //remaining diagonal entries
                  //need to do this column-wise so I do not need to compute the whole product
                  //similar to above, trace via squared norm
                  for (int j=0; j<blocksize; j++){
                    temp = colCorMatUI.col(j);
                    trace += (temp*SArmCDs(d)).pow(2).sum();
                  }
                  nll -= -0.5*trace*SprI(d,d);//no SArmRc because it has 1s on the diagonal
                  // nll -= -0.5*(SArmP*colCorMatUIBlocks(d).transpose()).trace()*SArmR(d,d)*SprI(d,d);
                }
              }
            }else{
              vector<Eigen::SparseMatrix<Type>> colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIs(d),colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
              }
              
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              //this is the same as in matrix normal, except the m by m matrix is d-specific
              matrix<Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
              vector<Type> temp(blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmPs(d).transpose()*colCorMatUIs(d);
                for (int d2=d+1; d2<ncov;d2++){
                  //reduced rank part: L
                  nll -= -(tempMat*(colCorMatUIs(d2).transpose()*SArmPs(d2))).trace()*SArmRc(d,d2)*SprI(d,d2);
                  // apparently faster than "just" computing the diagonal entries via an elementwise product
                  // I suppose Eigen can do some optimization with this
                  vector<Type> diags(blocksize);
                  diags.setZero();
                  for (int j=0; j<blocksize; j++){
                    diags += colCorMatUIs(d).col(j).cwiseProduct(colCorMatUIs(d2).col(j));
                  }
                  nll -= -(SArmCDs(d)*SArmCDs(d2)*diags).sum()*SArmRc(d,d2)*SprI(d,d2);
                }
                //reduced rank part
                nll -= -0.5*(tempMat).rowwise().squaredNorm().sum()*SprI(d,d);
                Type trace = 0;
                //remaining diagonal entries
                //need to do this column-wise so I do not need to compute the whole product
                //similar to above, trace via squared norm
                for (int j=0; j<blocksize; j++){
                  temp = colCorMatUIs(d).col(j);
                  trace += (temp*SArmCDs(d)).pow(2).sum();
                }
                nll -= -0.5*trace*SprI(d,d);//no SArmRc because it has 1s on the diagonal..
                // nll -= -0.5*(SArmP*colCorMatUIBlocks(d).transpose()).trace()*SArmR(d,d)*SprI(d,d);
              }
            }
          
          //species-specific submatrices in A
            matrix<Type>Aspp(ncov,ncov);
            for (int j=0; j<blocksize;j++){
              //need to build the dxd species-specific matrix in A
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  Aspp(d2,d) = ((SArmPs(d).row(j).cwiseProduct(SArmPs(d2).row(j))).sum()+(SArmCDs(d)(j)*SArmCDs(d2)(j)))*SArmRc(d2,d);
                  Aspp(d,d2) = Aspp(d2,d);
                }
                Aspp(d,d) = ((SArmPs(d).row(j).cwiseProduct(SArmPs(d).row(j))).sum()+(SArmCDs(d)(j)*SArmCDs(d)(j)));
              }
              cQ.col(sp+j).noalias() += 0.5*(xb*Aspp).cwiseProduct(xb).rowwise().sum();
            }
          //instead compute cholesky factors of the species-specific submatrices in A
            // matrix<Type>Aspp(ncov, ncov*blocksize);
            // Aspp.setZero();
            // 
            // Type sum_diag = 0;
            // Type sum_off_diag = 0;
            // 
            // for (int j=0; j<blocksize;j++){
            // for (int d = 0; d < ncov; d++) {
            //   Type SArmCDsdj = SArmCDs(d)(j);
            //   sum_diag = SArmPs(d).row(j).array().pow(2).sum() + pow(SArmCDsdj, 2);
            //   if (d > 0) {
            //     sum_diag -= (Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d).cwiseProduct(Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d))).sum();
            //   }
            //   Aspp.block(0,j*ncov,ncov,ncov)(d, d) = sqrt(sum_diag);
            //   
            //   vector<Type> SArmPsdj = SArmPs(d).row(j);
            //   matrix<Type> Asppdd = Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d);
            //   for (int d2 =d+1; d2 < ncov; d2++) {
            //       sum_off_diag = (SArmPs(d2).row(j).array() * SArmPsdj).sum() + SArmCDs(d2)(j) * SArmCDsdj;
            //       sum_off_diag *= SArmRc(d2, d);
            //       if (d > 0) {
            //         sum_off_diag -= (Aspp.block(0,j*ncov,ncov,ncov).row(d2).head(d).cwiseProduct(Asppdd)).sum();
            //       }
            //       Aspp.block(0,j*ncov,ncov,ncov)(d2, d) = sum_off_diag / Aspp.block(0,j*ncov,ncov,ncov)(d, d);
            //   }
            // }
            // }
          }else{
            for (int d=0; d<(ncov); d++){
              SArmPs(d).resize(blocksize,blocksize);
              SArmPs(d).setZero();
              
              for (int j=0; j<blocksize; j++){
                SArmPs(d)(j,j) = exp(Abb(sdcounter));
                sdcounter++;
              }
              for (int j=0; j<blocksize; j++){
                for (int r=j+1; r<blocksize; r++){
                  SArmPs(d)(r,j) =  Abb(covscounter);
                  covscounter++;
                }
              }
            }
            
            //determinant of this matrix is 2*p*sum(log(diag(SArmR)))+2*sum(log(diags(SArmCs)))
            //second part
            for (int d=0; d<ncov; d++){
              nll -= SArmPs(d).diagonal().array().log().sum();
            }
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                for (int d=0; d<ncov;d++){
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(colCorMatI*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(colCorMatI*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
                }
              }else if(nncolMat.rows()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb+1).cols()));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approximation to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                for (int d=0; d<ncov;d++){
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(colCorMatI*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(colCorMatI*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
                }
              }
            }else{
              vector<Eigen::SparseMatrix<Type>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  nll -= -(colCorMatUIs(d)*colCorMatUIs(d2).transpose()*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                }
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
              }
            }
            
            //remaining likelihood terms
            matrix<Type>tempMat(ncov,ncov);
            for (int j=0; j<blocksize;j++){
              //need to build the dxd species-specific matrix in A
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  tempMat(d2,d) = SArmPs(d).col(j).cwiseProduct(SArmPs(d2).row(j).transpose()).sum()*SArmRc(d2,d);
                  tempMat(d,d2) = tempMat(d2,d);
                }
                tempMat(d,d) = SArmPs(d).col(j).cwiseProduct(SArmPs(d).row(j).transpose()).sum();
              }
              for (int i=0; i<n;i++){
                cQ(i,j+sp) += 0.5*((xb.row(i)*tempMat).cwiseProduct(xb.row(i))).sum();
              }
            }
          }
          sp += blocksize;
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
      }else if(Abstruc == 5){
        // (Abb.size() == ((ncov*p)*(ncov*p)-(ncov*p)*((ncov*p)-1)/2)) || (Abb.size() == (ncov*p+(sum(nsp)*colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // matrix<Type> SArm(p*ncov,p*ncov);
        // Ab.struct == "unstructured". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        int sp = 0;
        
        typedef Eigen::Triplet<Type> T;
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          if(Abranks(cb)<blocksize){
            std::vector<T> tripletList;
            //can use sparse matrices
            Eigen::SparseMatrix<Type>SArm(blocksize*ncov,blocksize*ncov);
            //block structure :)
            
            for (int d=0; d<(blocksize*ncov); d++){ // diagonals of varcov
              tripletList.push_back(T(d,d,exp(Abb(sdcounter))));
              sdcounter++;
            }
            
            // unstructured block Var.cov
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(ncov*blocksize); r++){
                tripletList.push_back(T(r,j,Abb(covscounter)));
                covscounter++;
              }
            }
            
            SArm.setFromTriplets(tripletList.begin(),tripletList.end());
            
            nll -= SArm.diagonal().array().log().sum();
            
            matrix<Type> SArmb = SArm*SArm.transpose();
            
            //remaining likelihood terms
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type>SprMNI(blocksize*ncov, blocksize*ncov); // inverse of covariane
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                SprMNI = tmbutils::kronecker(colCorMatI,SprI);
                nll -= -0.5*((SprMNI*SArmb).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                Eigen::SparseMatrix<Type> SprMNI = tmbutils::kronecker(colCorMatI,tmbutils::asSparseMatrix(SprI));
                nll -= -0.5*((SprMNI*SArmb).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              matrix<Type>colCorMatUIs(ncov*blocksize,ncov*blocksize);
              colCorMatUIs.setZero();
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                Eigen::SparseMatrix<Type>colCorMatUI(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUI;
                colCorMatUIs.block(d*blocksize,d*blocksize,blocksize,blocksize) = colCorMatUI;
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //get colCorMat in same order as Ab
              vector<int> permVec(blocksize*ncov);
              Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>perm;
              int k = 0;
              for (int j=0; j<blocksize; j++){
                for (int d=0; d<ncov; d++){
                  permVec(k) = j+d*blocksize;
                  k++;
                }
              }
              perm.indices() = permVec;
              matrix<Type> Ip(blocksize,blocksize);
              Ip.setIdentity();
              Eigen::SparseMatrix<Type> kronSprI= tmbutils::kronecker(SprI,Ip).sparseView();
              Eigen::SparseMatrix<Type> SprMNI = tmbutils::asSparseMatrix(colCorMatUIs)*kronSprI*tmbutils::asSparseMatrix(colCorMatUIs).transpose();
              SprMNI = SprMNI.twistedBy(perm.transpose());
              nll -= -0.5*(SprMNI*SArmb).trace();
            }
            
            for (int j=0; j<blocksize;j++){
              //consume smore memory for some reason, but is faster in n
              cQ.col(j+sp) += 0.5*((xb*SArmb.block(j*ncov,j*ncov,ncov,ncov)).cwiseProduct(xb)).rowwise().sum();
            }
          }else{
            matrix<Type>SArm;
            //block structure :)
            SArm = matrix<Type>(blocksize*ncov,blocksize*ncov);
            SArm.setZero();
            for (int d=0; d<(blocksize*ncov); d++){ // diagonals of varcov
              SArm(d,d)=exp(Abb(sdcounter));
              sdcounter++;
            }
            
            // unstructured block Var.cov
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(ncov*colMatBlocksI(0)(cb+1)); r++){
                SArm(r,j)=Abb(covscounter);
                covscounter++;
              }
            }
            
            nll -= SArm.diagonal().array().log().sum();
            
            SArm = SArm*SArm.transpose();
            
            //remaining likelihood terms
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                matrix <Type> SprMNI = tmbutils::kronecker(colCorMatI,SprI);
                nll -= -0.5*((SprMNI*SArm).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                Eigen::SparseMatrix<Type>SprMNI = tmbutils::kronecker(colCorMatI,tmbutils::asSparseMatrix(SprI));
                nll -= -0.5*((SprMNI*SArm).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              matrix<Type>colCorMatUIs(ncov*blocksize,ncov*blocksize);
              colCorMatUIs.setZero();
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                Eigen::SparseMatrix<Type>colCorMatUI(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUI;
                colCorMatUIs.block(d*blocksize,d*blocksize,blocksize,blocksize) = colCorMatUI;
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //get colCorMat in same order as Ab
              vector<int> permVec(blocksize*ncov);
              Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>perm;
              int k = 0;
              for (int j=0; j<blocksize; j++){
                for (int d=0; d<ncov; d++){
                  permVec(k) = j+d*blocksize;
                  k++;
                }
              }
              perm.indices() = permVec;
              matrix<Type> Ip(blocksize,blocksize);
              Ip.setIdentity();
              Eigen::SparseMatrix<Type> kronSprI= tmbutils::kronecker(SprI,Ip).sparseView();
              Eigen::SparseMatrix<Type> SprMNI = tmbutils::asSparseMatrix(colCorMatUIs)*kronSprI*tmbutils::asSparseMatrix(colCorMatUIs).transpose();
              SprMNI = SprMNI.twistedBy(perm.transpose());
              nll -= -0.5*(SprMNI*SArm).trace();
            }
            
            //remaining likelihood terms
            for (int j=0; j<blocksize;j++){
              //consume smore memory for some reason, but is faster in n
              cQ.col(j+sp) += 0.5*((xb*SArm.block(j*ncov,j*ncov,ncov,ncov)).cwiseProduct(xb)).rowwise().sum();
            }
          }
          sp += blocksize;
          
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
      }
      
    }
