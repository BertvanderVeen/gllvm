#include <TMB.hpp>
#include <cmath>
#include "distrib.h"
#include "init.h"
#include "utils.h"
#include <R_ext/Error.h>

// Selkeyden vuoksi: nimeä perheiden koodit enumilla
enum Family : int {
  POISSON = 0,
    NEG_BINOMIAL = 1,
    BINOMIAL = 2,
    GAUSSIAN = 3,
    GAMMA = 4,
    TWEEDIE = 5,
    ZIP = 6,
    ORDINAL = 7,
    EXPONENTIAL = 8,
    BETA = 9,
    BETA_HURDLE = 10,
    ZINB = 11,
    ORDERED_BETA = 12,
    ZIB = 13,
    ZNIB = 14,
    BETA_BINOMIAL = 15,


    // lisää tarvittaessa muita perheitä: BINOMIAL=2, GAUSSIAN=3, ...
};

//--------------------------------------------------------
//GLLVM
//Authors: Jenni Niku, Bert van der Veen, Pekka Korhonen
//------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  //declares all data and parameters used
  DATA_MATRIX(y); // matrix of responses
  DATA_MATRIX(x); // matrix of covariates
  DATA_MATRIX(x_lv); // matrix of covariates for Reduced Rank and/or constrained ord
  DATA_IMATRIX(csb_lv);
  DATA_MATRIX(xr); // design matrix for fixed row effects
  DATA_MATRIX(xb); // design matrix for random species effects, (n, spnr)
  DATA_SPARSE_MATRIX(dr0); // design matrix for rows, ( n, nr)
  DATA_SPARSE_MATRIX(dLV); // design matrix for latent variables, (n, nu)
  DATA_IMATRIX(cs); // matrix with indices for correlations of random species effects
  // DATA_SPARSE_MATRIX(colMat); //lower cholesky of column similarity matrix
  DATA_STRUCT(colMatBlocksI, gllvmutils::dclist); //first entry is the number of species in colMat, rest are the blocks of colMat
  DATA_IMATRIX(nncolMat);
  DATA_VECTOR(Abranks);
  DATA_MATRIX(offset); //offset matrix
  DATA_IMATRIX(Ntrials);
  
  PARAMETER_MATRIX(r0f); // fixed site/row effects
  PARAMETER_MATRIX(r0r); // random site/row effects
  PARAMETER_MATRIX(b); // matrix of species specific intercepts and coefs
  // PARAMETER_MATRIX(bH); // matrix of species specific intercepts and coefs for beta hurdle model
  PARAMETER_MATRIX(B); // coefs of 4th corner model for TMBtrait, RE means for gllvm.TMB
  PARAMETER_MATRIX(Br); // random slopes for env species
  PARAMETER_MATRIX(b_lv); //slopes for RRR and constrained ord, VA means for random slopes
  //Left columns are for constrained ordination, Right for RRR
  PARAMETER_VECTOR(sigmaLV);//SD for LV
  PARAMETER_VECTOR(lambda); // lv loadings
  PARAMETER_MATRIX(lambda2);// quadratic lv loadings
  // PARAMETER_MATRIX(thetaH);// hurdle model lv loadings
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phiZINB);//extra param for ZINB
  PARAMETER_VECTOR(lg_phi); // dispersion params/extra zero probs for ZIP
  PARAMETER_VECTOR(sigmaB); // sds for random species effects
  PARAMETER_VECTOR(sigmab_lv); // sds for random slopes constr. ord.
  PARAMETER_VECTOR(sigmaij);// cov terms for random slopes covariance
  PARAMETER_VECTOR(log_sigma);// log(SD for row effect) and 
  PARAMETER_VECTOR(sigmaijr);// cors for row effect
  PARAMETER_MATRIX(rho_lvc);// correlation parameters for correlated LVs, matrix of q x 1 for corExp/corCS, qx2 for Matern, possibly q x times.cols() for corWithin
  
  DATA_INTEGER(num_lv); // number of lvs
  DATA_INTEGER(num_lv_c); //number of constrained lvs
  DATA_INTEGER(num_RR); //number of RRR dimensions
  DATA_INTEGER(num_corlv); //number of correlated lvs
  DATA_IVECTOR(family); // family index
  DATA_INTEGER(quadratic); // quadratic model, 0=no, 1=yes
  DATA_INTEGER(randomB); //0 = P,single,iid and 1 = LV
  PARAMETER_VECTOR(Au); // variational covariances for u
  PARAMETER_VECTOR(lg_Ar); // variational covariances for r0r
  PARAMETER_VECTOR(Abb);  // variational covariances for Br
  // PARAMETER_VECTOR(scaledc);// scale parameters for dc, of length of dc.cols()
  PARAMETER_VECTOR(Ab_lv); //variational covariances for b_lv
  PARAMETER_VECTOR(zeta); // ordinal family param

  PARAMETER(ePower);
  DATA_VECTOR(extra); // extra values, power of 
  DATA_INTEGER(method);// 0=VA, 1=LA, 2=EVA
  DATA_INTEGER(Abstruc); //0 = diagonal, blockdiagonal, 1 = MNdiagonal, MNunstructured, 2 = diagonalCL2, 3 = diagonalCL1, CL1, 4 = CL2, 5 = unstructured
  DATA_INTEGER(model);// which model, basic or 4th corner
  DATA_IVECTOR(random);//(0)1=random, (0)0=fixed row params, for Br: (1)1 = random slopes, (1)0 = fixed, for b_lv: (2)1 = random slopes, (2)0 = fixed slopes, for Br: (3) 1 = random
  DATA_INTEGER(zetastruc); //zeta param structure for ordinal model
  DATA_IMATRIX(trmsize); //2-row matrix. row 1: number of terms (LHS) in the random effect, row 2: number of groups (RHS) in  the random effect
  DATA_IMATRIX(csR); //2-column matrix. col 1: row number, col2: column number for correlation parameters of random row effects
  DATA_IMATRIX(times); //2 row matrix row 1: dim of LVs, row 2: LV number, used if multiple LVs of different structure/corwithin=TRUE
  DATA_IVECTOR(cstruc); //correlation structure for row.params 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm, 4= Matern
  DATA_STRUCT(proptoMats, gllvmutils::nesteddclist); //list of nested lists of length 2, first is the (inverse) matrix, second is the log determinant
  DATA_IVECTOR(cstruclv); //correlation structure for LVs 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm, 4= Matern
  DATA_STRUCT(dc, gllvmutils::dclist); //coordinates for sites, used for exponentially decaying cov. struc
  DATA_MATRIX(dc_lv); //coordinates for sites, used for exponentially decaying cov. struc
  DATA_INTEGER(Astruc); //Structure of the variational covariance, 0=diagonal, 1=RR, (2=sparse cholesky not implemented yet)
  DATA_IMATRIX(NN); //nearest neighbours,
  DATA_INTEGER(cw); //corWithin 0=FALSE, 1=TRUE,
  DATA_INTEGER(p_betaH); // number of bH columns, if non, zero
  
  int Klv = x_lv.cols();
  int n = y.rows();
  int p = y.cols();
  int truep = p-p_betaH; // For betaH
  // int nt =n;
  int nu =n; //CorLV
  
  // if(num_corlv>0){ //CorLV
  //   nu = dLV.cols();
  // }
  
  vector<Type> iphi = exp(lg_phi);
  
  // Set first row param to zero, if row effects are fixed
  // if(random(0)==0 && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){  r0f(0,0) = 0;}
  int nlvr = num_lv+num_lv_c;//treating constr. ord random slopes as a LV, to use existing infrastructure for integration
  
  matrix<Type> ucopy = u;
  // if(num_corlv>0){
  //   nlvr=0; num_lv=0; num_lv_c=0;
  //   quadratic=0;
  // }
  if(dLV.cols()>1){ //CorLV comb
    num_corlv = nlvr;
    u = (dLV*ucopy);
    // if(num_corlv>0){ //CorLV
    nu = dLV.cols();
    // }
    // nlvr=0; num_lv=0; num_lv_c=0;
    //Not yet combined with num_RR or quadratic
    quadratic=0;
    // num_RR=0;
  }
  
  // Distance matrix calculated from the coordinates for LVs
  matrix<Type> DiSc_lv(dc_lv.cols(),dc_lv.cols()); DiSc_lv.fill(0.0);
  matrix<Type> dc_scaled_lv(dc_lv.rows(),dc_lv.cols()); dc_scaled_lv.fill(0.0);
  // matrix<Type> DistM(dc.rows(),dc.rows());
  // if(((num_corlv>0) || (((random(0)>0) & (nlvr==(num_lv+num_lv_c))) & (rstruc>0))) & ((cstruc(0)==2) || (cstruc(0)>3))){
  //   matrix<Type> DiSc(dc.cols(),dc.cols());
  //   DiSc.setZero();
  // 
  //   for(int j=0; j<dc.cols(); j++){
  //     DiSc(j,j) += 1/exp(2*scaledc(j));
  //     // dc.col(j) *= 1/exp(scaledc(j));
  //   }
  //   // sigma_lvc(0,0) = 0;
  // 
  //   DistM.setZero();
  //   for (int d=0;d<dc.rows();d++) {
  //     for (int j=0;j<d;j++){
  //       DistM(d,j)=sqrt( ((dc.row(d)-dc.row(j))*DiSc*(dc.row(d)-dc.row(j)).transpose()).sum() ); // + extra(2);
  //   //     DistM(j,d)=DistM(d,j);
  //     }
  //   }
  // }
  
  matrix<Type> eta(n,p);
  eta.setZero();
  matrix<Type> lam(n,p);
  lam.setZero();
  // Type nll = 0;
  parallel_accumulator<Type> nll(this); // initial value of log-likelihood
  
  matrix<Type> RRgamma(num_RR,p);
  RRgamma.setZero();
  
  matrix <Type> Delta(nlvr,nlvr);
  Delta.setZero();
  
  matrix <Type> newlam(nlvr,p);
  newlam.setZero();  
  
  //K*K*d or d*d*K
  int sbl12 = Klv;
  int sbl3 = num_lv_c + num_RR;

  if(randomB>0){
    sbl12 = num_lv_c + num_RR;
    sbl3 = Klv;
  }
  
  vector<matrix<Type>> Sigmab_lv(sbl3);
  
  if(random(2)>0){
    for (int q=0; q<sbl3; q++){
      Sigmab_lv(q).resize(sbl12,sbl12);
      Sigmab_lv(q).setIdentity();
    }

    if(csb_lv.cols()<2){
    if(randomB<1){//Sigma_q = sigma_q I_klv
      //randomB="P","single","iid"
      if(sigmab_lv.size()>1){
      vector<Type>sigma1 = exp(sigmab_lv.head(x_lv.cols()));
      vector<Type>sigma2(num_lv_c+num_RR);
      sigma2.fill(1.0);
      sigma2.tail(num_lv_c+num_RR-1) = pow(exp(sigmab_lv.segment(x_lv.cols(), num_lv_c+num_RR-1)), 2);
      for (int q=0; q<(num_lv_c+num_RR); q++){
      Sigmab_lv(q).diagonal().array() = sigma1*sigma1;
      Sigmab_lv(q).diagonal().array() *= sigma2(q);
      }
      }else{
        for (int q=0; q<(num_lv_c+num_RR); q++){
          Sigmab_lv(q).diagonal().array() = exp(sigmab_lv(0))*exp(sigmab_lv(0));
        }
      }
      // matrix<Type> Sigmab_lvtemp(sbl12,sbl12);
      // Sigmab_lvtemp.setZero();
      // for (int q=0; q<sbl3; q++){
      //   Sigmab_lv(q) = Sigmab_lvtemp;
      //   Sigmab_lv(q).diagonal().array() = sigmab_lv(q);
      // }
    }else if(randomB>0){
      //randomB="LV"
        Sigmab_lv(0).diagonal().array() *= exp(2*sigmab_lv);
    }
    }else if((csb_lv.cols()==2) && (randomB<1)){
      matrix<Type> sds = Eigen::MatrixXd::Zero(x_lv.cols(),x_lv.cols());
      vector<Type>sigma1 = exp(sigmab_lv.head(x_lv.cols()));
      vector<Type>sigma2(num_lv_c+num_RR);
      sigma2.fill(1.0);
      sigma2.tail(num_lv_c+num_RR-1) = exp(sigmab_lv.segment(x_lv.cols(), num_lv_c+num_RR-1));
      sds.diagonal() = sigma1;
      vector<Type>corsb_lv((x_lv.cols()*x_lv.cols()-x_lv.cols())/2);
      corsb_lv.fill(0.0);
      matrix<Type>Sigmab_lvL(x_lv.cols(),x_lv.cols());
      Sigmab_lvL.setIdentity();
      if(csb_lv.cols()>1){
        //need a vector with covariances and zeros in the right places
        for(int i=0; i<csb_lv.rows(); i++){
          corsb_lv((csb_lv(i,0) - 1) * (csb_lv(i,0) - 2) / 2 + csb_lv(i,1)-1) = sigmab_lv(x_lv.cols()+num_lv_c+num_RR-1+i);
        }
        Sigmab_lvL = sds*gllvmutils::constructL(corsb_lv);
    }
      for (int q=0; q<(num_lv_c+num_RR); q++){
      Sigmab_lv(q) = Sigmab_lvL;
      Sigmab_lv(q) *= sigma2(q);
      }
    }else if((csb_lv.cols()==2) && (randomB>0)){
      Sigmab_lv.resize(num_lv_c+num_RR);//need to change this to same dimension as for randomB="P"
      
      for (int q=0; q<(num_lv_c+num_RR); q++){
        Sigmab_lv(q).resize(x_lv.cols(),x_lv.cols());
        Sigmab_lv(q).setIdentity();
      }
      // 
      // for (int q=0; q<(num_lv_c+num_RR); q++){
      //   Sigmab_lv(q).diagonal().array() = exp(sigmab_lv(q));
      // }
      // 
      vector<Type>corsb_lv((x_lv.cols()*x_lv.cols()-x_lv.cols())/2);
      corsb_lv.fill(0.0);
      matrix<Type>Sigmab_lvL(x_lv.cols(),x_lv.cols());
      Sigmab_lvL.setIdentity();
      if(csb_lv.cols()>1){
        //need a vector with covariances and zeros in the right places
        for(int i=0; i<csb_lv.rows(); i++){
          corsb_lv((csb_lv(i,0) - 1) * (csb_lv(i,0) - 2) / 2 + csb_lv(i,1)-1) = sigmab_lv(num_lv_c+num_RR+i);
        }
        Sigmab_lvL = gllvmutils::constructL(corsb_lv);
      }
        Sigmab_lv(0) = Sigmab_lvL;
    }
  }
  
  if((nlvr>0)||(num_RR>0)){
    
    if(nlvr>0){
      newlam.row(0).fill(1.0);
      if((num_lv+num_lv_c)>0){
        for (int d=0; d<nlvr; d++){
          Delta(d,d) = fabs(sigmaLV(d));
          // Delta(d,d) = fabs(sigmaLV(d));
        }
      }
    }
    //To create lambda as matrix Upper triangle
    // put LV loadings into a matrix
    if (num_lv>0){
      int tri = 0;
      if((num_lv_c+num_RR)>0){
        //because the lambdas for constrained and unconstrained LVs are separately identifiable and in the same vector
        tri += (num_lv_c+num_RR)*p-((num_lv_c+num_RR)*(num_lv_c+num_RR)-(num_lv_c+num_RR))/2-(num_lv_c+num_RR); //num_lv_c-1+p+(num_lv_c-1)*p-((num_lv_c-1)*(num_lv_c-1-1))/2-2*(num_lv_c-1); //number of elements for num_lv
      }
      for (int j=0; j<p; j++){
        for (int i=0; i<num_lv; i++){
          if(j<i){
            newlam(i+nlvr-num_lv-num_lv_c,j) = 0;
          }else if (j == i){
            newlam(i+nlvr-num_lv,j) = 1;
          }else if(j>i){
            newlam(i+nlvr-num_lv,j) = lambda(j+i*p-(i*(i-1))/2-2*i+tri-1);
          }
        }
      }
    }
    //species scores for constrained ordination and RRR
    if ((num_lv_c+num_RR)>0){
      for (int j=0; j<p; j++){
        for (int i=0; i<(num_lv_c+num_RR); i++){
          if(i<num_lv_c){
            if (j < i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = 0;
            } else if (j == i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = 1;
            }else if (j > i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = lambda(j+i*p-(i*(i-1))/2-2*i-1);//lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }else{
            if (j < i){
              RRgamma(i-num_lv_c,j) = 0;
            } else if (j == i){
              RRgamma(i-num_lv_c,j) = 1;
            }else if (j > i){
              RRgamma(i-num_lv_c,j) = lambda(j+i*p-(i*(i-1))/2-2*i-1);//lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }
          
        }
      }
    }
  }
  
  // Loadings for correlated latent variables //CorLV
  // matrix<Type> newlamCor;
  // matrix <Type> Delta_clv(num_corlv,num_corlv);
  // if((num_corlv)>0){
  //   newlamCor = matrix <Type> (num_corlv,p);
  //   //To create lambda as matrix Upper triangle
  //   // put LV loadings into a matrix
  //   for (int j=0; j<p; j++){
  //     for (int i=0; i<num_corlv; i++){
  //       if(j<i){
  //         newlamCor(i,j) = 0;
  //       }else if (j == i){
  //         newlamCor(i,j) = 1;
  //         // newlamCor(i,j) = exp(sigmaLV(i));
  //       }else if(j>i){
  //         newlamCor(i,j) = lambda(num_RR*p-num_RR*(num_RR+1)/2+j+i*p-(i*(i-1))/2-2*i-1);
  //       }
  //     }
  //   }
  //   for (int d=0; d<num_corlv; d++){
  //     // Delta_clv(d,d) = fabs(sigmaLV(d));
  //     newlamCor.row(d)*=fabs(sigmaLV(d));
  //     // newlamCor.row(d)*=exp(sigmaLV(d));
  //   }
  // }
  
  matrix<Type> mu(n,p);
  
  
  // Variational approximation
  if((method<1) || (method>1)){
    // add offset
    if(offset.rows()==n){
      eta += offset;
    }
    // add fixed row effects
    // if(r0f.size() == n && (random(0)==0) && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){
    //   eta += r0f.replicate(1,p);
    if(xr.rows()==n){
      eta += (xr*r0f).replicate(1,p);
    }
    
    matrix<Type> cQ(n,p);
    cQ.setZero();
    
    vector<matrix<Type>> A(n);
    
    if( (random(2)>0) && (num_RR>0) && (quadratic>0)){
      for(int i=0; i<n; i++){
        A(i).resize(nlvr+num_RR,nlvr+num_RR);
        A(i).setZero();
      }
    }else{
      for(int i=0; i<n; i++){
        A(i).resize(nlvr,nlvr);
        A(i).setZero();
      }
    }
    
    // lltOfB.matrixL() = A(0).template triangularView<Lower>;//wouuld be great if we could store A(i) each as a triangular matrix where the upper zeros are ignored
    // Set up variational covariance matrix for LVs 
    if((nlvr>0) & (num_corlv==0)){
      if((num_lv+num_lv_c)>0){
        // log-Cholesky parametrization for A_i:s
        // don't include num_RR for random slopes, comes in later
        for (int d=0; d<(num_lv+num_lv_c); d++){
          for(int i=0; i<n; i++){
            A(i)(d+(nlvr-num_lv-num_lv_c),d+(nlvr-num_lv-num_lv_c))=exp(Au(d*n+i));
            // A(d,d,i)=exp(Au(d*n+i));
          }
        }
        if(Au.size()>((num_lv+num_lv_c)*n)){
          int k=0;
          for (int c=0; c<(num_lv+num_lv_c); c++){
            for (int r=c+1; r<(num_lv+num_lv_c); r++){
              for(int i=0; i<n; i++){
                A(i)(r+(nlvr-num_lv-num_lv_c),c+(nlvr-num_lv-num_lv_c))=Au((num_lv+num_lv_c)*n+k*n+i);
                // A(r,c,i)=Au(nlvr*n+k*n+i);
                // A(c,r,i)=A(r,c,i);
              }
              k++;
            }}
        }
      }
      
      // Add VA terms to logL
      //Go this route if no random Bs
      matrix <Type> Atemp(nlvr,nlvr);
      for(int i=0; i<n; i++){
        Atemp = A(i).topLeftCorner(nlvr,nlvr);//to exlcude the 0 rows & columns for num_RR
        nll -= Atemp.diagonal().array().log().sum() - 0.5*((Atemp*Atemp.transpose()).trace()+(u.row(i)*u.row(i).transpose()).sum());
      }
      nll -= 0.5*n*nlvr;
      
      //scale LVs with standard deviations, as well as the VA covariance matrices
      u *= Delta;
      
      if((num_RR*random(2))>0 && (quadratic)>0){
        Delta.conservativeResize(nlvr+num_RR,nlvr+num_RR);
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          Delta.col(d).setZero();
          Delta.row(d).setZero();
        }
      }
      
      for (int i=0; i<n; i++) {
        A(i) = Delta*A(i);
      }
      
    } else if(num_corlv > 0) {
      u *= Delta;
      if((num_RR*random(2))>0 && (quadratic)>0){
        Delta.conservativeResize(nlvr+num_RR,nlvr+num_RR);
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          Delta.col(d).setZero();
          Delta.row(d).setZero();
        }
      }
    } // ad else for num_corlv to create ucopy & create D*u*= Delta;
    
    //random slopes for constr. ord.
    vector<matrix<Type>> Ab_lvcov;  //covariance of LVs due to random slopes
    if((random(2)>0) && ((num_RR+num_lv_c)>0)){
      // Variational covariance for random slopes
      vector<matrix<Type>> AB_lv(sbl3);
      for(int d=0; d<sbl3; d++){
        AB_lv(d).resize(sbl12,sbl12);
        AB_lv(d).setZero();
      }
      
      for(int d=0; d<sbl3; d++){
        for (int q=0; q<(sbl12); q++){
          AB_lv(d)(q,q)=exp(Ab_lv(q*sbl3+d));
        }
      }
      
      if(Ab_lv.size()>((sbl12)*sbl3)){
        int k=0;
        for (int c=0; c<(sbl12); c++){
          for (int r=c+1; r<(sbl12); r++){
            for(int d=0; d<sbl3; d++){
              AB_lv(d)(r,c)=Ab_lv((sbl12)*sbl3+k*sbl3+d);
            }
            k++;
          }}
      }

      //VA likelihood parts for random slope
      if(csb_lv.cols()<2){
      //randomB and no correlation
      for(int q=0; q<sbl3; q++){
        if(randomB<1){
          nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(Sigmab_lv(q).diagonal().cwiseInverse().array()*AB_lv(q).rowwise().squaredNorm().array()).sum()-0.5*(b_lv.col(q).transpose()*Sigmab_lv(q).diagonal().cwiseInverse().asDiagonal()*b_lv.col(q)).value());
          nll -= 0.5*(sbl12- Sigmab_lv(q).diagonal().array().log().sum());
        }
        if(randomB>0)nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(Sigmab_lv(0).diagonal().cwiseInverse().array()*AB_lv(q).rowwise().squaredNorm().array()).sum()-0.5*(b_lv.row(q)*Sigmab_lv(0).diagonal().cwiseInverse().asDiagonal()*b_lv.row(q).transpose()).value());// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
       }
      if(randomB>0)nll -= 0.5*(sbl3*sbl12- sbl3*Sigmab_lv(0).diagonal().array().log().sum());
      }else if((csb_lv.cols()==2) && (randomB<1)){
        //randomB="P" with predictor-wise correlation
          matrix<Type>Iblv = Eigen::MatrixXd::Identity(x_lv.cols(),x_lv.cols());
          matrix<Type>Sigmab_lvI = (Sigmab_lv(0)*Sigmab_lv(0).transpose()).ldlt().solve(Iblv);
          vector<Type>sigma2(num_lv_c+num_RR);
          sigma2.fill(1.0);
          sigma2.tail(num_lv_c+num_RR-1) = exp(-2*sigmab_lv.segment(x_lv.cols(), num_lv_c+num_RR-1));
          
          for(int q=0; q<(num_lv_c+num_RR); q++){
          nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(sigma2(q)*Sigmab_lvI*AB_lv(q)*AB_lv(q).transpose()).trace()-0.5*(b_lv.col(q).transpose()*(sigma2(q)*Sigmab_lvI)*b_lv.col(q)).sum());
          nll -= 0.5*Klv- Sigmab_lv(q).diagonal().array().log().sum();//Sigmab_lv already includes sigma2
        }
        
      }else if((csb_lv.cols()==2) && (randomB>0)){
        //randomB="LV" with predictor-wise correlation
        matrix<Type>Iblv = Eigen::MatrixXd::Identity(x_lv.cols(),x_lv.cols());
        //note that Sigmab_lv(0) is the cholesky of the correlation matrix
        matrix<Type>Sigmab_lvCI = Sigmab_lv(0).template triangularView<Eigen::Lower>().solve(Iblv);
        Sigmab_lvCI = Sigmab_lvCI.transpose() * Sigmab_lvCI;//inverse of correlation matrix via its cholesky
        for(int q=0; q<sbl3; q++){
          nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(exp(-2*sigmab_lv.head(num_lv_c+num_RR)).array()*Sigmab_lvCI(q,q)*(AB_lv(q)*AB_lv(q).transpose()).diagonal().array()).sum());//need to use sigmab_lv directly here, as Sigmab_lv is now of length x_lv.cols() for the correlation
        }
        
        for(int q=0; q<(num_lv_c+num_RR); q++){
          nll -= -0.5*(b_lv.col(q).transpose()*Sigmab_lvCI*b_lv.col(q)).sum()*exp(-2*sigmab_lv(q));
          nll -= 0.5*Klv- Klv*sigmab_lv(q)-Sigmab_lv(0).diagonal().array().log().sum();
        }
      }
      
      //resize ab_lvcov to correct size
      Ab_lvcov  = vector<matrix<Type>> (n);
      for(int i=0; i<n; i++){
        Ab_lvcov(i).resize(num_RR+num_lv_c,num_RR+num_lv_c);
        Ab_lvcov(i).setZero();
      }
      
      //fill ab_lvcov
      if(randomB>0){//variance per LV
        for(int klv=0; klv<Klv; klv++){
        matrix<Type>Ablv = AB_lv(klv)*AB_lv(klv).transpose();
        for(int i=0; i<n; i++){
            Ab_lvcov(i) += x_lv(i,klv)*x_lv(i,klv)*Ablv;
          }
        }
      }else{
      for(int q=0; q<(num_RR+num_lv_c); q++){
        matrix<Type>Ablv = AB_lv(q)*AB_lv(q).transpose();
        for(int i=0; i<n; i++){
            Ab_lvcov(i)(q,q) = (x_lv.row(i)*Ablv*x_lv.row(i).transpose()).sum();
          }
        }
      }
      
      if(num_lv_c>0){
        RRgamma.conservativeResize(num_RR+num_lv_c,Eigen::NoChange);
        if(num_RR>0)RRgamma.bottomRows(num_RR) = RRgamma.topRows(num_RR); 
        RRgamma.topRows(num_lv_c) = newlam.topRows(num_lv_c);
      }
      for(int i=0; i<n; i++){
        matrix<Type> Av = Ab_lvcov(i)*RRgamma; // reuse Ab_lvcov(i) across all species
        for(int j=0; j<p; j++){
          cQ(i,j) += 0.5*(RRgamma.col(j).transpose()*Av.col(j)).value();
        }
      }
      if(quadratic<1){
        eta += x_lv*b_lv*RRgamma;//for the quadratic model this component is added below
        
      }else if(quadratic>0){
        //now rebuild A and u with covariances for random slopes so that existing infrastructure below can be used
        //in essence, q(XBsigmab_lv + eDelta) ~ N(uDelta + \sum \limits^K X_ik b_lv_k , Delta A Delta + \sum \limits^K X_ik^2 AB_lv_k )
        //so build u and A accordingly (and note covariance due to Bs if num_lv_c and num_RR > 0)
        if(num_lv_c>0){
          u.leftCols(num_lv_c) += x_lv*b_lv.leftCols(num_lv_c);
        }
        //resize A, u, D, and add RRGamma to newlam.
        //add columns to u on the right for num_RR with random slopes
        if(num_RR>0){
          u.conservativeResize(n, nlvr + num_RR);
          //resize and fill newlam, we don't use RRgamma further with random Bs
          //easiest to do is slap RRgamma at the end of newlam
          //this makes the order of newlam, A, u, and D inconsistent with the R-side of things
          //nicer would be to have to same order as in R, but that isn't possible since
          //it requires going down the same route for fixed and random B
          //which would only work with diagonal of 0s in A
          //And that needs to be invertible for the quadratic case, so that is not possible
          newlam.conservativeResize(nlvr+num_RR,p);
          for(int d=nlvr; d<(nlvr+num_RR); d++){
            u.col(d).fill(0.0);
            newlam.row(d).fill(0.0);
          }
          nlvr += num_RR;
          newlam.bottomRows(num_RR) = RRgamma.bottomRows(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv.rightCols(num_RR);
          REPORT(RRgamma);
        }
        
        if((nlvr-num_RR-num_lv_c)>0){
          //rebuild Ab_lvcov to fit A below
          matrix<Type> tempRR(num_RR,num_RR);
          matrix<Type> tempCN(num_lv_c,num_lv_c);
          matrix<Type> tempRRCN(num_RR,num_lv_c);
          
          for(int i=0; i<n; i++){
            if(num_RR>0){
              tempRR = Ab_lvcov(i).bottomRightCorner(num_RR,num_RR);
            }
            if(num_lv_c>0){
              tempCN = Ab_lvcov(i).topLeftCorner(num_lv_c,num_lv_c);
            }
            if((num_RR+num_lv_c)>0){
              tempRRCN = Ab_lvcov(i).bottomLeftCorner(num_RR,num_lv_c);
            }
            //resize to fit A
            Ab_lvcov(i).resize(nlvr,nlvr);
            Ab_lvcov(i).setZero();
            
            //re-assign
            //place num_RR in back
            if(num_RR>0){
              Ab_lvcov(i).bottomRightCorner(num_RR,num_RR) = tempRR;
            }
            //num_lv_c is in front, but after a potential intercept
            if((num_lv_c)>0){
              Ab_lvcov(i).topLeftCorner(num_lv_c,num_lv_c) = tempCN;
            }
            //assign covariances of random slopes. There is no covariance if slb3==num_lv_c+num_RR, randomB == 0
            if((num_RR>0)&&(num_lv_c>0)&&(randomB>0)){
              Ab_lvcov(i).block(nlvr-num_RR,nlvr-num_lv_c-num_RR-num_lv,num_RR,num_lv_c) = tempRRCN;
              Ab_lvcov(i).block(nlvr-num_lv_c-num_RR-num_lv,nlvr-num_RR,num_lv_c,num_RR) = tempRRCN.transpose();
            }
            
          }
        }
      }
      
    }
    
    //components for reduced rank regression terms
    if((num_RR>0) && (random(2)<1)){
      //predictor coefficients RRR.  num_RR comes after num_lv_c
      //Since later dimensions are more likely to have less residual variance
      eta += x_lv*b_lv.rightCols(num_RR)*RRgamma;
      
      //quadratic terms for fixed-effects only RRR
      //-num_lv to ensure that we pick num_RR from the middle
      if(quadratic>0){
        matrix<Type> D_RR(num_RR,num_RR);
        D_RR.setZero();
        
        //quadratic coefficients for RRR
        if(lambda2.cols()==1){
          for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
            D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,0));
          }
          // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            eta.row(i).array() -=  (x_lv.row(i)*b_lv.rightCols(num_RR)*D_RR*(x_lv.row(i)*b_lv.rightCols(num_RR)).transpose()).value();
          }
          // }
        }else{
          for (int j=0; j<p;j++){
            D_RR.setZero();
            for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
              D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,j));
            }
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv.rightCols(num_RR)*D_RR*(x_lv.row(i)*b_lv.rightCols(num_RR)).transpose();
            }
            
          }
        }
        
      }
    }
    
    #include "species_effects.h"
    
    
    if(model<1){
      // basic gllvm, gllvm.TMB.R
      eta += x*b;
    } else {
      // Fourth corner model TMB.trait.R
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)*extra(p)+eta1(m,0); //extra(p)=0 if beta0comm=TRUE
          m++;
        }
      }
    }
    
    #include "row_effects.h"
    
    
    vector<Eigen::DiagonalMatrix<Type, Eigen::Dynamic>> D(p);
    
    // LVs in model:
    if(nlvr>0){
      matrix<Type> b_lv2(x_lv.cols(),nlvr);
      b_lv2.setZero();
      
      if((num_lv_c>0) && (random(2)<1)){
        //concurrent ordination terms
        //predictor coefficients for constrained ordination
        b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        eta += x_lv*b_lv2*newlam;
        
      }else if((nlvr>0) && (random(2)>0) && (quadratic > 0)){
        if(num_lv_c>0)b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        if(num_RR>0) b_lv2.rightCols(num_RR) = b_lv.rightCols(num_RR);
      }
      lam = u*newlam;
      
      // Update cQ for linear term
      //Binomial, Gaussian, Ordinal
      if(num_corlv==0){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            cQ(i,j) += 0.5*(newlam.col(j).transpose()*A(i)*A(i).transpose()*newlam.col(j)).value();
          }
        }
      } else if(num_corlv>0) { //CorLV // Correlated LVs
        int i,j,d;
        int arank = 2;
        matrix<Type> AQ(num_corlv,num_corlv);
        AQ.setZero(); AQ.diagonal().fill(1.0);
        
        // matrix <Type> newlamCor(num_corlv,p);
        // for (int d=0; d<num_corlv; d++){
        //   newlamCor.row(d)=newlam.row(d);
        //   // newlamCor.row(d)=newlam.row(d)*fabs(sigmaLV(d));
        // }
        // REPORT(nu);
        if(cw == 0){
            // eta += (dLV*ucopy)*newlamCor;
          matrix<Type> AAT; 
          
          if(cstruclv(0)==0){
            matrix<Type> DAATD;
            matrix<Type> Alvm(num_corlv,num_corlv);
            
            // Variational covariance: diagonal
            for (int d=0; d<(nu); d++){
              Alvm.setZero(num_corlv,num_corlv);
              for (int q=0; q<(num_corlv); q++){
                Alvm(q,q)=exp(Au(q*nu+d));
              }
            
              // Off diagonal
              if((Astruc>0) & (Au.size()>((num_corlv)*nu))){//unstructured cov
                int k=0;
                for (int c=0; c<(num_corlv); c++){
                  for (int r=c+1; r<(num_corlv); r++){
                      Alvm(r,c)=Au(nu*num_corlv+k*nu+d);
                    k++;
                  }}
              }
            
              nll -= Alvm.diagonal().array().log().sum() - (Alvm.array().square()).sum() - 0.5*((ucopy.row(d)*ucopy.row(d).transpose()).sum());  // nll -= atomic::logdet(Alvm.col(d).matrix()) + 0.5*( - (Alvm.col(d).matrix()*Alvm.col(d).matrix().transpose()).diagonal().sum() - (ucopy.row(d).matrix()*ucopy.row(d).matrix().transpose()).sum());
              
              AAT = Alvm*Alvm.transpose();
              DAATD = Delta * AAT * Delta.transpose();
              for (j=0; j<p;j++){
                cQ.col(j) += 0.5*dLV.col(d)*((newlam.col(j).transpose()*DAATD)*newlam.col(j));
              }
            }
            nll -= 0.5*(nu*num_corlv);
            
          } else {
            vector<matrix<Type> > Slv(num_corlv);
            
            matrix<Type> Slvinv;
            
            for(int q=0; q<num_corlv; q++){
              // site specific LVs, which are correlated between sites/groups
              if(cstruclv(0)==1){// AR1 covariance
                Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
              } else if(cstruclv(0)==3) {// Compound Symm  if(cstruclv==3)
                Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
              } else {
                
                // Slv(q) = exp(-dc_lv.array()*Type(1/exp(rho_lvc(q,0))) ).matrix()*Type(0.99);
                // Slv(q).diagonal().fill(1.0);
                DiSc_lv.fill(0.0);
                for(int j=0; j<dc_lv.cols(); j++){
                  DiSc_lv(j,j) += 1/exp(rho_lvc(q,0));
                }
                dc_scaled_lv = dc_lv*DiSc_lv;
                if(cstruclv(0)==2){// exp decaying
                  Slv(q) = gllvm::corExp(Type(1), Type(0), nu, dc_scaled_lv);
                } else if(cstruclv(0)==4) {// Matern
                  Slv(q) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,rho_lvc.cols()-1)), nu, dc_scaled_lv);
                }
              }
              nll -= 0.5*(nu - atomic::logdet(Slv(q)));
            }
            
            if(Astruc<3){
              
              for(int q=0; q<num_corlv; q++){
                
                // u^T*Sinv*u
                Slvinv = atomic::matinv(Slv(q));
                nll -= - 0.5*( ucopy.col(q).transpose()*(Slvinv*ucopy.col(q)) ).sum();
                
                if(Astruc==0 ){//diagonal A cov
                  vector<Type> Atemp = exp(Au.segment(q*nu, nu));
                  
                  vector<Type> AtempSq(nu);
                  for(int d=0; d<nu; d++) AtempSq[d] = Atemp[d] * Atemp[d];
                  
                  // summa tr(Sinv * A) = sum_i Sinv(ii) * AtempSq(i)
                  Type trSinvA = (Slvinv.diagonal().array() * AtempSq.array()).sum();
                  nll -= -0.5 * trSinvA;
                  
                  // for (d=0; d<(nu); d++){ // - tr(Sinv*A)
                  //   nll -= - 0.5*Slvinv(d,d)*pow(Atemp(d),2);
                  // }
                  
                  // 0.5*lambda_qj*A_qii*lambda_qj
                  for (j=0; j<p;j++){
                    Type sca = 0.5*pow(newlam(q,j),2)*pow(Delta(q,q),2);
                    cQ.col(j) += sca*(dLV*AtempSq.matrix());
                  }
                  // 0.5*logdet(A)
                  nll -= Atemp.array().log().sum();

                } else if((Astruc>0)){
                  matrix<Type> Atemp(nu, nu);
                  Atemp.setZero();
                  
                  // diagonal
                  Atemp.diagonal().array() = (Au.segment(q*nu, nu)).array().exp();
                  
                  int k=0;
                  if((Astruc==1) & (Au.size() > nu*num_corlv) ){ // unstructured variational covariance
                    for (d=0; d<nu; d++){
                      for (int r=d+1; r<(nu); r++){
                        Atemp(r,d)=Au(nu*num_corlv+k*num_corlv+q);
                        k++;
                      }
                    }
                  } else if((Astruc==2) & (Au.size() > nu*num_corlv)) { // bdNN variational covariance
                    arank = NN.rows();
                    for (int r=0; r<(arank); r++){
                      Atemp(NN(r,0)-1,NN(r,1)-1)=Au(nu*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                  // REPORT(k);
                  
                  // 0.5*lambda_qj*A_qii*lambda_qj
                  AAT = Atemp*Atemp.transpose();
                  for (j=0; j<p;j++){
                    Type sca = 0.5*pow(newlam(q,j),2)*pow(Delta(q,q),2);
                    cQ.col(j) += sca*(dLV*AAT.diagonal().matrix()); //this works
                  }
                  
                  // 0.5*logdet(A) -0.5*tr(Sinv*A)
                  nll -= Atemp.diagonal().array().log().sum() + 0.5*(- (Slvinv*AAT).diagonal().sum());
                }
                
              }
              
            } else if((num_corlv>1) & (Astruc<6)){
              // UNN/Kronecker variational covariance
              matrix<Type> Alvm = matrix<Type>::Zero(nu, nu);
              
              // diagonal
              Alvm.diagonal().array() = (Au.segment(0, nu)).array().exp();
              
              int k=0;
              arank = NN.rows();
              if(Au.size()>(nu+num_corlv*(num_corlv+1)/2)) {
                if(Astruc == 4) {
                  for (int r=0; r<(arank); r++){
                    Alvm(NN(r,0)-1,NN(r,1)-1)=Au(nu+k);
                    ++k;
                  }
                } else if(Astruc == 3) {
                  for (int i = 1; i < nu; ++i) {
                    for (int j = 0; j < i; ++j) {
                      Alvm(i, j) = Au(nu + k);
                      ++k;
                    }
                  }
                }
              }
              
              for (d=0; d<num_corlv; d++){
                AQ(d,d)=exp(Au(nu+k));
                ++k;
                for (int r=d+1; r<num_corlv; r++){
                  AQ(r,d)=Au(nu+k);
                  ++k;
                }
              }
              
              // logdet(A) for triang.mat = prod of diag. elements
              // Moved right after Slv initialization: + 0.5*num_corlv*nu;
              const Type logdet_Alvm = Alvm.diagonal().array().log().sum();
              const Type logdet_AQ   = AQ.diagonal().array().log().sum();
              nll -= num_corlv * logdet_Alvm + nu * logdet_AQ;
              // nll -= num_corlv*Alvm.diagonal().array().log().sum() + nu*AQ.diagonal().array().log().sum();
              // nll -= num_corlv*log(Alvm.determinant()) + nu*log(AQ.determinant()) + 0.5*num_corlv*nu;
              
              // Alvm *= Alvm.transpose();
              // AQ *= AQ.transpose();
              Alvm = Alvm*Alvm.transpose();
              AQ = AQ*AQ.transpose();
              
              
              // tr(Sinv*A) + u^T*Sinv*u
              for(int q=0; q<num_corlv; q++){
                const matrix<Type>& Slvinv = atomic::matinv(Slv(q));
                nll -= 0.5*(- AQ(q,q)*(Slvinv*Alvm).trace()-( ucopy.col(q).transpose()*(Slvinv*ucopy.col(q)) ).sum());
              }
              
              matrix<Type> AQt = Delta * AQ * Delta.transpose();
              // 0.5*lambda_qj*A_qii*lambda_qj
              for (j=0; j<p;j++){
                Type sca = 0.5 * (newlam.col(j).transpose() * AQt * newlam.col(j)).sum();
                cQ.col(j) += sca * (dLV * Alvm.diagonal());
                // cQ.col(j) += 0.5*(dLV*Alvm.diagonal())*((newlam.col(j).transpose()*(Delta*AQ*Delta.transpose()))*newlam.col(j));
              }
              
              REPORT(Alvm);
            }
            
          }
        } else {
          // Correlation within group
          // eta += ucopy*newlamCor;
          
          nu = times.row(0).size();
          int it_ind = 0;
          int nt = times.row(0).sum();
          vector<matrix<Type>> Slv(nu);        // Cor matrix
          // vector<matrix<Type>> Alvm(num_corlv);
          vector<matrix<Type>> Alvm(nu);
          matrix<Type> Slvinv;
          matrix<Type> AlvAlvT; 
          
          for (int i = 0; i < nu; i++) {
            Slv(i).setZero(times(0,i), times(0,i));
            // Alvm(i).setZero(times(i), times(i));
            Alvm(i) = matrix<Type>::Zero(times(0,i), times(0,i));
          }
          
          if(Astruc<3){
            for(int q=0; q<num_corlv; q++){
              
              // site specific LVs, which are correlated within groups
              
              // Variational covariance for row effects
              //diagonal
              it_ind = 0;
              for (int i = 0; i < nu; i++) {
                Alvm(i).setZero();                                 // nollataan vain tarvittaessa
                Alvm(i).diagonal().array() = (Au.segment(q*nt + it_ind, times(0,i))).array().exp();
                it_ind += times(0,i);
              }
              
              if((Astruc>0) && (Au.size() > nt*num_corlv)){//reduced rank cov
                int k=0;
                it_ind = 0;
                
                if(Astruc==1){
                  for(i=0; i<nu; i++){
                    for (d=0; (d<times(0,i)); d++){
                      for (int r=d+1; r<(times(0,i)); r++){
                        // Alvm(q)(it_ind+r,it_ind+d)=Au(nt*num_corlv+k*num_corlv+q);
                        Alvm(i)(r,d)=Au(nt*num_corlv+k*num_corlv+q);
                        k++;
                      }
                    }
                    it_ind += times(0,i); 
                  }
                } else if(Astruc==2) { //bdNN var cov
                  arank = NN.rows();
                    for (int r=0; r<(arank); r++){
                      Alvm(NN(r,2)-1)(NN(r,0)-1,NN(r,1)-1)=Au(nt*num_corlv+k*num_corlv+q);
                      k++;
                    }
                }
              }
              
              
              it_ind = 0;
              for(i=0; i<nu; i++){
                // Compute Alvm*Alvm'
                AlvAlvT = Alvm(i) * Alvm(i).transpose();
                
                // Update cQ with 0.5*gamma'A gamma
                for (j=0; j<p;j++){
                  Type sca = 0.5*pow(newlam(q,j),2)*pow(Delta(q,q),2);
                  cQ.col(j) += sca*(dLV.block(0,it_ind, dLV.rows(), times(0,i))*AlvAlvT.diagonal());
                }
                nll -= Alvm(i).diagonal().array().log().sum(); //log(Alvm(i).determinant());
                
                Slv(i).setZero();
                
                // Define covariance matrix
                int ics =0;
                if(cstruclv.size() >= nu) ics =i;
                if(cstruclv(ics)==1){// AR1 covariance
                  Slv(i) = gllvm::corAR1(Type(1), rho_lvc(q,i), times(0,i));
                } else if(cstruclv(ics)==3) {// Compound Symm  if(cstruclv==3)
                  Slv(i) = gllvm::corCS(Type(1), rho_lvc(q,i), times(0,i));
                } else {
                  DiSc_lv.setZero();
                  for(int j=0; j<dc_lv.cols(); j++){
                    DiSc_lv(j,j) += 1/exp(rho_lvc(q,i));
                  }
                  dc_scaled_lv = dc_lv.block(it_ind,0,times(0,i),dc_lv.cols())*DiSc_lv;
                  if(cstruclv(ics)==2){// exp decaying
                    Slv(i) = gllvm::corExp(Type(1), Type(0), times(0,i), dc_scaled_lv);
                  } else if(cstruclv(ics)==4) {// matern
                    Slv(i) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,rho_lvc.cols()-1)), times(0,i), dc_scaled_lv);
                  }
                }
                
                nll -= 0.5*(times(0,i) - atomic::logdet(Slv(i)));

                Slvinv = atomic::matinv(Slv(i));
                matrix<Type> ublock = ucopy.block(it_ind,q,times(0,i),1);
                nll -= 0.5*(- (Slvinv * AlvAlvT).trace()-(ublock.transpose() * Slvinv * ublock).sum());
                it_ind += times(0,i);
                
              }
            }
            
            // REPORT(Alvm);
          } else if(num_corlv>1){
            // Kron A=AQ*Alvm
            // Variational covariance
            
            //diagonal
            it_ind = 0;
            for(i=0; i<nu; i++){
              // Alvm(i).setZero();
              // Alvm(i).diagonal()=exp(Au.segment(it_ind, times(i)));
              Alvm(i).diagonal().array() = (Au.segment(it_ind, times(0,i))).array().exp();
              it_ind += times(0,i); 
            }
            
            int k=0;
            it_ind = 0;
            arank = NN.rows();
            if(Au.size()>(nt+num_corlv*(num_corlv+1)/2)) {
              if(Astruc == 4) {
                  for (int r=0; r<(arank); r++){
                    Alvm(NN(r,2)-1)(NN(r,0)-1,NN(r,1)-1)=Au(nt+k);
                    k++;
                  }
              } else if(Astruc == 3){
                for(i=0; i<nu; i++){
                  for (d=0; (d<times(0,i)); d++){
                    for (int r=d+1; r<(times(0,i)); r++){
                      Alvm(i)(r,d)=Au(nt+k);
                      k++;
                    }
                  }
                }
                
              }
            }

            for (d=0; d<num_corlv; d++){
              AQ(d,d)=exp(Au(nt+k));
              k++;
              for (int r=d+1; r<(num_corlv); r++){
                AQ(r,d)=Au(nt+k);
                k++;
              }
            }
            matrix<Type> AQAQT = AQ *AQ.transpose();
            
            it_ind = 0;
            for(i=0; i<nu; i++){
              // Compute Alvm*Alvm'
              AlvAlvT = Alvm(i) * Alvm(i).transpose();
              
              matrix<Type> DAQD = (Delta*AQAQT*Delta.transpose());
              for (j=0; j<p;j++){
                Type sca = 0.5 * (newlam.col(j).transpose() * DAQD * newlam.col(j)).sum();
                cQ.col(j) += sca*(dLV.block(0,it_ind, dLV.rows(), times(0,i))*AlvAlvT.diagonal());
              }
              
              nll -= num_corlv*Alvm(i).diagonal().array().log().sum() + times(0,i)*AQ.diagonal().array().log().sum(); //log(Alvm(i).determinant());
              
              for(int q=0; q<num_corlv; q++){
                Slv(i).setZero();

                int ics =0;
                if(cstruclv.size() >= nu) ics =i;
                // Define covariance matrix
                if(cstruclv(ics)==1){// AR1 covariance
                  Slv(i) = gllvm::corAR1(Type(1), rho_lvc(q,i), times(0,i));
                } else if(cstruclv(ics)==3) {// Compound Symm  if(cstruclv==3)
                  Slv(i) = gllvm::corCS(Type(1), rho_lvc(q,i), times(0,i));
                } else {
                  DiSc_lv.setZero();
                  for(int j=0; j<dc_lv.cols(); j++){
                    DiSc_lv(j,j) += 1/exp(rho_lvc(q,i));
                  }
                  dc_scaled_lv = dc_lv.block(it_ind,0,times(0,i),dc_lv.cols())*DiSc_lv;
                  if(cstruclv(ics)==2){// exp decaying
                    Slv(i) = gllvm::corExp(Type(1), Type(0), times(0,i), dc_scaled_lv);
                  } else if(cstruclv(ics)==4) {// matern
                    Slv(i) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,rho_lvc.cols()-1)), times(0,i), dc_scaled_lv);
                  }
                }
                
                nll -= 0.5*(times(0,i) - atomic::logdet(Slv(i)));
                // nll -= - 0.5*nu*atomic::logdet(Slv(i));
                Slvinv = atomic::matinv(Slv(i));
                
                matrix<Type> ublock = ucopy.block(it_ind,q,times(0,i),1);
                nll -=  0.5*(- AQAQT(q,q)*(Slvinv*AlvAlvT).trace() - (ublock.transpose()*(Slvinv*ublock)).sum());
                
              }
              it_ind += times(0,i); 
            }
            REPORT(Alvm);
          }
          
          
        }
        REPORT(AQ);
        // REPORT(newlam);
      }
      // REPORT(newlam);
      // REPORT(A);
      // REPORT(u);
      // REPORT(ucopy);
      // REPORT(Delta);
      
  
      eta += lam;
      
      if(((quadratic>0) && (nlvr>0)) || ((quadratic>0) && (num_RR>0))){
        //quadratic coefficients for ordination
        //if random rows, add quadratic coefficients for num_RR to D otherwise
        //they go into D_RR below
        //The ordering here is num_lv_c-num_lv-num_RR so that the code works for
        //fixed-effects B and random effects B
        //The order we need to pick them from lambda2 is
        //num_lv_c-num_RR-num_lv however, to ensure everything on the R-side works
        if(((num_lv+num_lv_c+num_RR*random(2))>0)){
          
          for (int j=0; j<p; j++){
            D(j).resize(nlvr);
            D(j).setZero();
          }
          
          if(num_lv_c>0){
            if(lambda2.cols()==1){
              for (int j=0; j<p; j++){
                for (int q=0; q<num_lv_c; q++){
                  D(j).diagonal()(q) = fabs(lambda2(q,0)); //common tolerances model
                }
              }
            }else{
              for (int j=0; j<p; j++){
                for (int q=0; q<num_lv_c; q++){
                  D(j).diagonal()(q) = fabs(lambda2(q,j)); //full quadratic model
                }
              }
            }
          }
          if((num_RR*random(2))>0){
            if(lambda2.cols()==1){
              //make sure that num_RR comes at the end..has to be
              //like this due to the difference between fixed and random Bs
              for (int j=0; j<p; j++){
                for (int q=(num_lv+num_lv_c); q<nlvr; q++){
                  D(j).diagonal()(q) = fabs(lambda2(q-num_lv,0)); //common tolerances model
                }
              }
            }else{
              for (int j=0; j<p; j++){
                for (int q=(num_lv+num_lv_c); q<nlvr; q++){
                  D(j).diagonal()(q) = fabs(lambda2(q-num_lv,j)); //full quadratic model
                }
              }
            }
          }
          if(num_lv>0){
            if(lambda2.cols()==1){
              //make sure that num_lv is taken from the middle even with num_RR
              for (int j=0; j<p; j++){
                for (int q=(num_lv_c+num_RR*random(2)); q<(num_lv_c+num_RR*random(2)+num_lv); q++){
                  D(j).diagonal()(q-num_RR*random(2)) = fabs(lambda2(q,0)); //common tolerances model
                }
              }
            }else{
              for (int j=0; j<p; j++){
                for (int q=(num_lv_c+num_RR*random(2)); q<(num_lv_c+num_RR*random(2)+num_lv); q++){
                  D(j).diagonal()(q-num_RR*random(2)) = fabs(lambda2(q,j)); //full quadratic model
                }
              }
            }
          }
        }
        if((num_lv_c>0) && (random(2)<1)){
          //quadratic reduced rank term for concurrent ordination
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv2*D(j)*(x_lv.row(i)*b_lv2).transpose();
            }
          }
        }
        
        // do not take this route not with quadratic model, (fixed-effect) constrained LVs and random row-effects.
        if(((nlvr > 0) && (num_lv+num_lv_c)>0) || ((quadratic>0) && (random(2) > 0))){
          //quadratic model approximation
          
          matrix <Type> Acov(nlvr,nlvr);
          matrix<Type> Id(nlvr,nlvr);
          Id.setIdentity();
          
          matrix<Type> B(nlvr,nlvr);
          matrix<Type> Binv(nlvr,nlvr);
          
          for (int j = 0; j < p; j++) {
            
            // group Poisson, ZIP, Tweedie, NB, gamma, exponential
            bool is_group1 =
              (family(j)==0 || family(j)==1 || family(j)==4 ||family(j)==5 || 
              family(j)==6 || family(j)==8 || family(j)==11);
            
            // group Binomial/Gaussian/Ordinal
            bool is_group2 =
              (family(j)==2 || family(j)==3 || family(j)==7);
            
            // sign only if Poisson/NB/gamma/exponential/ZIP/Tweedie
            int sign = 0;
            if (is_group1) {
              if ((family(j) > 0) && (family(j) != 6) && (family(j) != 5))
                sign = 1;   // NB, gamma, exponential, ZIP
              else
                sign = -1;  // Poisson, ZIP, Tweedie
            }
            
            //   //this implementation does not follow calculation from van der Veen et al. 2021
            //   //but prevents Acov^-1 via woodbury matrix identity
            //   //see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
            //   //uses the identity (2D + A^-1) = A - 2A(I+2DA)^-1DA
            for (int i = 0; i < n; i++) {
              // Precompute Acov
              Acov.noalias() = A(i) * A(i).transpose();
              if (random(2)>0 && (num_lv_c+num_RR)>0)
                Acov += Ab_lvcov(i);
              
              // group Poisson, ZIP, Tweedie, NB, gamma, exponential
              if (is_group1) {
                // B = I - sign * 2 * D(j) * Acov
                B = Id - sign*2*D(j)*Acov;
                
                Eigen::PartialPivLU<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> lu(B);
                Binv = lu.inverse();
                Type logdetC = -lu.matrixLU().diagonal().array().log().sum();
                //the calculation generally prevents having to explicitly invert A*A^t, or having to invert A(i).
                Type vBinvv = ( 2*sign*newlam.col(j).transpose() * Acov * Binv * D(j) * Acov * newlam.col(j)
                      - 4*u.row(i) * Binv * D(j) * Acov * newlam.col(j)
                      + 2*sign*u.row(i) * Binv * D(j) * u.row(i).transpose()
                  ).value();
                
                //extra cQ contribution for  XB,e cross term in concurrent model
                if((random(2)<1) && (num_lv_c>0)){
                  vBinvv += (-4*x_lv.row(i)*b_lv2*Binv*D(j)*Acov*newlam.col(j)+2*sign*x_lv.row(i)*b_lv2*Binv*D(j)*u.row(i).transpose()+2*sign*u.row(i)*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()+2*sign*x_lv.row(i)*b_lv2*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                  // get rid of extra terms that will occur due to the use of eta + cQ in the likelihood
                  cQ(i,j) -= (sign*2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose()+sign*x_lv.row(i)*b_lv2*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                }
                
                //-logdetA + logdetB = logdetQ + logdetB = logdetC = det(QB^-1)
                // partialPivLU because B is asymmetric, and it ensures that it is invertible.
                logdetC = Binv.partialPivLu().matrixLU().diagonal().array().log().sum();
                
                cQ(i,j) += 0.5*(vBinvv+logdetC);
                // get rid of extra terms that will occur due to the use of eta + cQ in the likelihood
                cQ(i,j) -=  sign*(D(j)*Acov).trace() + sign*(u.row(i)*D(j)*u.row(i).transpose()).sum();
              }
              
              // Binomial, Gaussian, Ordinal
              if (is_group2) {
                cQ(i,j) += (D(j)*Acov*D(j)*Acov).trace() +2*(u.row(i)*D(j)*Acov*D(j)*u.row(i).transpose()).value() - 2*(u.row(i)*D(j)*Acov*newlam.col(j)).value();
                if((num_lv_c>0) && (random(2)<1)){
                  //extra terms for concurrent ordination
                  cQ(i,j) += (2*x_lv.row(i)*b_lv2*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose() -2*x_lv.row(i)*b_lv2*D(j)*Acov*newlam.col(j)+4*u.row(i)*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                }
              }
              
              // Eta-update
              eta(i,j) +=-(u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*Acov).trace();
              if ((num_lv_c>0) && (random(2)<1)) {
                eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
              }
              
            } //end for i
          } //end for j
          
          // //Poisson, NB, gamma, exponential, ZIP
          // if((family==0)||(family==1)||(family==4)||(family==6)||(family==8)||(family==11)){
          //   int sign = 1;
          //   //sign controls the family
          //   if((family>0) && (family != 6) && (family != 5)){
          //     //NB, gamma, exponential, ZIP
          //     sign = 1;
          //   }else if((family==0)||(family==5)||(family==6)){
          //     //Poisson, ZIP, Tweedie
          //     sign = -1;
          //   }
          //   
          //   matrix<Type> Binv(nlvr,nlvr);
          //   // Type logdetC;
          //   matrix <Type> Id(nlvr,nlvr);
          //   Id.setZero();Id.diagonal().fill(1.0);
          //   //this implementation does not follow calculation from van der Veen et al. 2021
          //   //but prevents Acov^-1 via woodbury matrix identity
          //   //see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
          //   //uses the identity (2D + A^-1) = A - 2A(I+2DA)^-1DA
          //   matrix<Type>B(nlvr,nlvr);
          //   for (int i=0; i<n; i++) {
          //     Acov = A(i)*A(i).transpose();
          //     if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
          //     for (int j=0; j<p;j++){
          //       B = Id-sign*2*D(j)*Acov;
          //       Eigen::PartialPivLU<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> lu(B);
          //       Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Binv = lu.inverse();
          //       Type logdetC = -lu.matrixLU().diagonal().array().log().sum();
          //       //the calculation generally prevents having to explicitly invert A*A^t, or having to invert A(i).
          //       Type vBinvv = (2*sign*newlam.col(j).transpose()*Acov*Binv*D(j)*Acov*newlam.col(j)-4*u.row(i)*Binv*D(j)*Acov*newlam.col(j) +2*sign*u.row(i)*Binv*D(j)*u.row(i).transpose()).value();
          //       
          //       //extra cQ contribution for  XB,e cross term in concurrent model
                // if((random(2)<1) && (num_lv_c>0)){
                //   vBinvv += (-4*x_lv.row(i)*b_lv2*Binv*D(j)*Acov*newlam.col(j)+2*sign*x_lv.row(i)*b_lv2*Binv*D(j)*u.row(i).transpose()+2*sign*u.row(i)*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()+2*sign*x_lv.row(i)*b_lv2*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                //   // get rid of extra terms that will occur due to the use of eta + cQ in the likelihood
                //   cQ(i,j) -= (sign*2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose()+sign*x_lv.row(i)*b_lv2*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                // }
          //       
          //       //-logdetA + logdetB = logdetQ + logdetB = logdetC = det(QB^-1)
          //       // partialPivLU because B is asymmetric, and it ensures that it is invertible.
          //       logdetC = Binv.partialPivLu().matrixLU().diagonal().array().log().sum();
          //       cQ(i,j) += 0.5*(vBinvv+logdetC);
          //       // get rid of extra terms that will occur due to the use of eta + cQ in the likelihood
          //       cQ(i,j) -=  sign*(D(j)*Acov).trace() + sign*(u.row(i)*D(j)*u.row(i).transpose()).sum();
          //     }
          //   }
          // }
          // Binomial, Gaussian, Ordinal
          // if((family==2)||(family==3)||(family==7)){
          //   for (int i=0; i<n; i++) {
          //     Acov = A(i)*A(i).transpose();
          //     if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
          //     for (int j=0; j<p;j++){
          //       cQ(i,j) += (D(j)*Acov*D(j)*Acov).trace() +2*(u.row(i)*D(j)*Acov*D(j)*u.row(i).transpose()).value() - 2*(u.row(i)*D(j)*Acov*newlam.col(j)).value();
          //       if((num_lv_c>0) && (random(2)<1)){
          //         //extra terms for concurrent ordination
          //         cQ(i,j) += (2*x_lv.row(i)*b_lv2*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose() -2*x_lv.row(i)*b_lv2*D(j)*Acov*newlam.col(j)+4*u.row(i)*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
          //       }
          //     }
          //   }
          // }
          
          // for (int i=0; i<n; i++) {
          //   Acov = A(i)*A(i).transpose();
          //   if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
          //   for (int j=0; j<p;j++){
          //     eta(i,j) += - (u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*A(i)*A(i).transpose()).trace();
          //     if((num_lv_c>0) && (random(2)<1)){
          //       eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
          //     }
          //   }
          // }
        }
      }
    }
    
    
    #include "family_nll.h"

  } // LA end
  return nll;
}
