#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double NRMSE(arma::mat output, arma::mat target) {
  arma::mat diff = output - target;
  arma::mat sqdiff = pow(diff, 2);
  arma::rowvec msqdiff = arma::mean(sqdiff, 0);
  arma::rowvec colvars = 0.5 * arma::var(output, 0, 0) + 0.5 * arma::var(target, 0, 0);
  arma::rowvec rmsqdiff = pow(msqdiff / colvars, 0.5);
  double avgNRMSE = arma::mean(rmsqdiff);
  return avgNRMSE;
}

// [[Rcpp::export]]
arma::cube runRNN(arma::mat input, arma::cube W, float RSInitial, float bscale, float iscale) {
  float N = W.n_rows;
  float nfits = W.n_slices;
  float time_steps = input.n_rows;
  float dimen = input.n_cols;

  arma::cube RS = RSInitial * arma::cube(N, nfits, 1, arma::fill::ones);
  arma::cube ResStates = arma::join_slices(RS, arma::cube(N, nfits, time_steps));

  for(int t = 0; t < time_steps; t++) {
    for(int i = 0; i < nfits; i++) {
      ResStates.slice(t + 1).col(i) = tanh(bscale * W.slice(i).col(0) + iscale * W.slice(i).cols(1, dimen) * trans(input.row(t)) + W.slice(i).cols(dimen + 1, dimen + N) * ResStates.slice(t).col(i));
    }
  }

  return ResStates;
}

// [[Rcpp::export]]
arma::cube conceptorCalc(arma::cube RS, float aperture) {
  float N = RS.n_rows;
  float nfits = RS.n_cols;
  float time_steps = RS.n_slices;

  arma::cube C = arma::cube(N, N, nfits);

  for(int i = 0; i < nfits; i++) {
    arma::mat RStemp = RS.col(i);
    arma::mat R = RStemp * trans(RStemp) / time_steps;
    arma::mat U, V, O;
    arma::vec D;
    arma::svd(U, D, V, R);
    arma::mat S = arma::diagmat(D) * inv(arma::diagmat(D) + pow(aperture, -2) * O.eye(N, N));
    C.slice(i) = U * S * trans(U);
  }
  return C;
}

// [[Rcpp::export]]
arma::cube runCRNN(arma::mat input, arma::cube W, arma::cube C, float RSInitial, int washoutL, float bscale, float iscale) {
  float N = W.n_rows;
  float nfits = W.n_slices;
  float time_steps = input.n_rows;
  float dimen = input.n_cols;

  arma::cube RS = RSInitial * arma::cube(N, nfits, 1, arma::fill::ones);
  arma::cube ResStates = arma::join_slices(RS, arma::cube(N, nfits, time_steps));
  arma::cube CResStates = ResStates;

  for(int t = 0; t < washoutL; t++) {
    for(int i = 0; i < nfits; i++) {
      ResStates.slice(t + 1).col(i) = tanh(bscale * W.slice(i).col(0) + iscale * W.slice(i).cols(1, dimen) * trans(input.row(t)) + W.slice(i).cols(dimen + 1, dimen + N) * ResStates.slice(t).col(i));
      CResStates.slice(t + 1).col(i) = ResStates.slice(t + 1).col(i);
    }
  }

  for(int t = washoutL; t < time_steps; t++) {
    for(int i = 0; i < nfits; i++) {
      ResStates.slice(t + 1).col(i) = tanh(bscale * W.slice(i).col(0) + iscale * W.slice(i).cols(1, dimen) * trans(input.row(t)) + W.slice(i).cols(dimen + 1, dimen + N) * ResStates.slice(t).col(i));
      CResStates.slice(t + 1).col(i) = C.slice(i) * ResStates.slice(t + 1).col(i);
    }
  }
  arma::cube CRNN = arma::join_slices(ResStates, CResStates);

  return CRNN;
}

// [[Rcpp::export]]
arma::cube WoutCalc(arma::mat input, arma::cube ResStates, float regular) {
  float N = ResStates.n_rows;
  float nfits = ResStates.n_cols;
  float dimen = input.n_cols;
  arma::cube Wout = arma::cube(N, dimen, nfits);

  for(int j = 0; j < nfits; j++) {
    arma::mat RS = ResStates.col(j);
    arma::mat WS = inv(RS * trans(RS) + regular * arma::mat(N, N, arma::fill::eye)) * RS * input;
    Wout.slice(j) = WS;
  }

  return Wout;
}

// [[Rcpp::export]]
arma::cube outputCalc(arma::cube Wout, arma::cube CResStates) {
  float nfits = Wout.n_slices;
  float ts = CResStates.n_slices;
  float dimen = Wout.n_cols;
  arma::cube output = arma::cube(ts, dimen, nfits);

  for(int j = 0; j < nfits; j++) {
    arma::mat WS = Wout.slice(j);
    arma::mat CRS = CResStates.col(j);
    arma::mat os = trans(WS) * CRS;
    output.slice(j) = trans(os);
  }

  return output;
}

// [[Rcpp::export]]
arma::mat RNNParamFit(arma::mat input, arma::cube W, float RSInitial, int washoutL, int trainL, arma::vec bscales, arma::vec iscales) {
  float nfits = W.n_slices;
  float nb = bscales.n_elem;
  float ni = iscales.n_elem;
  float regular = 1e-4;

  arma::mat trainErrors = arma::mat(nb, ni);
  for(int b = 0; b < nb; b++) {
    for(int i = 0; i < ni; i++) {
      arma::cube rnn0 = runRNN(input, W, 0, bscales(b), iscales(i));
      arma::cube Wout = WoutCalc(input.rows(washoutL, washoutL + trainL - 1), rnn0.slices(washoutL + 1, washoutL + trainL), regular);
      arma::cube output = outputCalc(Wout, rnn0.slices(washoutL + 1, washoutL + trainL));

      arma::vec errors = arma::vec(nfits);
      for(int j = 0; j < nfits; j++) {
        errors(j) = NRMSE(output.slice(j), input.rows(washoutL, washoutL + trainL - 1));
      }
      trainErrors(b, i) = arma::mean(errors);
    }
  }
  return trainErrors;
}

// [[Rcpp::export]]
arma::mat angleCalc(arma::cube RS, arma::cube CRS) {
  float time_steps = RS.n_slices;
  float nfits = RS.n_cols;

  arma::mat angles = arma::mat(time_steps, nfits);

  for(int i = 0; i < nfits; i++) {
    arma::mat CRSt = CRS.col(i);
    arma::mat RSt = RS.col(i);
    angles.col(i) = diagvec(trans(CRSt) * RSt)/ sqrt(diagvec(trans(CRSt) * CRSt)) / sqrt(diagvec(trans(RSt) * RSt));
  }
  return angles;
}

// [[Rcpp::export]]
arma::vec KSstatCalc(arma::vec angles) {
  float timepts = angles.n_elem;
  arma::vec timevector = arma::linspace(1, timepts, timepts);
  arma::uvec asort = arma::sort_index(angles) + 1;

  arma::mat KSmat = arma::mat(timepts, timepts);
  for(int t1 = 0; t1 < timepts; t1++) {
    for(int t2 = 0; t2 < timepts; t2++) {
      if(asort(t1) <= timevector(t2)) {
        KSmat(t1, t2) = 1 / timevector(t2);
      }else{
        KSmat(t1, t2) = -1 / (timepts - timevector(t2));
      }
    }
  }

  KSmat = abs(arma::cumsum(KSmat));
  arma::vec cmaxes = arma::vec(timepts);
  for(int t1 = 0; t1 < timepts; t1++) {
    cmaxes(t1) = arma::max(KSmat.col(t1));
  }


  arma::vec KSseries = arma::vec(timepts);
  for(int t1 = 0; t1 < timepts; t1++){
    double coef1 = timevector(t1) * (timepts - timevector(t1)) / pow(timepts, 1.5);
    arma::vec qscale = {pow(timevector(t1) / timepts, 0.5) * pow(1 - timevector(t1) / timepts, 0.5), 0.01};
    double coef2 = arma::max(qscale);
    KSseries(t1) = coef1 * cmaxes(t1) / coef2;
  }

  return KSseries;
}

// [[Rcpp::export]]
arma::vec CRNNBootstrap(arma::cube bootinput, arma::cube W, arma::cube C, float RSInitial, int washoutL, int trainL, float bscale, float iscale) {
  float nboots = bootinput.n_slices;
  float L = bootinput.n_rows;
  arma::vec mbbKS = arma::vec(nboots);

  for(int i = 0; i < nboots; i++) {
    arma::cube crnn = runCRNN(bootinput.slice(i), W, C, RSInitial, washoutL, bscale, iscale);
    arma::mat Angles = angleCalc(crnn.slices(1, L), crnn.slices(L + 2, 2 * L + 1));
    arma::vec aAngles = arma::mean(Angles, 1);
    arma::vec anglesKS = aAngles.rows(washoutL + trainL, L-1);
    arma::vec KS = KSstatCalc(anglesKS);
    float maxKS = max(KS);
    mbbKS(i) = maxKS;
  }

  return mbbKS;
}
