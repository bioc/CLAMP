// [[Rcpp::depends(RcppArmadillo)]]
//#define ARMA_DONT_USE_CXX11
//#define ARMA_USE_CXX11_RNG
#include <RcppArmadillo.h>




// [[Rcpp::export]]
void updateZcpp(SEXP ZSEXP,
                   SEXP YSEXP,
                   SEXP BSEXP,
                   SEXP B2SEXP,
                   const arma::vec& bk2_all,
                   SEXP YBtSEXP,
                   SEXP QSEXP,
                   SEXP ZmultiplierSEXP,
                   double L1,
                   double multiplier,
                   int inner_iter,
                   int iter,
                   bool Y_is_fbm = false,
                   int n = 0,
                   int p = 0) {

  // writable views for Z and Q (no copy, modifies R object directly)
  Rcpp::NumericMatrix Zr(ZSEXP);
  arma::mat Z(Zr.begin(), Zr.nrow(), Zr.ncol(), false);

  Rcpp::NumericMatrix Qr(QSEXP);
  arma::mat Q(Qr.begin(), Qr.nrow(), Qr.ncol(), false);

  // read-only matrices (copied safely)
  arma::mat B           = Rcpp::as<arma::mat>(BSEXP);
  arma::mat B2          = Rcpp::as<arma::mat>(B2SEXP);
  arma::mat YBt         = Rcpp::as<arma::mat>(YBtSEXP);
  arma::mat Zmultiplier = Rcpp::as<arma::mat>(ZmultiplierSEXP);

  // handle Y: normal or FBM
  arma::mat Y;
  if (Y_is_fbm) {
    void* Y_ptr = R_ExternalPtrAddr(YSEXP);
    if (Y_ptr == nullptr)
      Rcpp::stop("Invalid FBM pointer");
    Y = arma::mat(static_cast<double*>(Y_ptr), n, p, false, true);
  } else {
    Y = Rcpp::as<arma::mat>(YSEXP);
  }

  // main update loop
  int k = Z.n_cols;
  for (int inner = 0; inner < inner_iter; ++inner) {
//     arma::arma_rng::set_seed(inner);
// //    std::mt19937 rng(inner);
//     Rcpp::Rcout << "RNG check (after seed " << inner << "): "
//                 << arma::randu() << "\n";
//     arma::uvec order = arma::shuffle(arma::linspace<arma::uvec>(0, k - 1, k));

    arma::uvec order = arma::linspace<arma::uvec>(0, k - 1, k);
    std::mt19937 rng(static_cast<uint32_t>(inner*iter));
    std::shuffle(order.begin(), order.end(), rng);


    // print update order
    // Rcpp::Rcout << "inner=" << inner << " order: ";
    // for (unsigned int jj : order) Rcpp::Rcout << jj << " ";
    // Rcpp::Rcout << "\n";

    for (unsigned int j : order) {
      arma::vec gene_var = L1 / (multiplier * Zmultiplier.col(j) + 1.0);
      double bk2 = bk2_all[j];
      arma::vec denom = bk2 + gene_var;
      arma::vec num = YBt.col(j) - Q.col(j) + Z.col(j) * bk2;
      arma::vec newZk = num / denom;
      arma::vec delta = newZk - Z.col(j);
      Z.col(j) = newZk;
      Q += delta * B2.row(j);
    }
  }
}


// [[Rcpp::export]]
void updateZcppOld(arma::mat& Z,
                const arma::mat& Y,
                const arma::mat& B,
                const arma::mat& B2,
                const arma::vec& bk2_all,
                const arma::mat& YBt,
                arma::mat& Q,
                const arma::mat& Zmultiplier,
                double L1,
                double multiplier,
                int inner_iter,
                int iter) {

  int n = Z.n_rows, k = Z.n_cols;
  for (int inner = 0; inner < inner_iter; ++inner) {
    arma::arma_rng::set_seed(inner * iter);
    arma::uvec order = arma::shuffle(arma::linspace<arma::uvec>(0, k - 1, k));
    for (unsigned int j : order) {
      arma::vec gene_var = L1 / (multiplier * Zmultiplier.col(j) + 1.0);
      double bk2 = bk2_all[j];
      arma::vec denom = bk2 + gene_var;
      arma::vec num = YBt.col(j) - Q.col(j) + Z.col(j) * bk2;
      arma::vec newZk = num / denom;
      arma::vec delta = newZk - Z.col(j);
      Z.col(j) = newZk;
      Q += delta * B2.row(j);
    }
  }
}
