// [[Rcpp::depends(RcppArmadillo)]]
#include "ctmc_common.h"
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// ===========================================================================
// Joint (multi-individual) CTMC likelihood with support tracking.
//
// The observations of all individuals are stacked in time order *within*
// individual.  `seg_start` gives the (0-based) first observation of every
// independent chain ("segment": an individual, or a burst of an individual
// split at a long gap).  All segments share one rate matrix Q, so this is a
// single joint model with common parameters:
//
//      -2 log L(theta) = sum_k  -2 log L_k(theta)
//
// Each chain starts from `delta` at its first observation.
//
// Two things make this much faster than calling the single-track forward
// algorithm repeatedly:
//
//  (1) Support tracking.  After the observation update the filtered
//      distribution phi is zero outside the support of L_i (the error
//      ellipse).  Each uniformisation mat-vec can only spread mass to the
//      neighbours of currently non-zero cells, so we keep an explicit
//      active set and only touch those cells.  This is EXACT (no
//      approximation): cells outside the active set are structurally 0.
//      An optional `prune` threshold drops filtered probabilities below
//      prune * max(phi) (default 0 = exact).
//
//  (2) Segments are independent given theta, so they are evaluated in
//      parallel with OpenMP when available (`nthreads`).
//
// The uniformisation matrix P = I + Q/alpha_0 is built once per call, and
// its transpose is stored in CSC form so that row i of P (the cells that
// cell i can move to) is contiguous in memory.
// ===========================================================================

namespace {

struct Workspace {
  std::vector<double> vc, vn, res;          // current, next, accumulated
  std::vector<unsigned int> mark_n, mark_r; // stamps for active sets
  std::vector<arma::uword> act, nact, ract; // active index lists
  unsigned int stamp_n = 0, stamp_r = 0;
  explicit Workspace(arma::uword ns)
      : vc(ns, 0.0), vn(ns, 0.0), res(ns, 0.0), mark_n(ns, 0), mark_r(ns, 0) {
    act.reserve(1024); nact.reserve(1024); ract.reserve(1024);
  }
};

// poisson_window(): see ctmc_common.h (shared with v_exp_Q_t.cpp so results match).

// phi (given by values ws.vc on index list ws.act) -> phi * exp(Q t)
// stored in ws.res on index list ws.ract.  On exit ws.vc/ws.act are cleared.
void propagate(Workspace& ws, const arma::uword* Pcp, const arma::uword* Pri,
               const double* Pv, double alpha_0, double t, double tol) {
  ws.ract.clear();
  ++ws.stamp_r;
  if (t <= 0.0 || alpha_0 <= 0.0) {
    for (arma::uword i : ws.act) {
      ws.res[i] = ws.vc[i];
      ws.mark_r[i] = ws.stamp_r;
      ws.ract.push_back(i);
      ws.vc[i] = 0.0;
    }
    ws.act.clear();
    return;
  }
  double alpha = alpha_0 * t;
  int m_lower, m;
  poisson_window(alpha, tol, m_lower, m);
  double la = std::log(alpha);

  for (int n = 0; n <= m; ++n) {
    if (n >= m_lower) {
      double pn = std::exp(n * la - alpha - std::lgamma(n + 1));
      for (arma::uword i : ws.act) {
        if (ws.mark_r[i] != ws.stamp_r) {
          ws.mark_r[i] = ws.stamp_r;
          ws.res[i] = 0.0;
          ws.ract.push_back(i);
        }
        ws.res[i] += pn * ws.vc[i];
      }
    }
    if (n < m) {
      // v_next = v_curr * P restricted to the active set.
      // Pt is CSC of P^T, so column i of Pt = row i of P.
      ++ws.stamp_n;
      ws.nact.clear();
      for (arma::uword i : ws.act) {
        double vi = ws.vc[i];
        for (arma::uword k = Pcp[i]; k < Pcp[i + 1]; ++k) {
          arma::uword j = Pri[k];
          if (ws.mark_n[j] != ws.stamp_n) {
            ws.mark_n[j] = ws.stamp_n;
            ws.vn[j] = 0.0;
            ws.nact.push_back(j);
          }
          ws.vn[j] += vi * Pv[k];
        }
      }
      for (arma::uword i : ws.act) ws.vc[i] = 0.0;
      std::swap(ws.vc, ws.vn);
      std::swap(ws.act, ws.nact);
    }
  }
  for (arma::uword i : ws.act) ws.vc[i] = 0.0;
  ws.act.clear();
}

// One independent chain over observations [i0, i1).
double segment_n2ll(arma::uword i0, arma::uword i1, const arma::uword* Lcp,
                    const arma::uword* Lri, const double* Lv,
                    const arma::vec& dt, const arma::rowvec& delta, double p,
                    arma::uword ns, const arma::uword* Pcp,
                    const arma::uword* Pri, const double* Pv, double alpha_0,
                    double tol, double prune, bool skip_zero,
                    std::vector<int>& zobs, Workspace& ws) {
  double ll = 0.0;
  bool full = (p > 0.0);  // outlier mixture => support is everything we reach

  // ---- first observation: phi = delta * L_i (or delta if L_i empty) -----
  {
    ws.act.clear();
    double u = 0.0;
    arma::uword k0 = Lcp[i0], k1 = Lcp[i0 + 1];
    if (k1 > k0 && !full) {
      for (arma::uword k = k0; k < k1; ++k) {
        arma::uword j = Lri[k];
        double val = delta[j] * Lv[k];
        if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; }
      }
    } else {
      // empty L row, or outlier mixture: start from all of delta
      std::vector<double> Ld;
      if (k1 > k0) { Ld.assign(ns, 0.0); for (arma::uword k = k0; k < k1; ++k) Ld[Lri[k]] = Lv[k]; }
      for (arma::uword j = 0; j < ns; ++j) {
        double val = delta[j];
        if (k1 > k0) val = val * ((1.0 - p) * Ld[j]) + (p / ns) * val;
        if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; }
      }
    }
    if (!(u > 0.0)) return R_PosInf;
    ll += std::log(u);
    for (arma::uword j : ws.act) ws.vc[j] /= u;
  }

  // ---- forward recursion ------------------------------------------------
  for (arma::uword i = i0 + 1; i < i1; ++i) {
    propagate(ws, Pcp, Pri, Pv, alpha_0, dt(i), tol);   // -> ws.res on ws.ract
    arma::uword k0 = Lcp[i], k1 = Lcp[i + 1];
    double u = 0.0;
    ws.act.clear();
    if (k1 == k0) {
      // no information at this time: v = phi exp(Qt)
      for (arma::uword j : ws.ract) { double val = ws.res[j]; if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; } }
    } else if (!full) {
      // v = phi exp(Qt) % L_i, on supp(L_i) intersect supp(res)
      for (arma::uword k = k0; k < k1; ++k) {
        arma::uword j = Lri[k];
        if (ws.mark_r[j] != ws.stamp_r) continue;
        double val = ws.res[j] * Lv[k];
        if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; }
      }
    } else {
      // outlier mixture: v = res % ((1-p) L + p/ns)
      for (arma::uword j : ws.ract) ws.vn[j] = (p / ns) * ws.res[j];
      for (arma::uword k = k0; k < k1; ++k) {
        arma::uword j = Lri[k];
        if (ws.mark_r[j] == ws.stamp_r) ws.vn[j] += (1.0 - p) * Lv[k] * ws.res[j];
      }
      for (arma::uword j : ws.ract) { double val = ws.vn[j]; ws.vn[j] = 0.0; if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; } }
    }
    if (!(u > 0.0)) {
      // The observation is impossible under the current Q (its error
      // support does not overlap the propagated distribution).
      for (arma::uword j : ws.act) ws.vc[j] = 0.0;
      ws.act.clear();
      zobs.push_back((int)i + 1);   // 1-based observation index
      if (!skip_zero) {
        for (arma::uword j : ws.ract) ws.res[j] = 0.0;
        return R_PosInf;
      }
      // treat as missing: carry the propagated distribution forward
      u = 0.0;
      for (arma::uword j : ws.ract) { double val = ws.res[j]; if (val > 0.0) { ws.vc[j] = val; ws.act.push_back(j); u += val; } }
      for (arma::uword j : ws.ract) ws.res[j] = 0.0;
      if (!(u > 0.0)) return R_PosInf;
      for (arma::uword j : ws.act) ws.vc[j] /= u;   // u ~ 1, not added to ll
      continue;
    }
    for (arma::uword j : ws.ract) ws.res[j] = 0.0;
    ll += std::log(u);
    double mx = 0.0;
    for (arma::uword j : ws.act) { ws.vc[j] /= u; if (ws.vc[j] > mx) mx = ws.vc[j]; }
    if (prune > 0.0) {
      double thr = prune * mx, kept = 0.0;
      ws.nact.clear();
      for (arma::uword j : ws.act) {
        if (ws.vc[j] >= thr) { ws.nact.push_back(j); kept += ws.vc[j]; }
        else ws.vc[j] = 0.0;
      }
      std::swap(ws.act, ws.nact);
      for (arma::uword j : ws.act) ws.vc[j] /= kept;
    }
  }
  for (arma::uword j : ws.act) ws.vc[j] = 0.0;
  ws.act.clear();
  return -2.0 * ll;
}

}  // namespace

//' Joint multi-segment CTMC -2 log-likelihood (internal)
//' @return list(n2ll = per-segment -2 log L, zero_obs = 1-based indices of
//'   observations with zero predictive probability)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ctmc_n2ll_joint_arma(
    const arma::sp_mat& L,          // ns x N, one column per observation
    const arma::vec& dt,
    const arma::uvec& seg_start,    // 0-based, increasing, seg_start[0] == 0
    int ns,
    const arma::umat& from_to,
    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
    double p,
    const arma::rowvec& delta,
    const arma::vec& hij,
    double eq_prec,
    int link_r, double a_r, double l_r, double u_r,
    int link_m, double a_m,
    int form, double k, bool norm, double clip,
    int nthreads = 1, double prune = 0.0, bool skip_zero = false)
{
  arma::uword N = dt.n_elem;
  if (L.n_cols != N) Rcpp::stop("ncol(L) must equal length(dt)");
  if ((int)L.n_rows != ns) Rcpp::stop("nrow(L) must equal ns");
  check_seg_start(seg_start, N);

  arma::sp_mat Q = build_Q(from_to, Xb_q_r, Xb_q_m, hij, ns, link_r, a_r, l_r, u_r,
                           link_m, a_m, form, k, norm, clip);
  arma::sp_mat P;
  double alpha_0 = uniformization_P(Q, P);
  arma::sp_mat Pt = P.t();
  Pt.sync();
  arma::sp_mat Lc(L);
  Lc.sync();

  const arma::uword* Pcp = Pt.col_ptrs;
  const arma::uword* Pri = Pt.row_indices;
  const double* Pv = Pt.values;
  const arma::uword* Lcp = Lc.col_ptrs;
  const arma::uword* Lri = Lc.row_indices;
  const double* Lv = Lc.values;

  arma::uword S = seg_start.n_elem;
  arma::vec out(S, arma::fill::zeros);
  std::vector< std::vector<int> > zo(S);

  // Thread count: never more than requested, than segments, or than OpenMP allows.
  int nt = 1;
#ifdef _OPENMP
  nt = std::max(1, nthreads);
  nt = std::min(nt, omp_get_thread_limit());
  nt = std::min(nt, (int)std::max<arma::uword>(1, S));
#else
  (void)nthreads;
#endif

  // All allocation happens here, on the main thread, BEFORE the parallel
  // region: a std::bad_alloc is then caught by Rcpp and becomes an R error.
  // Inside the region nothing allocates except push_back() on the per-thread
  // index vectors / zero-obs lists, which is guarded below so that an
  // exception can never escape the OpenMP region (that would terminate R).
  // No R API is called inside the region.
  std::vector<Workspace> wss;
  wss.reserve(nt);
  for (int t = 0; t < nt; ++t) wss.emplace_back(ns);
  for (arma::uword s = 0; s < S; ++s) zo[s].reserve(4);
  int failed = 0;

  // Work distribution: a plain parallel region in which threads take the next
  // segment from a shared atomic counter (dynamic load balancing without
  // `schedule(dynamic)`). We deliberately avoid `omp for schedule(dynamic)`:
  // clang >= 19 compiles it to a call to __kmpc_dispatch_deinit, which older
  // libomp runtimes (e.g. some Homebrew / CRAN-shipped macOS builds) do not
  // export, so walk.so then fails to load ("symbol not found in flat
  // namespace '___kmpc_dispatch_deinit'"). This form needs only
  // __kmpc_fork_call. Each segment writes only its own slot out(s), and the
  // slots are summed later in a fixed order, so results do not depend on the
  // number of threads or on which thread ran which segment.
  long long next_seg = 0;
  const long long S_ll = (long long)S;
#ifdef _OPENMP
  #pragma omp parallel num_threads(nt)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    for (;;) {
      long long sl;
#ifdef _OPENMP
      #pragma omp atomic capture
#endif
      sl = next_seg++;
      if (sl >= S_ll) break;
      arma::uword s = (arma::uword)sl;
      arma::uword i0 = seg_start(s), i1 = (s + 1 < S) ? seg_start(s + 1) : N;
      try {
        out(s) = segment_n2ll(i0, i1, Lcp, Lri, Lv, dt, delta, p, ns, Pcp, Pri, Pv,
                              alpha_0, eq_prec, prune, skip_zero, zo[s], wss[tid]);
      } catch (...) {
#ifdef _OPENMP
        #pragma omp atomic write
#endif
        failed = 1;
        out(s) = R_NaN;
      }
    }
  }
  if (failed) Rcpp::stop("ctmc_n2ll_joint_arma: out of memory (or other C++ error) in a worker thread.");

  std::vector<int> zall;
  for (arma::uword s = 0; s < S; ++s) zall.insert(zall.end(), zo[s].begin(), zo[s].end());
  return Rcpp::List::create(Rcpp::Named("n2ll") = out,
                            Rcpp::Named("zero_obs") = Rcpp::wrap(zall),
                            Rcpp::Named("alpha_0") = alpha_0);
}

//' Report whether walk was compiled with OpenMP
//' @keywords internal
// [[Rcpp::export]]
bool walk_has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
