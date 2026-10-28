// Helpers for exclude_obs(): bounding boxes of the likelihood supports and the vmax filter.
// Both work directly on the compressed-column storage of L (no triplet copy of L).
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
using namespace Rcpp;

//' Bounding box of the non-zero cells of each row of L (internal)
//' @param Lp,Li,Lx slots `p`, `i`, `x` of `L` (dgCMatrix, rows = observations, columns = cells)
//' @param nrow number of rows of `L`
//' @param x,y cell coordinates, one per column of `L`
//' @return `nrow x 4` matrix (xmin, xmax, ymin, ymax); `NA` for rows without non-zero cells
//' @noRd
// [[Rcpp::export]]
NumericMatrix row_support_bbox(const IntegerVector& Lp, const IntegerVector& Li, const NumericVector& Lx,
                               int nrow, const NumericVector& x, const NumericVector& y) {
  const int nc = Lp.size() - 1;
  if (x.size() != nc || y.size() != nc) stop("'x' and 'y' must have one value per column of L.");
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> x0(nrow, inf), x1(nrow, -inf), y0(nrow, inf), y1(nrow, -inf);
  for (int j = 0; j < nc; ++j) {
    const double xj = x[j], yj = y[j];
    for (int k = Lp[j]; k < Lp[j + 1]; ++k) {
      if (!(Lx[k] > 0)) continue;
      const int i = Li[k];
      if (xj < x0[i]) x0[i] = xj;
      if (xj > x1[i]) x1[i] = xj;
      if (yj < y0[i]) y0[i] = yj;
      if (yj > y1[i]) y1[i] = yj;
    }
  }
  NumericMatrix out(nrow, 4);
  for (int i = 0; i < nrow; ++i) {
    if (x0[i] == inf) { out(i, 0) = out(i, 1) = out(i, 2) = out(i, 3) = NA_REAL; continue; }
    out(i, 0) = x0[i]; out(i, 1) = x1[i]; out(i, 2) = y0[i]; out(i, 3) = y1[i];
  }
  return out;
}

namespace {
struct SpeedCtx {
  const NumericMatrix& bb; const NumericVector& tt; double buffer;
  // minimum speed needed between observations a and b (0-based rows)
  double spd(int a, int b) const {
    const double dx = std::max(0.0, std::max(bb(a, 0) - bb(b, 1), bb(b, 0) - bb(a, 1)));
    const double dy = std::max(0.0, std::max(bb(a, 2) - bb(b, 3), bb(b, 2) - bb(a, 3)));
    const double d = std::max(0.0, std::sqrt(dx * dx + dy * dy) - buffer);
    const double t = std::fabs(tt[b] - tt[a]);
    if (d == 0) return 0.0;
    return t > 0 ? d / t : R_PosInf;
  }
};
}

//' Maximum-speed filter (internal)
//' @param bb output of `row_support_bbox()`
//' @param tt cumulative time (`cumsum(dt)` with `dt = 0` at segment starts)
//' @param seg segment of each row
//' @param keep rows that take part (retained observations; rows without cells are skipped)
//' @param vmax,buffer maximum speed and distance tolerance
//' @return 1-based rows to exclude, in the order they are removed. Within each segment
//' (in order of first appearance) the fastest step between consecutive retained rows is
//' found (the first one if tied); if it exceeds `vmax`, whichever of its two rows leaves the
//' slower bypass between its neighbours is removed (ties: the earlier row, except on the
//' last step of the segment, where the last row is removed), until no step exceeds `vmax`.
//' @noRd
// [[Rcpp::export]]
IntegerVector speed_filter_cpp(const NumericMatrix& bb, const NumericVector& tt, const IntegerVector& seg,
                               const LogicalVector& keep, double vmax, double buffer) {
  const int N = seg.size();
  if (tt.size() != N || keep.size() != N || bb.nrow() != N) stop("Inconsistent lengths in speed_filter_cpp().");
  SpeedCtx S{bb, tt, buffer};
  // rows of each segment, segments in order of first appearance
  std::vector<int> seg_order; std::vector<std::vector<int>> rows;
  {
    std::vector<int> lab; std::vector<int> pos;  // small linear map (few segments)
    for (int i = 0; i < N; ++i) {
      if (keep[i] != TRUE || ISNAN(bb(i, 0))) continue;   // rows without cells never take part
      int s = seg[i], g = -1;
      for (size_t q = 0; q < lab.size(); ++q) if (lab[q] == s) { g = (int)q; break; }
      if (g < 0) { g = (int)lab.size(); lab.push_back(s); rows.emplace_back(); }
      rows[g].push_back(i);
    }
  }
  std::vector<int> flagged;
  for (auto& r : rows) {
    const int n = r.size();
    if (n < 2) continue;
    // doubly linked list over positions 0..n-1 of the retained rows; sp[q] = speed of the step
    // from position q to its successor
    std::vector<int> nxt(n), prv(n); std::vector<double> sp(n, 0.0);
    for (int q = 0; q < n; ++q) { nxt[q] = q + 1; prv[q] = q - 1; }
    for (int q = 0; q + 1 < n; ++q) sp[q] = S.spd(r[q], r[q + 1]);
    int head = 0, m = n;
    while (m >= 2) {
      // first maximum along the list (as which.max)
      int jq = -1; double best = -1.0;
      for (int q = head; nxt[q] < n; q = nxt[q]) if (sp[q] > best) { best = sp[q]; jq = q; }
      if (!(best > vmax)) break;
      const int aq = jq, bq = nxt[jq];                  // the step aq -> bq
      const int pa = prv[aq], nb = nxt[bq];               // neighbours (-1 / n if none)
      const double by_a = (pa >= 0) ? S.spd(r[pa], r[bq]) : 0.0;     // remove aq
      const double by_b = (nb < n) ? S.spd(r[aq], r[nb]) : 0.0;      // remove bq
      int dq;
      if (by_b < by_a) dq = bq; else if (by_a < by_b) dq = aq;
      else dq = (nb >= n) ? bq : aq;                      // tie: the earlier row, except on the last step (the last row)
      // unlink dq and update the one step that changes
      const int p = prv[dq], q2 = nxt[dq];
      if (p >= 0) { nxt[p] = q2; sp[p] = (q2 < n) ? S.spd(r[p], r[q2]) : 0.0; } else head = q2;
      if (q2 < n) prv[q2] = p;
      flagged.push_back(r[dq] + 1);
      --m;
    }
  }
  return wrap(flagged);
}

//' Remove rows from a dgCMatrix (and any stored zeros) without a matrix product (internal)
//' @return list(p, i, x) for the new matrix with the same dimensions
//' @noRd
// [[Rcpp::export]]
List drop_rows_csc(const IntegerVector& Lp, const IntegerVector& Li, const NumericVector& Lx,
                   const LogicalVector& drop) {
  const int nc = Lp.size() - 1;
  const R_xlen_t nnz = Lx.size();
  R_xlen_t keepn = 0;
  for (R_xlen_t k = 0; k < nnz; ++k) if (drop[Li[k]] != TRUE && Lx[k] != 0) ++keepn;
  IntegerVector p(nc + 1), i(keepn); NumericVector x(keepn);
  R_xlen_t o = 0;
  for (int j = 0; j < nc; ++j) {
    for (int k = Lp[j]; k < Lp[j + 1]; ++k) {
      if (drop[Li[k]] == TRUE || Lx[k] == 0) continue;
      i[o] = Li[k]; x[o] = Lx[k]; ++o;
    }
    p[j + 1] = (int)o;
  }
  return List::create(_["p"] = p, _["i"] = i, _["x"] = x);
}
