#include "milp/heuristics/lmip.hpp"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>

#include "milp/core/interrupt.hpp"
#include "milp/heuristics/fj.hpp"

namespace milp {

namespace {

constexpr double kTol = 1e-7;  // row violation tolerance (gap = a^T x - b); try_incumbent checks 1e-6
constexpr double kHuge = 1e12;  // moves to |values| beyond this are not considered

// BMS sample sizes and tabu tenure: the Local-MIP 2.0 defaults.
// MILP_LMIP_BMS="unsat rows,unsat ops,sat rows,sat ops,easy columns,random ops" overrides them (research).
struct Bms {
  int unsat_rows = 10, unsat_ops = 2250, sat_rows = 1, sat_ops = 80, easy = 5, random_ops = 250;
  Bms() {
    if (const char* v = std::getenv("MILP_LMIP_BMS")) {
      std::sscanf(v, "%d,%d,%d,%d,%d,%d", &unsat_rows, &unsat_ops, &sat_rows, &sat_ops, &easy, &random_ops);
    }
  }
};
const Bms& bms() {
  static const Bms b;
  return b;
}
constexpr int kTabuBase = 4, kTabuVar = 7;
// Breakthrough moves: at most this many objective columns per neighbourhood (a window at a random
// offset; Local-MIP takes all of them, which is O(#objective columns) per step).
constexpr int kBreakCols = 1000;

class LocalMip {
 public:
  LocalMip(const MipProblem& p, const LocalMipLimits& lim)
      : p_(p), lim_(lim), rng_(lim.seed * 2654435761u + 97u) {
    n_ = p.num_cols();
    // Normalised rows a^T x <= b: the upper side as is, the lower side negated.
    const Index m0 = p.num_rows();
    std::vector<Index> up(m0, -1), dn(m0, -1);
    for (Index i = 0; i < m0; ++i) {
      if (is_finite(p.row_ub[i])) {
        up[i] = m_++;
        b_.push_back(p.row_ub[i]);
        eq_half_.push_back(p.row_lb[i] == p.row_ub[i]);
      }
      if (is_finite(p.row_lb[i])) {
        dn[i] = m_++;
        b_.push_back(-p.row_lb[i]);
        eq_half_.push_back(p.row_lb[i] == p.row_ub[i]);
      }
    }
    cbeg_.assign(n_ + 1, 0);
    std::vector<Index> rcount(m_, 0);
    for (Index j = 0; j < n_; ++j) {
      for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
        const double a = p.A.value[k];
        if (a == 0.0) continue;
        const Index i = p.A.index[k];
        if (up[i] >= 0) {
          crow_.push_back(up[i]);
          cval_.push_back(a);
          ++rcount[up[i]];
        }
        if (dn[i] >= 0) {
          crow_.push_back(dn[i]);
          cval_.push_back(-a);
          ++rcount[dn[i]];
        }
      }
      cbeg_[j + 1] = static_cast<Index>(crow_.size());
    }
    rbeg_.assign(m_ + 1, 0);
    for (Index i = 0; i < m_; ++i) rbeg_[i + 1] = rbeg_[i] + rcount[i];
    rcol_.resize(crow_.size());
    rval_.resize(crow_.size());
    std::vector<Index> fill(rbeg_.begin(), rbeg_.end() - 1);
    for (Index j = 0; j < n_; ++j) {
      for (Index k = cbeg_[j]; k < cbeg_[j + 1]; ++k) {
        rcol_[fill[crow_[k]]] = j;
        rval_[fill[crow_[k]]++] = cval_[k];
      }
    }
    objpos_.assign(n_, -1);
    for (Index j = 0; j < n_; ++j) {
      if (p.cost[j] != 0.0 && p.col_lb[j] < p.col_ub[j]) {
        objpos_[j] = static_cast<Index>(objvars_.size());
        objvars_.push_back(j);
      }
      if (p.col_lb[j] < p.col_ub[j]) nonfixed_.push_back(j);
    }
    obj_integral_ = true;
    for (Index j = 0; j < n_ && obj_integral_; ++j) {
      if (p.cost[j] != 0.0) obj_integral_ = p.is_int[j] && p.cost[j] == std::round(p.cost[j]);
    }
    stamp_.assign(n_, 0);
  }

  bool run(const std::function<void(const std::vector<double>&, double)>& on_solution, LocalMipStats* stats) {
    t0_ = std::chrono::steady_clock::now();
    x_.resize(n_);
    if (lim_.start && static_cast<Index>(lim_.start->size()) == n_) {
      for (Index j = 0; j < n_; ++j) x_[j] = clamp_var(j, (*lim_.start)[j]);
    } else {
      for (Index j = 0; j < n_; ++j) x_[j] = clamp_var(j, 0.0);
    }
    if (is_finite(lim_.cutoff)) {
      // A known solution value: aim below it from the start.
      found_ = true;
      best_obj_ = lim_.cutoff;
      obj_bound_ = lim_.cutoff - improvement(lim_.cutoff);
    }
    w_.assign(m_, 1);
    allow_inc_.assign(n_, 0);
    allow_dec_.assign(n_, 0);
    last_move_.assign(n_, 0);
    lift_delta_.assign(objvars_.size(), 0.0);
    recompute();
    long long progress_effort = 0;
    min_unsat_ = static_cast<long long>(unsat_.size()) + 1;
    long long next_check = 0;
    while (effort_ < lim_.effort) {
      if (effort_ >= next_check) {
        next_check = effort_ + 100000;
        if ((lim_.stop && lim_.stop->load(std::memory_order_relaxed)) || interrupted()) break;
        if (lim_.checkpoint && !lim_.checkpoint(effort_)) break;
        if (lim_.time_limit < kInf &&
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count() > lim_.time_limit) {
          break;
        }
      }
      if (effort_ - progress_effort > lim_.stall_effort) {
        if (!lim_.stall_restart || found_) break;
        restart();
        progress_effort = effort_;
      }
      if (step_ > last_improve_ + lim_.restart_steps) restart();
      if (moves_since_refresh_ >= 100000) recompute();
      if (unsat_.empty()) {
        // The first feasible point is reported at once; later ones after the lift moves (each
        // report costs a pass over the matrix).
        if (!found_) {
          const int r = new_best(on_solution);
          if (r < 0) continue;  // drift: activities recomputed, not feasible after all
          if (r > 0) progress_effort = effort_;
          if (objvars_.empty()) break;  // feasibility problem: done
        }
        ++step_;
        if (lift_move()) continue;
        if (obj_ <= obj_bound_) {
          const int r = new_best(on_solution);
          if (r < 0) continue;
          if (r > 0) progress_effort = effort_;
        }
      } else if (static_cast<long long>(unsat_.size()) < min_unsat_) {
        min_unsat_ = static_cast<long long>(unsat_.size());
        if (solutions_ == 0) {
          progress_effort = effort_;
          if (stats) {
            stats->closest = x_;
            effort_ += n_;
          }
        }
      }
      explore();
      ++step_;
    }
    // Budget used in the middle of a chain of lift moves: report the point reached.
    const bool stopped = (lim_.stop && lim_.stop->load(std::memory_order_relaxed)) || interrupted();
    if (!stopped && found_ && unsat_.empty() && obj_ <= obj_bound_) new_best(on_solution);
    if (stats) {
      stats->effort = effort_;
      stats->steps = step_;
      stats->lift_moves = lifts_;
      stats->restarts = restarts_;
      stats->solutions = solutions_;
      stats->min_violated = solutions_ > 0 ? 0 : min_unsat_;
      stats->first_effort = first_effort_;
    }
    return solutions_ > 0;
  }

 private:
  double clamp_var(Index j, double v) const {
    v = std::clamp(v, p_.col_lb[j], p_.col_ub[j]);
    if (!std::isfinite(v)) v = 0.0;
    return p_.is_int[j] ? std::round(v) : v;
  }
  double improvement(double best) const {
    return obj_integral_ ? 1.0 - 1e-6 : std::max(1e-6, 1e-5 * std::abs(best));
  }

  // --- state -------------------------------------------------------------------------------
  void set_status(Index i) {
    const bool sat = act_[i] - b_[i] <= kTol;
    if (!sat && upos_[i] < 0) {
      remove(sat_, spos_, i);
      upos_[i] = static_cast<Index>(unsat_.size());
      unsat_.push_back(i);
    } else if (sat && upos_[i] >= 0) {
      remove(unsat_, upos_, i);
      spos_[i] = static_cast<Index>(sat_.size());
      sat_.push_back(i);
    }
  }
  static void remove(std::vector<Index>& list, std::vector<Index>& pos, Index i) {
    if (pos[i] < 0) return;
    const Index last = list.back();
    list[pos[i]] = last;
    pos[last] = pos[i];
    list.pop_back();
    pos[i] = -1;
  }
  void recompute() {
    act_.assign(m_, 0.0);
    obj_ = 0.0;
    for (Index j = 0; j < n_; ++j) {
      const double xj = x_[j];
      if (xj == 0.0) continue;
      for (Index k = cbeg_[j]; k < cbeg_[j + 1]; ++k) act_[crow_[k]] += cval_[k] * xj;
      obj_ += p_.cost[j] * xj;
    }
    effort_ += static_cast<long long>(crow_.size()) + n_;
    unsat_.clear();
    sat_.clear();
    upos_.assign(m_, -1);
    spos_.assign(m_, -1);
    for (Index i = 0; i < m_; ++i) {
      if (act_[i] - b_[i] > kTol) {
        upos_[i] = static_cast<Index>(unsat_.size());
        unsat_.push_back(i);
      } else {
        spos_[i] = static_cast<Index>(sat_.size());
        sat_.push_back(i);
      }
    }
    keep_feas_ = false;
    moves_since_refresh_ = 0;
  }
  void apply(Index j, double d) {
    if (d == 0.0) return;
    x_[j] += d;
    for (Index k = cbeg_[j]; k < cbeg_[j + 1]; ++k) {
      act_[crow_[k]] += cval_[k] * d;
      set_status(crow_[k]);
    }
    obj_ += p_.cost[j] * d;
    effort_ += 1 + cbeg_[j + 1] - cbeg_[j];
    const long long tenure = step_ + kTabuBase + static_cast<long long>(rng_() % kTabuVar);
    if (d > 0) allow_dec_[j] = tenure;
    else allow_inc_[j] = tenure;
    last_move_[j] = step_;
    last_j_ = j;
    last_d_ = d;
    ++moves_since_refresh_;
  }
  bool tabu(Index j, double d) const {
    return d < 0 ? step_ < allow_dec_[j] : step_ < allow_inc_[j];
  }
  // Only the reverse of the previous step (the last neighbourhood's weaker tabu).
  bool tabu_latest(Index j, double d) const {
    return j == last_j_ && step_ > 0 && last_move_[j] == step_ - 1 && (d > 0) != (last_d_ > 0);
  }

  // Tight move of column j (coefficient a) on a row with gap g = a^T x - b: the value that makes
  // the row tight (satisfied, for integers), clamped to the bounds (violated row: towards
  // satisfying it; satisfied row: towards using its slack).
  double tight(Index j, double a, double g, bool violated) const {
    double d;
    if (p_.is_int[j]) {
      const double t = (-g + kTol) / a;
      d = a > 0 ? std::floor(t + 1e-9) : std::ceil(t - 1e-9);
    } else {
      d = -g / a;
    }
    const double v = x_[j] + d;
    if (!(v >= p_.col_lb[j] && v <= p_.col_ub[j])) {
      const bool down = violated ? a > 0 : a < 0;
      d = (down ? p_.col_lb[j] : p_.col_ub[j]) - x_[j];
    }
    if (!std::isfinite(d) || std::abs(x_[j] + d) > kHuge) return 0.0;
    return d;
  }

  // --- moves -------------------------------------------------------------------------------
  void add_op(Index j, double d) {
    if (d == 0.0 || std::abs(d) < 1e-9) return;
    op_j_.push_back(j);
    op_d_.push_back(d);
  }
  void add_row_ops(Index i, bool latest_only) {
    const double g = act_[i] - b_[i];
    const bool violated = g > kTol;
    for (Index k = rbeg_[i]; k < rbeg_[i + 1]; ++k) {
      const Index j = rcol_[k];
      const double d = tight(j, rval_[k], g, violated);
      if (d == 0.0 || (latest_only ? tabu_latest(j, d) : tabu(j, d))) continue;
      add_op(j, d);
    }
    effort_ += rbeg_[i + 1] - rbeg_[i];
  }
  void add_breakthrough_ops(bool latest_only) {
    if (!found_ || obj_ <= obj_bound_ || objvars_.empty()) return;
    const Index no = static_cast<Index>(objvars_.size());
    const Index cnt = std::min<Index>(no, kBreakCols);
    const Index off = no > cnt ? static_cast<Index>(rng_() % static_cast<uint32_t>(no)) : 0;
    const double g = obj_ - obj_bound_;
    for (Index t = 0; t < cnt; ++t) {
      const Index j = objvars_[(off + t) % no];
      const double d = tight(j, p_.cost[j], g, true);
      if (d == 0.0 || (latest_only ? tabu_latest(j, d) : tabu(j, d))) continue;
      add_op(j, d);
    }
    effort_ += cnt;
  }
  void sample_ops(int k) {
    const size_t avail = op_j_.size();
    if (avail <= static_cast<size_t>(k)) return;
    for (size_t t = 0; t < static_cast<size_t>(k); ++t) {
      const size_t r = t + rng_() % (avail - t);
      std::swap(op_j_[t], op_j_[r]);
      std::swap(op_d_[t], op_d_[r]);
    }
    op_j_.resize(k);
    op_d_.resize(k);
    effort_ += k;
  }
  // Progress score with the breakthrough bonus and age as tie-breakers.
  void score(Index j, double d) {
    if (p_.is_int[j] && p_.col_lb[j] == 0.0 && p_.col_ub[j] == 1.0) {
      if (stamp_[j] == token_) return;  // a binary has one move
      stamp_[j] = token_;
    }
    long long s = 0, bonus = 0;
    for (Index k = cbeg_[j]; k < cbeg_[j + 1]; ++k) {
      const Index i = crow_[k];
      const double g = act_[i] - b_[i], ng = g + cval_[k] * d;
      const bool pre = g <= kTol, now = ng <= kTol;
      const long long w = w_[i];
      if (!pre && now) s += 2 * w;
      else if (pre && !now) s -= 2 * w;
      else if (!pre) s += ng < g ? w : -w;
    }
    effort_ += 1 + cbeg_[j + 1] - cbeg_[j];
    const double c = p_.cost[j];
    if (found_ && c != 0.0) {
      const double nobj = obj_ + c * d;
      s += nobj < obj_ ? 2 * w_obj_ : -2 * w_obj_;
      if (nobj < best_obj_) bonus = w_obj_;
    }
    const long long age = last_move_[j];
    if (s > best_s_ || (s == best_s_ && (bonus > best_bonus_ || (bonus == best_bonus_ && age < best_age_)))) {
      best_s_ = s;
      best_bonus_ = bonus;
      best_age_ = age;
      best_j_ = j;
      best_d_ = d;
    }
  }
  void score_ops() {
    for (size_t t = 0; t < op_j_.size(); ++t) score(op_j_[t], op_d_[t]);
  }
  void reset_best(bool positive) {
    if (++token_ == 0) {
      std::fill(stamp_.begin(), stamp_.end(), 0);
      token_ = 1;
    }
    best_s_ = positive ? 0 : LLONG_MIN;
    best_bonus_ = LLONG_MIN;
    best_age_ = LLONG_MAX;
    best_j_ = -1;
    best_d_ = 0.0;
  }
  // Distinct random entries of `list` (k small).
  void sample_rows(const std::vector<Index>& list, int k, std::vector<Index>& out) {
    out.clear();
    const size_t sz = list.size();
    if (sz <= static_cast<size_t>(k)) {
      out.assign(list.begin(), list.end());
      return;
    }
    while (static_cast<int>(out.size()) < k) {
      const Index i = list[rng_() % sz];
      if (std::find(out.begin(), out.end(), i) == out.end()) out.push_back(i);
    }
  }

  void explore() {
    reset_best(true);
    // 1. Mixed tight moves on a sample of violated rows, breakthrough moves.
    op_j_.clear();
    op_d_.clear();
    sample_rows(unsat_, bms().unsat_rows, rows_);
    for (Index i : rows_) add_row_ops(i, false);
    add_breakthrough_ops(false);
    sample_ops(bms().unsat_ops);
    score_ops();
    if (best_s_ > 0) return take();
    // 2. Tight moves on a satisfied (inequality) row: use its slack.
    if (found_ && !sat_.empty()) {
      op_j_.clear();
      op_d_.clear();
      sample_rows(sat_, bms().sat_rows, rows_);
      for (Index i : rows_) {
        if (!eq_half_[i]) add_row_ops(i, false);
      }
      sample_ops(bms().sat_ops);
      score_ops();
      if (best_s_ > 0) return take();
    }
    // 3. "Easy" moves: a few columns to 0 / a bound / the middle of their range.
    if (!nonfixed_.empty()) {
      op_j_.clear();
      op_d_.clear();
      for (int t = 0; t < bms().easy; ++t) {
        const Index j = nonfixed_[rng_() % nonfixed_.size()];
        const double lo = p_.col_lb[j], hi = p_.col_ub[j], xj = x_[j];
        auto easy = [&](double v) {
          const double d = v - xj;
          if (std::abs(v) <= kHuge && !tabu(j, d)) add_op(j, d);
        };
        easy(lo > 0.0 ? lo : hi < 0.0 ? hi : 0.0);
        if (!p_.is_int[j] && is_finite(lo) && is_finite(hi)) easy(0.5 * (lo + hi));
        if (is_finite(lo) && lo < 0.0) easy(lo);
        if (is_finite(hi) && hi > 0.0) easy(hi);
      }
      effort_ += bms().easy;
      score_ops();
      if (best_s_ > 0) return take();
    }
    // 4. Local minimum: bump weights, then the best move of one random violated row.
    update_weights();
    reset_best(false);
    op_j_.clear();
    op_d_.clear();
    if (!unsat_.empty()) add_row_ops(unsat_[rng_() % unsat_.size()], true);
    add_breakthrough_ops(true);
    sample_ops(bms().random_ops);
    score_ops();
    if (best_j_ >= 0) take();
  }
  void take() {
    apply(best_j_, best_d_);
    keep_feas_ = false;
  }
  void update_weights() {
    for (Index i : unsat_) ++w_[i];
    effort_ += static_cast<long long>(unsat_.size());
    if (found_ && unsat_.empty()) ++w_obj_;
  }

  // --- feasible mode -------------------------------------------------------------------------
  // Largest objective-improving step of objective column j that keeps every row satisfied.
  double lift_delta(Index j) {
    double lo = p_.col_lb[j] - x_[j], hi = p_.col_ub[j] - x_[j];
    for (Index k = cbeg_[j]; k < cbeg_[j + 1] && lo < hi; ++k) {
      const Index i = crow_[k];
      const double a = cval_[k], g = act_[i] - b_[i];
      if (g >= 0.0) {
        if (a > 0) hi = 0.0;
        else lo = 0.0;
        continue;
      }
      double d = -g / a;
      if (a > 0) {
        if (p_.is_int[j]) d = std::floor(d + 1e-9);
        hi = std::min(hi, d);
      } else {
        if (p_.is_int[j]) d = std::ceil(d - 1e-9);
        lo = std::max(lo, d);
      }
    }
    effort_ += 1 + cbeg_[j + 1] - cbeg_[j];
    const double d = p_.cost[j] > 0 ? lo : hi;
    if (!std::isfinite(d) || std::abs(x_[j] + d) > kHuge) return 0.0;
    return d;
  }
  bool lift_move() {
    if (objvars_.empty()) return false;
    const Index no = static_cast<Index>(objvars_.size());
    if (!keep_feas_) {
      for (Index t = 0; t < no; ++t) lift_delta_[t] = lift_delta(objvars_[t]);
    }
    const double thresh = 1e-9 * std::max(1.0, std::abs(obj_));
    Index bj = -1;
    double bs = thresh, bd = 0.0;
    long long bage = LLONG_MAX;
    for (Index t = 0; t < no; ++t) {
      const Index j = objvars_[t];
      const double s = -p_.cost[j] * lift_delta_[t];
      if (s > bs + 1e-12 || (bj >= 0 && s >= bs && last_move_[j] < bage)) {
        bj = j;
        bs = s;
        bd = lift_delta_[t];
        bage = last_move_[j];
      }
    }
    effort_ += no;
    if (bj < 0) {
      keep_feas_ = false;
      return false;
    }
    apply(bj, bd);
    ++lifts_;
    // Recompute the lift steps of the objective columns sharing a row with bj.
    ++token_;
    if (token_ == 0) {
      std::fill(stamp_.begin(), stamp_.end(), 0);
      token_ = 1;
    }
    lift_delta_[objpos_[bj]] = lift_delta(bj);
    stamp_[bj] = token_;
    for (Index k = cbeg_[bj]; k < cbeg_[bj + 1]; ++k) {
      const Index i = crow_[k];
      for (Index q = rbeg_[i]; q < rbeg_[i + 1]; ++q) {
        const Index j = rcol_[q];
        if (stamp_[j] == token_) continue;
        stamp_[j] = token_;
        if (objpos_[j] >= 0) lift_delta_[objpos_[j]] = lift_delta(j);
      }
      effort_ += rbeg_[i + 1] - rbeg_[i];
    }
    keep_feas_ = unsat_.empty();
    return true;
  }

  // All rows satisfied and better than the best: verify and report. Returns 1 if reported, 0 if
  // feasible but not reported (not below the cutoff), -1 if the activities had drifted.
  int new_best(const std::function<void(const std::vector<double>&, double)>& on_solution) {
    recompute();
    if (!unsat_.empty()) return -1;
    if (found_ && obj_ > obj_bound_) return 0;
    effort_ += static_cast<long long>(crow_.size());
    if (p_.max_violation(x_) > 1e-6) {
      // Tolerances disagree (should not happen): treat as violated, keep going.
      for (Index i = 0; i < m_; ++i) w_[i] += 1;
      return 0;
    }
    double obj = 0.0;
    for (Index j = 0; j < n_; ++j) obj += p_.cost[j] * x_[j];
    found_ = true;
    best_obj_ = obj;
    best_x_ = x_;
    obj_bound_ = obj - improvement(obj);
    last_improve_ = step_;
    if (obj < lim_.cutoff) {
      if (first_effort_ < 0) first_effort_ = effort_;
      ++solutions_;
      on_solution(x_, obj);
      return 1;
    }
    return 0;
  }

  void restart() {
    ++restarts_;
    if (!best_x_.empty()) {
      x_ = best_x_;
    } else {
      for (Index j = 0; j < n_; ++j) {
        const double lo = p_.col_lb[j], hi = p_.col_ub[j];
        if (!p_.is_int[j] || !is_finite(lo) || !is_finite(hi) || hi - lo > 1e9) continue;
        x_[j] = lo + static_cast<double>(rng_() % static_cast<uint64_t>(hi - lo + 1));
      }
    }
    std::fill(w_.begin(), w_.end(), 1);
    w_obj_ = 1;
    std::fill(allow_inc_.begin(), allow_inc_.end(), 0);
    std::fill(allow_dec_.begin(), allow_dec_.end(), 0);
    last_improve_ = step_;
    recompute();
  }

  const MipProblem& p_;
  LocalMipLimits lim_;
  std::mt19937 rng_;
  std::chrono::steady_clock::time_point t0_;
  Index n_ = 0, m_ = 0;
  std::vector<double> b_;
  std::vector<char> eq_half_;
  std::vector<Index> cbeg_, crow_, rbeg_, rcol_;
  std::vector<double> cval_, rval_;
  std::vector<Index> objvars_, objpos_, nonfixed_;
  bool obj_integral_ = false;
  // state
  std::vector<double> x_, act_;
  double obj_ = 0.0;
  std::vector<Index> unsat_, upos_, sat_, spos_;
  std::vector<long long> w_;
  long long w_obj_ = 1;
  std::vector<long long> allow_inc_, allow_dec_, last_move_;
  Index last_j_ = -1;  // column of the last move and its direction (tabu_latest)
  double last_d_ = 0.0;
  bool found_ = false, keep_feas_ = false;
  double best_obj_ = kInf, obj_bound_ = kInf;
  std::vector<double> best_x_;
  std::vector<double> lift_delta_;
  // move selection
  std::vector<Index> op_j_, rows_;
  std::vector<double> op_d_;
  std::vector<uint32_t> stamp_;
  uint32_t token_ = 0;
  long long best_s_ = 0, best_bonus_ = 0, best_age_ = 0;
  Index best_j_ = -1;
  double best_d_ = 0.0;
  // counters
  long long effort_ = 0, step_ = 0, last_improve_ = 0, lifts_ = 0, restarts_ = 0, min_unsat_ = 0;
  long long moves_since_refresh_ = 0, first_effort_ = -1;
  int solutions_ = 0;
};

}  // namespace

bool local_mip(const MipProblem& p, const LocalMipLimits& limits,
               const std::function<void(const std::vector<double>&, double)>& on_solution, LocalMipStats* stats) {
  if (p.num_cols() == 0) return false;
  LocalMip ls(p, limits);
  return ls.run(on_solution, stats);
}

bool local_search_bench(const MipProblem& p) {
  const char* v = std::getenv("MILP_LS_BENCH");
  if (v == nullptr) return false;
  const long long budget = std::max(1LL, std::atoll(v));
  const char* tv = std::getenv("MILP_LS_BENCH_TIME");
  const double tmax = tv ? std::atof(tv) : 120.0;
  const char* sv = std::getenv("MILP_LS_BENCH_SEED");
  const unsigned seed = sv ? static_cast<unsigned>(std::atoi(sv)) : 0u;
  const char* mv = std::getenv("MILP_LS_BENCH_METHODS");
  const std::string methods = mv ? mv : "fj,lmip";
  using Clock = std::chrono::steady_clock;
  auto secs = [](Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); };
  auto disp = [&](double obj) { return p.display_sense * (obj + p.offset); };
  std::printf("LSBENCH problem rows %d cols %d nnz %d budget %lld\n", static_cast<int>(p.num_rows()),
              static_cast<int>(p.num_cols()), static_cast<int>(p.A.nnz()), budget);
  if (methods.find("fj") != std::string::npos) {
    // FJ as the solver runs it, but restarted with new seeds until the whole budget is used.
    const Clock::time_point t0 = Clock::now();
    double best = kInf, first_t = -1.0;
    long long used = 0, first_e = -1, min_viol = -1;
    int sols = 0, attempts = 0;
    while (used < budget && secs(t0) < tmax) {
      FjLimits lim;
      lim.effort = budget - used;
      lim.stall_effort = 300000000LL;
      lim.seed = seed + 1000003u * static_cast<unsigned>(attempts++);
      lim.time_limit = tmax - secs(t0);
      FjStats st;
      const long long base = used;
      feasibility_jump(p, lim, [&](const std::vector<double>& x, double obj) {
        if (p.max_violation(x) > 1e-6) std::printf("LSBENCH ERROR fj infeasible point\n");
        if (first_t < 0) first_t = secs(t0);
        if (obj < best) {
          best = obj;
          ++sols;
        }
      }, &st);
      if (first_e < 0 && st.first_effort >= 0) first_e = base + st.first_effort;
      used += std::max(1LL, st.effort);
      if (st.solutions == 0 && (min_viol < 0 || st.min_violated < min_viol)) min_viol = st.min_violated;
    }
    std::printf("LSBENCH fj found %d first_t %.3f first_e %lld best %.10g sols %d effort %lld time %.3f "
                "minviol %lld attempts %d\n", best < kInf, first_t, first_e, best < kInf ? disp(best) : 0.0, sols, used,
                secs(t0), min_viol, attempts);
  }
  if (methods.find("lmip") != std::string::npos) {
    const Clock::time_point t0 = Clock::now();
    LocalMipLimits lim;
    lim.effort = budget;
    const char* st_env = std::getenv("MILP_LS_BENCH_STALL");
    lim.stall_effort = st_env ? std::atoll(st_env) : budget;
    lim.stall_restart = true;
    lim.seed = seed;
    lim.time_limit = tmax;
    const char* rs = std::getenv("MILP_LS_BENCH_RESTART");
    if (rs) lim.restart_steps = std::atoll(rs);
    double best = kInf, first_t = -1.0;
    int sols = 0;
    LocalMipStats st;
    local_mip(p, lim, [&](const std::vector<double>& x, double obj) {
      if (p.max_violation(x) > 1e-6) std::printf("LSBENCH ERROR lmip infeasible point\n");
      if (first_t < 0) first_t = secs(t0);
      best = obj;
      ++sols;
    }, &st);
    std::printf("LSBENCH lmip found %d first_t %.3f first_e %lld best %.10g sols %d effort %lld time %.3f "
                "minviol %lld steps %lld lifts %lld restarts %lld\n", best < kInf, first_t, st.first_effort,
                best < kInf ? disp(best) : 0.0, sols, st.effort, secs(t0), st.min_violated, st.steps, st.lift_moves,
                st.restarts);
  }
  std::fflush(stdout);
  return true;
}

}  // namespace milp
