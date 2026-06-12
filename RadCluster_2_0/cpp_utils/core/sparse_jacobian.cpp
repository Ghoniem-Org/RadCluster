/**
 * sparse_jacobian.cpp — Colored finite-difference sparse Jacobian for KLU.
 *
 * See sparse_jacobian.h for interface notes.
 */
#include "sparse_jacobian.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Sparsity pattern
// ─────────────────────────────────────────────────────────────────────────────
//
// State layout (full_CD, see parameters.h):
//   y[0..I-1]            SIA c_n,           n=1..I
//   y[I..I+V-1]          VAC c_m,           m=1..V
//   y[I+V..N_eq-6]       He extras:
//       Case 2 (he_mode=0): Q_tot at I+V, c_h at I+V+1 (dynamic)
//       Case 1 (he_mode=1): Q_m at I+V..I+2V-1, c_h at I+2V (dynamic)
//   y[N_eq-5..N_eq-1]    conservation accounting (J_SIA_fixed, J_SIA_mutual,
//                        J_VAC_fixed, J_VAC_mutual, J_He_sink)
//
// Mobile-cluster columns couple to many rows via global sums and are treated
// as dense columns:
//   SIA mobile : 0..i_mobile-1
//   VAC mobile : I..I+v_mobile-1
//   He extras  : I+V..N_eq-6
//
// Non-mobile cluster columns have local (banded) row coupling driven by
// coalescence/annihilation/emission. Conservation columns have no incoming
// coupling (only the diagonal is stored, KLU requires a structurally non-
// zero diagonal).

int build_sparsity_pattern_full_CD(Parameters& P)
{
    if (P.physics_option != 0 && P.physics_option != 1) return -1;

    const int I       = P.I;
    const int V       = P.V;
    const int N_eq    = P.N_eq;
    const int i_mob   = std::max(1, P.i_mobile);
    const int v_mob   = std::max(1, P.v_mobile);
    const int he_off  = I + V;
    // cons_off is stable under the optional appended ⟨100⟩ block (P.cons_off =
    // N_eq-5 when loop conversion is off).  n_he is the He-extra count between
    // the VAC block and the conservation block.
    const int cons_off = P.cons_off;                // 5 conservation entries
    const int n_he    = cons_off - he_off;          // Q vars + (maybe) c_h

    // Build pattern as list of nonzero rows per column (CSC).
    std::vector<std::vector<int>> col_rows(N_eq);

    auto add = [&](int row, int col) {
        if (row < 0 || row >= N_eq || col < 0 || col >= N_eq) return;
        col_rows[col].push_back(row);
    };

    // Mark column j as a global "dense" column: every row depends on it.
    auto mark_dense_col = [&](int col) {
        for (int r = 0; r < N_eq; ++r) add(r, col);
    };

    // ── Mobile SIA columns (0..i_mobile-1): dense ───────────────────────────
    for (int j = 0; j < std::min(i_mob, I); ++j) mark_dense_col(j);

    // ── Mobile VAC columns (I..I+v_mobile-1): dense ─────────────────────────
    for (int j = I; j < I + std::min(v_mob, V); ++j) mark_dense_col(j);

    // ── He-extra columns (Q_tot, Q_m, c_h): dense ───────────────────────────
    for (int j = he_off; j < he_off + n_he; ++j) mark_dense_col(j);

    // ── SIA rows: local band coupling among non-mobile SIA columns ─────────
    // Row r = n - 1 (size n = 1..I).  Couples to:
    //   self
    //   n - 1 (thermal emission predecessor) and n + 1 (successor)
    //   coalescence reach: backward by i_mobile (gain), forward 0
    //   V-I annihilation reach: forward by v_mobile (gain c_{m'} c_{n+m'})
    // Non-mobile cluster columns are added here; mobile cluster cols and
    // VAC band cols are handled below or via mark_dense_col above.
    for (int r = 0; r < I; ++r) {
        const int j_lo = std::max(0, r - i_mob);
        const int j_hi = std::min(I - 1, r + v_mob);
        for (int j = j_lo; j <= j_hi; ++j) {
            if (j < i_mob) continue;            // already dense
            add(r, j);
        }
        // Thermal emission band of width 1 beyond the above (defensive)
        if (r - 1 >= i_mob) add(r, r - 1);
        if (r + 1 < I)      add(r, r + 1);

        // Mobile SIA rows additionally couple to all VAC cols (cluster–
        // cavity capture) — already covered: VAC mobile cols are dense and
        // the V-I annihilation gain reaches non-mobile VAC cols below.
        if (r < i_mob) {
            for (int j = I + v_mob; j < I + V; ++j) add(r, j);
        }
    }

    // ── VAC rows: local band coupling among non-mobile VAC columns ─────────
    // Row r = I + m - 1 (size m = 1..V). Couples to:
    //   self
    //   V-V coalescence backward by v_mobile (gain) and forward 0
    //   SIA-induced shrinkage gain forward by i_mobile (c_n c_{m+n})
    //   thermal emission ±1
    for (int m = 1; m <= V; ++m) {
        const int r = I + m - 1;
        const int m_lo = std::max(1, m - v_mob);
        const int m_hi = std::min(V, m + i_mob);
        for (int mm = m_lo; mm <= m_hi; ++mm) {
            const int j = I + mm - 1;
            if (j < I + v_mob) continue;        // already dense
            add(r, j);
        }
        if (m - 1 >= 1 + v_mob) add(r, I + m - 2);
        if (m + 1 <= V)         add(r, I + m);
    }

    // ── He-extra rows: depend on all VAC columns and all He extras ──────────
    // (Mobile SIA cols and self are already dense.)
    for (int r = he_off; r < he_off + n_he; ++r) {
        for (int j = I; j < I + V; ++j) add(r, j);
        for (int j = he_off; j < he_off + n_he; ++j) add(r, j);
    }

    // ── Conservation rows (the 5 at cons_off): aggregates over all clusters ──
    // Bounded by cons_off+5 (NOT N_eq): the appended ⟨100⟩ block sits after the
    // conservation block, so r<N_eq would spill across all ⟨100⟩ rows.
    for (int r = cons_off; r < cons_off + 5; ++r) {
        for (int j = 0; j < I; ++j)             add(r, j);   // SIA
        for (int j = I; j < I + V; ++j)         add(r, j);   // VAC
        for (int j = he_off; j < he_off + n_he; ++j) add(r, j);
    }

    // ── Loop-conversion ⟨100⟩ block (appended at sia100_off) ────────────────
    // Mobile SIA columns (0..i_mob-1) and mobile VAC columns (I..I+v_mob-1) are
    // already dense over ALL rows (incl. the appended ⟨100⟩ rows), which covers
    // the junction / absorption / growth / shrink dependencies on the mobile
    // ½⟨111⟩ and vacancy monomer pools.  Here we add the remaining structure:
    //   ⟨100⟩ self band (ladders ±1, absorption product reach +i_mob backward),
    //   the unary same-size ½⟨111⟩→⟨100⟩ diagonal coupling, and the ⟨100⟩
    //   columns feeding the monomer / mobile-½⟨111⟩ / conservation rows.
    if (P.loop_conversion) {
        const int s0    = P.sia100_off;
        const int i_mob_c = std::min(i_mob, I);
        for (int r = 0; r < I; ++r) {
            const int r100 = s0 + r;                 // ⟨100⟩ row, size r+1
            add(r100, r100);                         // self (diagonal)
            if (r - 1 >= 0) add(r100, s0 + (r - 1)); // growth predecessor
            if (r + 1 < I)  add(r100, s0 + (r + 1)); // shrink/emission predecessor
            for (int n = 1; n <= i_mob && r - n >= 0; ++n)
                add(r100, s0 + (r - n));             // absorption product reach
            add(r100, r);                            // unary: dci100[r] <- c_i[r]

            // ⟨100⟩ column r100 feeds: ½⟨111⟩ & vacancy monomer rows, the mobile
            // ½⟨111⟩ rows (absorption consumes them), and the conservation rows.
            add(0, r100);                            // ½⟨111⟩ monomer row
            add(I, r100);                            // vacancy monomer row
            for (int n = 0; n < i_mob_c; ++n) add(n, r100);
            for (int rc = cons_off; rc < cons_off + 5; ++rc) add(rc, r100);
        }
    }

    // ── Diagonal (structural, KLU requirement) ──────────────────────────────
    for (int j = 0; j < N_eq; ++j) add(j, j);

    // ── Deduplicate + sort each column's rows (CSC requires row indices
    //    sorted within each column) ───────────────────────────────────────
    int nnz = 0;
    for (int j = 0; j < N_eq; ++j) {
        auto& v = col_rows[j];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        nnz += static_cast<int>(v.size());
    }

    // ── Flatten into CSC arrays ────────────────────────────────────────────
    P.jac_col_ptr.assign(N_eq + 1, 0);
    P.jac_row_idx.assign(nnz, 0);
    int pos = 0;
    for (int j = 0; j < N_eq; ++j) {
        P.jac_col_ptr[j] = pos;
        for (int r : col_rows[j]) P.jac_row_idx[pos++] = r;
    }
    P.jac_col_ptr[N_eq] = pos;

    return nnz;
}

// ─────────────────────────────────────────────────────────────────────────────
// Greedy CPR coloring.  Two columns share a color iff they share no row.
// Greedy first-fit ordering by descending column degree gives a near-optimal
// number of colors for banded + dense-border patterns.
// ─────────────────────────────────────────────────────────────────────────────
void color_columns_greedy(Parameters& P)
{
    const int N = static_cast<int>(P.jac_col_ptr.size()) - 1;
    P.jac_colors.assign(N, -1);

    // For each row, build the list of columns that hit it (CSR-like).
    int nnz = P.jac_col_ptr[N];
    std::vector<int> row_count(N, 0);
    for (int k = 0; k < nnz; ++k) row_count[P.jac_row_idx[k]]++;
    std::vector<int> row_ptr(N + 1, 0);
    for (int i = 0; i < N; ++i) row_ptr[i + 1] = row_ptr[i] + row_count[i];
    std::vector<int> row_cols(nnz);
    std::vector<int> row_fill(N, 0);
    for (int j = 0; j < N; ++j) {
        for (int k = P.jac_col_ptr[j]; k < P.jac_col_ptr[j + 1]; ++k) {
            int i = P.jac_row_idx[k];
            row_cols[row_ptr[i] + row_fill[i]++] = j;
        }
    }

    // Process columns in descending-degree order so dense columns get colored
    // first (each takes its own color), and band columns share colors after.
    std::vector<int> order(N);
    for (int j = 0; j < N; ++j) order[j] = j;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        int da = P.jac_col_ptr[a + 1] - P.jac_col_ptr[a];
        int db = P.jac_col_ptr[b + 1] - P.jac_col_ptr[b];
        return da > db;
    });

    std::vector<int> forbidden;     // colors used by neighbours of current col
    int n_colors = 0;
    for (int j : order) {
        forbidden.assign(n_colors, 0);
        for (int k = P.jac_col_ptr[j]; k < P.jac_col_ptr[j + 1]; ++k) {
            int i = P.jac_row_idx[k];
            // Any column that shares row i with j is a "neighbour" in the
            // intersection graph — its color is forbidden for j.
            for (int p = row_ptr[i]; p < row_ptr[i + 1]; ++p) {
                int jj = row_cols[p];
                int c  = P.jac_colors[jj];
                if (c >= 0) {
                    if (c >= static_cast<int>(forbidden.size()))
                        forbidden.resize(c + 1, 0);
                    forbidden[c] = 1;
                }
            }
        }
        int chosen = -1;
        for (int c = 0; c < static_cast<int>(forbidden.size()); ++c)
            if (!forbidden[c]) { chosen = c; break; }
        if (chosen < 0) chosen = n_colors++;
        P.jac_colors[j] = chosen;
        if (chosen >= n_colors) n_colors = chosen + 1;
    }
    P.jac_n_colors = n_colors;
}

void build_color_groups(Parameters& P)
{
    const int N  = static_cast<int>(P.jac_colors.size());
    const int nc = P.jac_n_colors;

    P.jac_color_offsets.assign(nc + 1, 0);
    for (int j = 0; j < N; ++j) P.jac_color_offsets[P.jac_colors[j] + 1]++;
    for (int c = 0; c < nc; ++c)
        P.jac_color_offsets[c + 1] += P.jac_color_offsets[c];

    P.jac_color_groups.assign(N, 0);
    std::vector<int> fill(nc, 0);
    for (int j = 0; j < N; ++j) {
        int c = P.jac_colors[j];
        P.jac_color_groups[P.jac_color_offsets[c] + fill[c]++] = j;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sparse Jacobian via colored finite differences.
//
// For each color group, build a perturbation vector with one increment per
// column in the group. Within a group no two columns share a row, so the
// resulting (f(y+δ) - f(y))_i corresponds unambiguously to the unique column
// of that color which has i in its sparsity pattern.
// ─────────────────────────────────────────────────────────────────────────────
int sparse_fd_jac(sunrealtype /*t*/, N_Vector y, N_Vector fy,
                  SUNMatrix Jac, void* user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    UserData* ud = static_cast<UserData*>(user_data);
    Parameters& P = *(ud->P);
    const int  N  = P.N_eq;

    double* y_p   = N_VGetArrayPointer(y);
    double* f0    = N_VGetArrayPointer(fy);
    double* y_pert = N_VGetArrayPointer(tmp1);
    double* f1    = N_VGetArrayPointer(tmp2);
    double* delta = N_VGetArrayPointer(tmp3);

    // ── CSC arrays of the SUNSparseMatrix ───────────────────────────────────
    // We assume the sparsity pattern is already installed; only fill values.
    sunindextype* col_ptr = SUNSparseMatrix_IndexPointers(Jac);
    sunindextype* row_idx = SUNSparseMatrix_IndexValues  (Jac);
    sunrealtype*  data    = SUNSparseMatrix_Data         (Jac);

    SUNMatZero(Jac);

    // Per-column FD perturbation size, scaled with current state.
    const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());

    for (int c = 0; c < P.jac_n_colors; ++c) {
        // Reset perturbation vector to zero.
        std::fill(delta, delta + N, 0.0);

        // Build delta for all columns in this color group.
        // FD step is scaled by max(|y|, atol): using atol (not an unrelated
        // 1e-12 floor) avoids catastrophic cancellation for near-floor columns
        // and keeps this consistent with the Woodbury preconditioner's
        // fd_delta() in rate_equations.cpp.
        const int g_lo = P.jac_color_offsets[c];
        const int g_hi = P.jac_color_offsets[c + 1];
        for (int g = g_lo; g < g_hi; ++g) {
            int j = P.jac_color_groups[g];
            double scale = std::max(std::abs(y_p[j]),
                                    std::max(P.atol, 1e-30));
            delta[j] = sqrt_eps * scale;
        }

        // y_pert = y + delta
        for (int i = 0; i < N; ++i) y_pert[i] = y_p[i] + delta[i];

        // f1 = f(y_pert)
        N_Vector y_pert_vec = tmp1;
        N_Vector f1_vec     = tmp2;
        if (ud->rhs_fn(0.0, y_pert_vec, f1_vec, user_data) != 0)
            return -1;

        // For each column in the group, fill its sparse entries.
        for (int g = g_lo; g < g_hi; ++g) {
            int j = P.jac_color_groups[g];
            double inv_d = 1.0 / delta[j];
            for (sunindextype k = col_ptr[j]; k < col_ptr[j + 1]; ++k) {
                int i = static_cast<int>(row_idx[k]);
                data[k] = (f1[i] - f0[i]) * inv_d;
            }
        }
    }

    return 0;
}
