// Function bodies for common/flow.hpp.

#include "common/flow.hpp"

#include <cmath>

namespace {

constexpr double pi = M_PI;

}  // namespace

FlowCoeffs make_flow_coeffs(double nu, double dx, double dy) {
    return FlowCoeffs{1 / (12 * dx * dy), nu / (dx * dx), nu / (dy * dy)};
}

std::vector<double> wind_profile(const Config& c, const Block& b) {
    const double dy = c.ly / c.ny;
    std::vector<double> profile(static_cast<std::size_t>(b.ny) + 1);
    for (int i = 0; i <= b.ny; ++i) {
        const double y = corner_y(b.row_offset + i, dy, c.ly);
        profile[static_cast<std::size_t>(i)] = std::sin(2 * pi * y / c.ly);
    }
    return profile;
}

double wind_strength(const Config& c, double t) {
    if (c.wind == Wind::steady) return c.wind_amplitude;
    const double omega = 2 * pi / c.wind_period;
    return c.wind_amplitude * (1 + c.wind_alpha * std::sin(omega * t));
}

// The compass neighbours of corner (i, j): N is the row above (i - 1), E the
// column to the right (j + 1).
void vorticity_rhs_row(const Field& psi, const Field& zeta, Field& rhs, const FlowCoeffs& k, double forcing,
                       const double* extra, int i, int j0, int j1, bool advection) {
    const double* pn = psi.row(i - 1);
    const double* p = psi.row(i);
    const double* ps = psi.row(i + 1);
    const double* zn = zeta.row(i - 1);
    const double* z = zeta.row(i);
    const double* zs = zeta.row(i + 1);
    double* out = rhs.row(i);
    for (int j = j0; j <= j1; ++j) {
        const double pE = p[j + 1], pW = p[j - 1], pN = pn[j], pS = ps[j];
        const double pNE = pn[j + 1], pNW = pn[j - 1], pSE = ps[j + 1], pSW = ps[j - 1];
        const double zE = z[j + 1], zW = z[j - 1], zN = zn[j], zS = zs[j];
        const double zNE = zn[j + 1], zNW = zn[j - 1], zSE = zs[j + 1], zSW = zs[j - 1];

        double value = forcing + k.nu_dx2 * (zE - 2 * z[j] + zW) + k.nu_dy2 * (zN - 2 * z[j] + zS);
        if (advection) {
            // Arakawa Jacobian: the average of J++, J+x and Jx+.
            const double jpp = (pE - pW) * (zN - zS) - (pN - pS) * (zE - zW);
            const double jpx = pE * (zNE - zSE) - pW * (zNW - zSW) - pN * (zNE - zNW) + pS * (zSE - zSW);
            const double jxp = zN * (pNE - pNW) - zS * (pSE - pSW) - zE * (pNE - pSE) + zW * (pNW - pSW);
            value += (jpp + jpx + jxp) * k.inv_12dxdy;
        }
        if (extra) value += extra[j];
        out[j] = value;
    }
}

AbWeights ab_weights(long s) {
    if (s == 0) return AbWeights{1.0, 0.0, 0.0};
    if (s == 1) return AbWeights{1.5, -0.5, 0.0};
    return AbWeights{23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0};
}

void ab_update_row(Field& zeta, const Field& r0, const Field& r1, const Field& r2, const AbWeights& w, double dt,
                   int i, int j0, int j1) {
    double* z = zeta.row(i);
    const double* a = r0.row(i);
    const double* b = r1.row(i);
    const double* c = r2.row(i);
    for (int j = j0; j <= j1; ++j) {
        z[j] = z[j] + dt * (w.w0 * a[j] + w.w1 * b[j] + w.w2 * c[j]);
    }
}

void faces_from_psi_row(const Field& psi, FaceVelocities& vel, double dx, double dy, int i) {
    const double* above = psi.row(i - 1);
    const int nx = vel.nx();
    if (i <= vel.ny()) {
        const double* here = psi.row(i);
        double* u = vel.east_row(i);
        for (int j = 0; j <= nx; ++j) u[j] = (above[j] - here[j]) / dy;
    }
    double* v = vel.north_row(i);
    for (int j = 1; j <= nx; ++j) v[j - 1] = -(above[j] - above[j - 1]) / dx;
}

void faces_from_psi(const Field& psi, FaceVelocities& vel, double dx, double dy) {
    for (int i = 1; i <= vel.ny() + 1; ++i) faces_from_psi_row(psi, vel, dx, dy, i);
}

void local_energy_enstrophy(const Field& psi, const Field& zeta, const Block& b, double dx, double dy,
                            double& energy, double& enstrophy) {
    const CornerRange r = unknown_corners(b);
    double e = 0.0;
    double z2 = 0.0;
    for (int i = r.i_first; i <= r.i_last; ++i) {
        const double* p = psi.row(i);
        const double* z = zeta.row(i);
        for (int j = r.j_first; j <= r.j_last; ++j) {
            e += p[j] * z[j];
            z2 += z[j] * z[j];
        }
    }
    energy = 0.5 * dx * dy * e;
    enstrophy = 0.5 * dx * dy * z2;
}
