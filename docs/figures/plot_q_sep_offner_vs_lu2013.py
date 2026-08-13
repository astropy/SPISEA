#!/usr/bin/env python
"""
Generate q / separation comparison figures vs Lu et al. (2013):

    docs/figures/q_offner_vs_lu2013.png
    docs/figures/sep_offner_vs_lu2013.png
    docs/figures/meanq_offner_vs_lu2013.png

Two-panel layout matching ``plot_mf_offner_vs_lu2013.py``: brown-dwarf
zoom and BD through O. Offner curves are evaluated from the
multiplicity objects so they cannot drift from the code.

Run from the repository root::

    python docs/figures/plot_q_sep_offner_vs_lu2013.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Allow running without installing the package.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from spisea.imf import multiplicity


# Table 1 γ_trunc (1–100 au / 1–102 au) with 1σ. Masses are geometric
# means of the tabulated M1 intervals. The L/early-T γ=2.5 text value
# is an interpolation knot only and is not plotted as a Table 1 point.
_TABLE1_GAMMA = [
    # name, M_lo, M_hi, gamma, gamma_err
    ('Fontanive+2018', 0.019, 0.058, 4.8, 2.2),
    ('Close+2003', 0.080, 0.095, 3.3, 1.2),
    ('Allen+2007', 0.06, 0.15, 1.7, 0.5),
    ('Winters mid-M', 0.15, 0.30, 0.7, 0.5),
    ('Winters early-M', 0.3, 0.6, 0.1, 0.4),
    ('Raghavan+2010', 0.75, 1.25, 0.2, 0.4),
    ('De Rosa A', 1.6, 2.4, -1.3, 0.4),
    ('MDS 3-5', 3.0, 5.0, -1.0, 0.5),
    ('MDS 5-8', 5.0, 8.0, -1.7, 0.5),
    ('MDS 8-17', 8.0, 17.0, -1.6, 0.5),
    ('Sana O', 17.0, 50.0, -1.4, 0.4),
]

# Table 1 ã_all (au) with 1σ.
_TABLE1_A_ALL = [
    ('Fontanive+2018', 0.019, 0.058, 2.9, 1.1),
    ('Close+2003', 0.080, 0.095, 3.7, 1.3),
    ('Allen+2007', 0.06, 0.15, 6.9, 1.4),
    ('Winters late-M', 0.075, 0.15, 3.9, 1.2),
    ('Winters mid-M', 0.15, 0.30, 10.0, 3.0),
    ('Winters early-M', 0.3, 0.6, 26.0, 4.0),
    ('Raghavan+2010', 0.75, 1.25, 49.0, 6.0),
    ('Tokovinin 2014b', 0.85, 1.5, 31.0, 5.0),
    ('Moe & Kratter', 1.6, 2.4, 32.0, 8.0),
    ('MDS 3-5', 3.0, 5.0, 28.0, 7.0),
    ('MDS 5-8', 5.0, 8.0, 25.0, 7.0),
    ('MDS 8-17', 8.0, 17.0, 23.0, 7.0),
    ('Sana O', 17.0, 50.0, 19.0, 6.0),
]

# Table 2 lognormal μ (au) at the three published bins.
_TABLE2_MU = [
    ('late-M', 0.075, 0.15, 4.0),
    ('early-M', 0.3, 0.6, 25.0),
    ('FGK', 0.75, 1.25, 40.0),
]

_OFFNER_COLOR = '#8b3a2a'
_TABLE1_COLOR = '#2f6db3'
_BD_SHADE = '#e8d5b5'
_FIGSIZE = (11.2, 4.6)
_DPI = 160
_BD_XLIM = (0.012, 0.20)
_FULL_XLIM = (0.015, 40.0)


def _geom(lo, hi):
    return float(np.sqrt(lo * hi))


def _table_xy(rows, y_idx=3, e_idx=4):
    m = np.array([_geom(r[1], r[2]) for r in rows])
    y = np.array([r[y_idx] for r in rows], dtype=float)
    err = np.array([r[e_idx] for r in rows], dtype=float)
    return m, y, err


def _offner_mean_a_au(resolved, mass):
    """
    Characteristic μ(a) in AU: same log-M interpolation as
    ``MultiplicityResolvedOffner2023.log_semimajoraxis``, without drawing.
    """
    mass = np.atleast_1d(np.asarray(mass, dtype=float))
    mass_clip = np.clip(mass, resolved.sep_mass[0], resolved.sep_mass[-1])
    log_a_mean = np.interp(
        np.log10(mass_clip),
        np.log10(resolved.sep_mass),
        np.log10(resolved.sep_mu_au),
    )
    return 10.0 ** log_a_mean


def _lu2013_dk_mean_a_au(dk, mass):
    """
    Duchêne & Kraus mean a in AU from ``MultiplicityResolvedDK`` coefficients,
    plus the BD log-a interpolation and sigmoid blend in ``log_semimajoraxis``.
    """
    mass = np.atleast_1d(np.asarray(mass, dtype=float))
    logm = np.log10(mass)
    x = mass / dk.a_break
    a_star = np.where(
        mass < dk.a_break,
        dk.a_amp * np.power(x, -dk.a_slope1),
        dk.a_amp * np.power(x, -dk.a_slope2),
    )
    a_star = np.maximum(a_star, 1e-30)
    log_a_mean_star = np.log10(a_star)
    log_a_mean_bd = np.interp(
        logm,
        [np.log10(0.01), np.log10(0.08)],
        [np.log10(2.5), np.log10(8.0)],
    )
    w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
    log_a_mean = (1.0 - w) * log_a_mean_bd + w * log_a_mean_star
    return 10.0 ** log_a_mean


def _mean_q_from_gamma(gamma, q_min=0.01):
    """⟨q⟩ for P(q) ∝ q^γ on [q_min, 1]."""
    g = np.asarray(gamma, dtype=float)
    qmin = float(q_min)
    out = np.empty(np.shape(g), dtype=float)
    g_flat = np.atleast_1d(g).astype(float)
    out_flat = np.empty(g_flat.shape, dtype=float)
    near_m1 = np.abs(g_flat + 1.0) < 1e-12
    near_m2 = np.abs(g_flat + 2.0) < 1e-12
    ok = ~near_m1 & ~near_m2
    if np.any(near_m1):
        out_flat[near_m1] = (1.0 - qmin) / (-np.log(qmin))
    if np.any(near_m2):
        out_flat[near_m2] = -np.log(qmin) / (1.0 / qmin - 1.0)
    if np.any(ok):
        gp = g_flat[ok]
        num = (1.0 - np.power(qmin, gp + 2.0)) / (gp + 2.0)
        den = (1.0 - np.power(qmin, gp + 1.0)) / (gp + 1.0)
        out_flat[ok] = num / den
    out = out_flat.reshape(np.shape(g))
    return float(out) if np.isscalar(gamma) else out


def _bd_shade(ax, xlim):
    ax.axvspan(xlim[0], 0.08, color=_BD_SHADE, alpha=0.55, zorder=0)
    ax.axvline(0.08, color='#c4a574', ls='--', lw=1.2, zorder=1)


def _finish_panel(ax, xlim, ylim, title, ylabel, ylog=False):
    ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r'Primary mass $M_1$ ($M_\odot$)')
    ax.set_ylabel(ylabel)
    ax.tick_params(which='both', direction='in', top=True, right=True)


def _two_axes(suptitle):
    fig, axes = plt.subplots(1, 2, figsize=_FIGSIZE, gridspec_kw={'wspace': 0.28})
    fig.suptitle(suptitle, fontsize=13, y=1.02)
    return fig, axes


def _save(fig, filename):
    out = os.path.join(os.path.dirname(__file__), filename)
    fig.savefig(out, dpi=_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Wrote', out)


def _lu_gamma_step_masses():
    """Dense sampling so the 0.08 Msun γ step renders as a vertical jump."""
    return np.concatenate([
        np.logspace(np.log10(0.012), np.log10(0.07999), 300),
        np.array([0.08, 0.08001]),
        np.logspace(np.log10(0.081), np.log10(40.0), 300),
    ])


def plot_gamma(offner, lu):
    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    g_off = offner.q_power_at_mass(m_wide)
    m_step = _lu_gamma_step_masses()
    g_lu = lu.q_power_at_mass(m_step)
    m_tab, g_tab, err_tab = _table_xy(_TABLE1_GAMMA)

    fig, axes = _two_axes(r'Mass-ratio index $\gamma$ vs primary mass')
    ylabel = r'$\gamma$  ($P(q)\propto q^{\gamma}$)'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        _bd_shade(ax, xlim)
        ax.axhline(0.0, color='0.55', ls=':', lw=1.1, zorder=2)
        ax.plot(m_step, g_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, g_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        ax.errorbar(m_tab, g_tab, yerr=err_tab, fmt='o', color=_TABLE1_COLOR,
                    ms=5.5, mfc='white', mew=1.3, elinewidth=1.1, capsize=2.5,
                    zorder=5)
        _finish_panel(ax, xlim, (-2.2, 7.0), title, ylabel)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'Lu+2013  $\gamma=6.1$ (BD), $-0.4$ (stellar)'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner Table 1 $\gamma_\mathrm{trunc}$ interp'),
        Line2D([0], [0], marker='o', color=_TABLE1_COLOR, ls='none',
               mfc='white', mew=1.3, ms=6, label=r'Table 1 $\gamma_\mathrm{trunc}$'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    _save(fig, 'q_offner_vs_lu2013.png')


def plot_sep(resolved, dk):
    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    a_off = _offner_mean_a_au(resolved, m_wide)
    a_lu = _lu2013_dk_mean_a_au(dk, m_wide)
    m_tab, a_tab, err_tab = _table_xy(_TABLE1_A_ALL)
    m_t2 = np.array([_geom(r[1], r[2]) for r in _TABLE2_MU])
    a_t2 = np.array([r[3] for r in _TABLE2_MU], dtype=float)

    fig, axes = _two_axes(r'Characteristic separation vs primary mass')
    ylabel = r'$\mu(a)$ (AU)'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        _bd_shade(ax, xlim)
        ax.plot(m_wide, a_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, a_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        ax.errorbar(m_tab, a_tab, yerr=err_tab, fmt='o', color=_TABLE1_COLOR,
                    ms=5.5, mfc='white', mew=1.3, elinewidth=1.1, capsize=2.5,
                    zorder=5)
        ax.plot(m_t2, a_t2, 's', color=_OFFNER_COLOR, mfc=_OFFNER_COLOR,
                ms=7, zorder=6, mew=0.6)
        _finish_panel(ax, xlim, (1.0, 400.0), title, ylabel, ylog=True)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'Lu+2013 DK mean $a$'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner $\mu(a)$ (Table 1/2 interp)'),
        Line2D([0], [0], marker='o', color=_TABLE1_COLOR, ls='none',
               mfc='white', mew=1.3, ms=6, label=r'Table 1 $\tilde{a}_\mathrm{all}$'),
        Line2D([0], [0], marker='s', color=_OFFNER_COLOR, ls='none',
               mfc=_OFFNER_COLOR, ms=7, label=r'Table 2 $\mu$ knots'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[0].legend(handles=legend_handles, loc='upper left', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    _save(fig, 'sep_offner_vs_lu2013.png')


def plot_meanq(offner, lu):
    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    q_off = _mean_q_from_gamma(offner.q_power_at_mass(m_wide),
                              q_min=offner.q_min)
    m_step = _lu_gamma_step_masses()
    q_lu = _mean_q_from_gamma(lu.q_power_at_mass(m_step), q_min=lu.q_min)

    fig, axes = _two_axes(r'Mean mass ratio $\langle q\rangle$ vs primary mass')
    ylabel = r'$\langle q\rangle$  on $[0.01,\,1]$'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        _bd_shade(ax, xlim)
        ax.plot(m_step, q_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, q_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        _finish_panel(ax, xlim, (0.0, 1.0), title, ylabel)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'Lu+2013  from $\gamma$ step'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner  from $\gamma_\mathrm{trunc}(M)$'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    _save(fig, 'meanq_offner_vs_lu2013.png')


def main():
    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()
    resolved = multiplicity.MultiplicityResolvedOffner2023()
    dk = multiplicity.MultiplicityResolvedDK()
    plot_gamma(offner, lu)
    plot_sep(resolved, dk)
    plot_meanq(offner, lu)


if __name__ == '__main__':
    main()
