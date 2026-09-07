"""psmgrating -- PSM, PSM*, Kogelnik and TCWT two-wave models and a rigorous
coupled-wave reference for slanted, absorbing planar volume gratings, in the
conventions of Brotherton-Ratcliffe and Vivian Sureshkumar (2026).

Independent implementation validated against the paper's printed values; the
authors' own deposit (release v1.1) is a separate code base.
"""
from . import conventions, twowave, rcwa, chain, metrics, presets  # noqa: F401
from .twowave import efficiency, sweep, MODELS, LABELS  # noqa: F401
from .rcwa import sweep as rcwa_sweep  # noqa: F401

__version__ = "0.1.0"
