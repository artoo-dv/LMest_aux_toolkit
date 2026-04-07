# Lmest Auxillary functions

This document proposes a set of additions to the LMest package to support
convergence diagnostics, model output interpretation, and parallel estimation.
All functions are illustrated using the built-in `PSIDlong` dataset. Each
proposed function is self-contained: it takes a fitted `lmest` model object as
its primary input and returns results consistent with existing LMest conventions.
