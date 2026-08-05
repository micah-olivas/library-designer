"""Integrations with external, separately-installed tools.

Each integration runs its tool as a subprocess over files, so
library_designer never imports or bundles the other tool's code. This keeps
library_designer's own dependency set lean and its MIT license unaffected by tools
that carry a different license (see ``omega``).
"""
