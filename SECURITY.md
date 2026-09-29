# Security policy

Report vulnerabilities privately to the repository maintainers rather than opening a public
issue. The modern package never evaluates configuration as Python and never deletes external
GIS workspaces. CI has read-only repository permission, uses no project secrets, pins actions
and the GRASS image by immutable digest, and runs dependency and Bandit audits.

Legacy code is archival and contains known unsafe patterns (`eval`, `shell=True`, and recursive
workspace deletion). Do not execute it on untrusted input or outside a disposable environment.
