"""Drug/target evidence: Open Targets and ChEMBL clients and prioritisation logic.

Provides typed clients for the Open Targets GraphQL API and ChEMBL REST API,
and logic for converting raw tractability/bioactivity evidence into the
druggability evidence schema consumed by ``oncocartograph.scoring``.
"""
