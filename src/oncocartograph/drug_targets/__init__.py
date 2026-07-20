"""Drug/target evidence: Open Targets and ChEMBL clients and prioritisation logic.

Provides typed clients for the Open Targets GraphQL API and ChEMBL REST API,
and logic for converting raw tractability/bioactivity evidence into the
DruggabilityEvidence schema consumed by ``oncocartograph.scoring``.

Note: as of this work package, only RNA-seq/copy-number (Ensembl ID) and
mutation (gene symbol, resolved via Open Targets ``mapIds``) candidates
are supported. Methylation candidates are CpG probe IDs, not gene
identifiers, and require a probe-to-gene mapping (e.g. the Illumina 450K
manifest) that is out of scope here -- see
``docs/adr/0008-druggability-evidence-sources.md``.
"""

from oncocartograph.drug_targets.chembl_client import ChEMBLClient, ChEMBLRequestError
from oncocartograph.drug_targets.druggability import build_druggability_evidence
from oncocartograph.drug_targets.open_targets_client import (
    OpenTargetsClient,
    OpenTargetsRequestError,
)
from oncocartograph.drug_targets.tractability import score_tractability

__all__ = [
    "ChEMBLClient",
    "ChEMBLRequestError",
    "OpenTargetsClient",
    "OpenTargetsRequestError",
    "build_druggability_evidence",
    "score_tractability",
]
