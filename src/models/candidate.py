from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime
from enum import Enum

class ConfidenceTier(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNCERTAIN = "uncertain"
    UNSCORED = "unscored"

class CandidateStatus(str, Enum):
    ACTIVE = "active"
    DEPRIORITIZED = "deprioritized"
    DISCARDED = "discarded"

class EpitopeType(str, Enum):
    CTL = "CTL"  # Cytotoxic T-lymphocyte (CD8+ T-cell)
    HTL = "HTL"  # Helper T-lymphocyte (CD4+ T-cell)
    B_CELL_LINEAR = "B-cell-linear"
    B_CELL_CONFORMATIONAL = "B-cell-conformational"

class EpitopeResult(BaseModel):
    """Individual epitope prediction with all associated scores."""
    sequence: str = Field(description="Amino acid sequence of the epitope")
    epitope_type: EpitopeType
    hla_allele: Optional[str] = Field(None, description="HLA allele if applicable (for T-cell epitopes)")

    # Binding predictions
    ic50_nm: Optional[float] = Field(None, description="IC50 binding affinity in nanomolar")
    percentile_rank: Optional[float] = Field(None, description="Percentile rank (lower is better)")

    # Biological properties
    conservation_score: Optional[float] = Field(None, description="Sequence conservation across strains")
    allergenicity_safe: Optional[bool] = Field(None, description="Passed allergenicity screening")
    toxicity_safe: Optional[bool] = Field(None, description="Passed toxicity screening")

    # Meta
    confidence_tier: ConfidenceTier = ConfidenceTier.UNSCORED
    tool_outputs: Dict[str, Any] = Field(default_factory=dict, description="Raw tool outputs")
    created_at: datetime = Field(default_factory=datetime.now)

class CandidateProtein(BaseModel):
    """A protein candidate moving through the vaccine pipeline."""

    # Basic identification
    protein_id: str = Field(description="UniProt ID or other identifier")
    protein_name: str = Field(description="Human-readable protein name")
    sequence: str = Field(description="Full amino acid sequence")
    source: str = Field(description="Where this protein came from (uniprot, ncbi, user_input)")

    # Pipeline tracking
    stage: str = Field(description="Current pipeline stage")
    status: CandidateStatus = CandidateStatus.ACTIVE
    confidence_tier: ConfidenceTier = ConfidenceTier.UNSCORED
    flags: List[str] = Field(default_factory=list, description="Warning flags or special notes")

    # Antigen screening scores (filled by early pipeline stages)
    psortb_localization: Optional[str] = Field(None, description="Cellular localization prediction")
    tmhmm_helices: Optional[int] = Field(None, description="Number of transmembrane helices")
    vaxijen_score: Optional[float] = Field(None, description="Antigenicity score")
    blast_human_identity: Optional[float] = Field(None, description="Max identity to human proteins")

    # Structure information
    structure_source: Optional[str] = Field(None, description="alphafold_db, colabfold, or none")
    structure_pdb_path: Optional[str] = Field(None, description="Local path to PDB file")

    # Epitope predictions (filled by prediction stages)
    ctl_epitopes: List[EpitopeResult] = Field(default_factory=list)
    htl_epitopes: List[EpitopeResult] = Field(default_factory=list)
    bcell_epitopes: List[EpitopeResult] = Field(default_factory=list)

    # Conservation analysis
    conservation_profile: Optional[Dict[str, Any]] = Field(None, description="Strain conservation data")

    # Coverage analysis
    hla_coverage_global: Optional[float] = Field(None, description="Global population HLA coverage")
    hla_coverage_africa: Optional[float] = Field(None, description="African population HLA coverage")

    # Decision audit trail - FIXED
    decisions: List[Dict[str, Any]] = []

    # Metadata
    created_at: datetime = Field(default_factory=datetime.now)
    updated_at: datetime = Field(default_factory=datetime.now)

    def add_decision(self, stage: str, decision: str, reasoning: str, **kwargs):
        """Add a decision to the audit trail."""
        self.decisions.append({
            "stage": stage,
            "decision": decision,
            "reasoning": reasoning,
            "timestamp": datetime.now().isoformat(),
            **kwargs
        })
        self.updated_at = datetime.now()

    def get_total_epitopes(self) -> int:
        """Count all predicted epitopes."""
        return len(self.ctl_epitopes) + len(self.htl_epitopes) + len(self.bcell_epitopes)

    def get_high_confidence_epitopes(self) -> int:
        """Count high-confidence epitopes only."""
        all_epitopes = self.ctl_epitopes + self.htl_epitopes + self.bcell_epitopes
        return len([e for e in all_epitopes if e.confidence_tier == ConfidenceTier.HIGH])

class PipelineRun(BaseModel):
    """Represents a complete pipeline execution."""
    run_id: str
    pathogen_name: str
    input_type: str  # 'uniprot_id', 'fasta', 'proteome'
    raw_input: str

    # User constraints
    target_populations: List[str] = Field(default_factory=lambda: ["global"])
    coverage_threshold: float = 0.70
    max_candidates_output: int = 20

    # Pipeline state
    current_stage: str = "initialization"
    candidates: List[CandidateProtein] = Field(default_factory=list)
    active_candidate_count: int = 0

    # Coverage loop control
    coverage_loop_count: int = 0
    coverage_met: bool = False

    # Final outputs
    construct_sequence: Optional[str] = None
    construct_properties: Optional[Dict[str, Any]] = None
    experiment_plan: Optional[str] = None
    final_report: Optional[Dict[str, Any]] = None

    # Error handling - FIXED
    errors: List[Dict[str, Any]] = []
    warnings: List[str] = []

    # Metadata
    created_at: datetime = Field(default_factory=datetime.now)
    completed_at: Optional[datetime] = None
    status: str = "running"  # running, completed, failed

    def add_error(self, stage: str, error: str, **kwargs):
        """Add an error to the error log."""
        self.errors.append({
            "stage": stage,
            "error": error,
            "timestamp": datetime.now().isoformat(),
            **kwargs
        })

    def add_warning(self, warning: str):
        """Add a warning message."""
        self.warnings.append(f"{datetime.now().isoformat()}: {warning}")

    def get_active_candidates(self) -> List[CandidateProtein]:
        """Get only active candidates."""
        return [c for c in self.candidates if c.status == CandidateStatus.ACTIVE]