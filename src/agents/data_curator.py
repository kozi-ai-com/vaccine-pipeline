import logging
import uuid
from typing import List, Dict, Any
from src.models.candidate import CandidateProtein, PipelineRun, CandidateStatus, ConfidenceTier

logger = logging.getLogger(__name__)

class DataCuratorAgent:
    """
    Agent responsible for:
    1. Fetching protein sequences (UniProt, FASTA, or proteome search)
    2. Initial antigen screening (location, antigenicity, safety)
    3. Claude-based filtering decisions
    """

    def __init__(self):
        self.stage_name = "data_curation"
        # Import tools lazily to avoid circular imports
        self._uniprot = None
        self._vaxijen = None
        self._psortb = None
        self._claude = None
        self._db = None

    @property
    def uniprot(self):
        if self._uniprot is None:
            from src.tools.uniprot_client import uniprot
            self._uniprot = uniprot
        return self._uniprot

    @property
    def vaxijen(self):
        if self._vaxijen is None:
            from src.tools.vaxijen_client import vaxijen
            self._vaxijen = vaxijen
        return self._vaxijen

    @property
    def psortb(self):
        if self._psortb is None:
            from src.tools.psortb_client import psortb
            self._psortb = psortb
        return self._psortb

    @property
    def claude(self):
        if self._claude is None:
            from src.tools.claude_client import claude
            self._claude = claude
        return self._claude

    @property
    def db(self):
        if self._db is None:
            from src.storage.supabase_client import db
            self._db = db
        return self._db

    def run(self, pipeline_run: PipelineRun) -> PipelineRun:
        """
        Execute data curation stage.

        Args:
            pipeline_run: Current pipeline state

        Returns:
            Updated pipeline with curated candidates
        """
        logger.info(f"Starting data curation for {pipeline_run.pathogen_name}")

        try:
            # Step 1: Fetch protein sequences
            proteins = self._fetch_proteins(pipeline_run)
            logger.info(f"Fetched {len(proteins)} proteins")

            if not proteins:
                pipeline_run.add_error("data_curation", "No proteins found for input")
                return pipeline_run

            # Step 2: Convert to candidates and run screening
            candidates = []
            for i, protein_data in enumerate(proteins):
                logger.info(f"Processing protein {i+1}/{len(proteins)}: {protein_data.get('protein_name', 'Unknown')}")

                candidate = self._create_candidate(protein_data)

                # Run antigen screening tools
                self._run_antigen_screening(candidate)

                # Claude decision: advance or filter out
                decision = self._make_claude_decision(candidate)
                self._apply_decision(candidate, decision)

                candidates.append(candidate)

                # Save to database after each candidate
                try:
                    self.db.save_candidate(pipeline_run.run_id, candidate)
                except Exception as e:
                    logger.warning(f"Failed to save candidate to database: {e}")

            # Step 3: Update pipeline state
            pipeline_run.candidates = candidates
            pipeline_run.active_candidate_count = len([c for c in candidates if c.status == CandidateStatus.ACTIVE])
            pipeline_run.current_stage = "antigen_screening_complete"

            # Update run in database
            try:
                self.db.update_run_stage(pipeline_run.run_id, "antigen_screening_complete")
            except Exception as e:
                logger.warning(f"Failed to update run in database: {e}")

            logger.info(f"Data curation complete. {pipeline_run.active_candidate_count} active candidates")
            return pipeline_run

        except Exception as e:
            logger.error(f"Data curation failed: {e}")
            pipeline_run.add_error("data_curation", str(e))
            return pipeline_run

    def _fetch_proteins(self, pipeline_run: PipelineRun) -> List[Dict[str, Any]]:
        """Fetch protein sequences based on input type."""

        input_type = pipeline_run.input_type
        raw_input = pipeline_run.raw_input.strip()

        logger.info(f"Fetching proteins: type={input_type}, input={raw_input[:50]}...")

        if input_type == "uniprot_id":
            # Single protein by UniProt ID
            result = self.uniprot.fetch_protein_by_id(raw_input)
            return [result] if result else []

        elif input_type == "fasta":
            # Parse FASTA sequences
            return self.uniprot.parse_fasta_input(raw_input)

        elif input_type == "proteome":
            # Full proteome search (limited for MVP)
            return self.uniprot.search_proteome(raw_input, max_proteins=20)

        else:
            logger.error(f"Unknown input type: {input_type}")
            return []

    def _create_candidate(self, protein_data: Dict[str, Any]) -> CandidateProtein:
        """Convert protein data to candidate object."""

        return CandidateProtein(
            protein_id=protein_data["protein_id"],
            protein_name=protein_data["protein_name"],
            sequence=protein_data["sequence"],
            source=protein_data["source"],
            stage=self.stage_name,
            status=CandidateStatus.ACTIVE,
            confidence_tier=ConfidenceTier.UNSCORED,
            flags=[]
        )

    def _run_antigen_screening(self, candidate: CandidateProtein):
        """Run all antigen screening tools on a candidate."""

        logger.info(f"Screening candidate: {candidate.protein_name}")

        try:
            # 1. Protein localization (PSORTb)
            localization_result = self.psortb.predict_localization(candidate.sequence)
            candidate.psortb_localization = localization_result["localization"]

            # 2. Antigenicity prediction (VaxiJen)
            antigenicity_score = self.vaxijen.predict_antigenicity(candidate.sequence)
            candidate.vaxijen_score = antigenicity_score

            # 3. Basic safety flags
            self._check_basic_safety(candidate)

            logger.info(f"Screening complete for {candidate.protein_name}")
            logger.info(f"  Location: {candidate.psortb_localization}")
            logger.info(f"  Antigenicity: {candidate.vaxijen_score:.3f}")
            logger.info(f"  Flags: {candidate.flags}")

        except Exception as e:
            logger.error(f"Screening failed for {candidate.protein_name}: {e}")
            candidate.flags.append("screening_failed")

    def _check_basic_safety(self, candidate: CandidateProtein):
        """Basic safety checks before Claude decision."""

        sequence = candidate.sequence

        # Check sequence length
        if len(sequence) < 50:
            candidate.flags.append("too_short")
        elif len(sequence) > 2000:
            candidate.flags.append("very_long")

        # Check for unusual amino acids
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        if not set(sequence).issubset(valid_aa):
            candidate.flags.append("invalid_amino_acids")

        # Check for highly repetitive sequences
        if self._is_highly_repetitive(sequence):
            candidate.flags.append("repetitive_sequence")

        # Check for potential signal peptide
        if self._has_potential_signal_peptide(sequence):
            candidate.flags.append("signal_peptide")

    def _is_highly_repetitive(self, sequence: str) -> bool:
        """Check if sequence has high repetitive content."""
        if len(sequence) < 20:
            return False

        # Simple check: look for runs of the same amino acid
        max_run = 0
        current_run = 1

        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
            else:
                max_run = max(max_run, current_run)
                current_run = 1
        max_run = max(max_run, current_run)

        return max_run > 10 or max_run > len(sequence) * 0.3

    def _has_potential_signal_peptide(self, sequence: str) -> bool:
        """Rough check for signal peptide."""
        if len(sequence) < 25:
            return False

        n_terminal = sequence[:25]
        hydrophobic = sum(1 for aa in n_terminal if aa in "AILMFWYV")
        return hydrophobic / len(n_terminal) > 0.5

    def _make_claude_decision(self, candidate: CandidateProtein) -> Dict[str, Any]:
        """Get Claude's decision on whether to advance this candidate."""

        # Prepare data for Claude
        protein_data = {
            "protein_name": candidate.protein_name,
            "protein_id": candidate.protein_id,
            "sequence_length": len(candidate.sequence),
            "psortb_localization": candidate.psortb_localization,
            "vaxijen_score": candidate.vaxijen_score,
            "flags": candidate.flags,
            "source": candidate.source
        }

        # Get Claude's decision
        return self.claude.make_antigen_decision(protein_data)

    def _apply_decision(self, candidate: CandidateProtein, decision: Dict[str, Any]):
        """Apply Claude's decision to the candidate."""

        # Record the decision
        candidate.add_decision(
            stage="antigen_screening",
            decision=decision["decision"],
            reasoning=decision["reasoning"]
        )

        # Update candidate status
        if decision["decision"] == "advance":
            candidate.status = CandidateStatus.ACTIVE
            candidate.confidence_tier = ConfidenceTier(decision.get("confidence", "medium"))
        elif decision["decision"] == "deprioritize":
            candidate.status = CandidateStatus.DEPRIORITIZED
            candidate.confidence_tier = ConfidenceTier.LOW
        elif decision["decision"] == "discard":
            candidate.status = CandidateStatus.DISCARDED
            candidate.confidence_tier = ConfidenceTier.UNCERTAIN

        # Add any new flags
        if decision.get("flags"):
            candidate.flags.extend(decision["flags"])

# Global instance
data_curator = DataCuratorAgent()