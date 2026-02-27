import logging
from typing import Dict, Any

logger = logging.getLogger(__name__)

class PSORTbClient:
    def __init__(self):
        pass

    def predict_localization(self, sequence: str, organism_type: str = "gram_negative") -> Dict[str, Any]:
        """
        Predict subcellular localization for a bacterial protein.

        Args:
            sequence: Amino acid sequence
            organism_type: 'gram_positive', 'gram_negative', or 'archaea'

        Returns:
            Dictionary with localization prediction and confidence
        """
        try:
            return self._predict_localization_rules(sequence, organism_type)
        except Exception as e:
            logger.error(f"PSORTb prediction failed: {e}")
            return {
                "localization": "unknown",
                "confidence": 0.0,
                "surface_score": 0.0,
                "scores": {}
            }

    def _predict_localization_rules(self, sequence: str, organism_type: str) -> Dict[str, Any]:
        """
        Rule-based localization prediction using signal peptides and transmembrane regions.
        """
        # Check for signal peptides (first 30 amino acids)
        n_terminal = sequence[:30] if len(sequence) > 30 else sequence

        # Signal peptide indicators
        has_signal = self._check_signal_peptide(n_terminal)

        # Check for transmembrane domains
        tm_regions = self._count_transmembrane_regions(sequence)

        # Check for specific localization signals
        has_lipoprotein_signal = self._check_lipoprotein_signal(sequence)

        # Basic prediction logic
        if has_lipoprotein_signal:
            localization = "outer_membrane"
            confidence = 0.8
        elif tm_regions > 2:
            localization = "inner_membrane"
            confidence = 0.7
        elif has_signal and tm_regions == 0:
            if organism_type == "gram_positive":
                localization = "extracellular"
            else:
                localization = "periplasmic"
            confidence = 0.6
        elif tm_regions == 1:
            localization = "inner_membrane"
            confidence = 0.5
        else:
            localization = "cytoplasmic"
            confidence = 0.7

        # Surface accessibility score (higher is better for vaccines)
        if localization in ["outer_membrane", "extracellular"]:
            surface_score = 1.0
        elif localization == "periplasmic":
            surface_score = 0.6
        elif localization == "inner_membrane":
            surface_score = 0.3
        else:  # cytoplasmic
            surface_score = 0.1

        return {
            "localization": localization,
            "confidence": confidence,
            "surface_score": surface_score,
            "scores": {
                "cytoplasmic": confidence if localization == "cytoplasmic" else 0.2,
                "inner_membrane": confidence if localization == "inner_membrane" else 0.1,
                "periplasmic": confidence if localization == "periplasmic" else 0.1,
                "outer_membrane": confidence if localization == "outer_membrane" else 0.1,
                "extracellular": confidence if localization == "extracellular" else 0.1
            },
            "features": {
                "signal_peptide": has_signal,
                "transmembrane_regions": tm_regions,
                "lipoprotein_signal": has_lipoprotein_signal
            }
        }

    def _check_signal_peptide(self, n_terminal: str) -> bool:
        """Check for signal peptide characteristics."""
        if len(n_terminal) < 15:
            return False

        # Signal peptides typically have:
        # - Positively charged N-region (first ~5 aa)
        # - Hydrophobic H-region (next ~10-15 aa)

        n_region = n_terminal[:5]
        h_region = n_terminal[5:15] if len(n_terminal) >= 15 else n_terminal[5:]

        # Check for positive charges in N-region
        positive_charges = sum(1 for aa in n_region if aa in 'KR')

        # Check for hydrophobic residues in H-region
        hydrophobic = sum(1 for aa in h_region if aa in 'AILMFWYV')
        hydrophobic_ratio = hydrophobic / len(h_region) if h_region else 0

        return positive_charges >= 1 and hydrophobic_ratio > 0.4

    def _count_transmembrane_regions(self, sequence: str) -> int:
        """Simple transmembrane domain prediction."""
        hydrophobic_aa = 'AILMFWYV'
        tm_count = 0

        # Scan sequence with sliding window
        window_size = 20
        hydrophobic_threshold = 0.65

        i = 0
        while i < len(sequence) - window_size + 1:
            window = sequence[i:i + window_size]
            hydrophobic_ratio = sum(1 for aa in window if aa in hydrophobic_aa) / window_size

            if hydrophobic_ratio >= hydrophobic_threshold:
                tm_count += 1
                # Skip ahead to avoid counting the same region multiple times
                i += window_size
            else:
                i += 1

        return tm_count

    def _check_lipoprotein_signal(self, sequence: str) -> bool:
        """Check for lipoprotein signal sequence pattern."""
        if len(sequence) < 20:
            return False

        n_terminal = sequence[:20]
        # Look for cysteine with appropriate preceding amino acids
        for i in range(1, len(n_terminal)):
            if n_terminal[i] == 'C':
                if i >= 2:
                    pattern = n_terminal[i-2:i+1]
                    if (pattern[0] in 'LVI' and
                        pattern[1] in 'ASTVI' and
                        pattern[2] == 'C'):
                        return True
        return False

# Global instance
psortb = PSORTbClient()