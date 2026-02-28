import requests
import time
import logging
from typing import List, Optional, Dict, Any
from Bio import SeqIO
from io import StringIO

logger = logging.getLogger(__name__)


class UniProtClient:
    def __init__(self):
        self.base_url = "https://rest.uniprot.org"
        self.session = requests.Session()
        # UniProt doesn't require API key but appreciates User-Agent
        self.session.headers.update({
            "User-Agent": "Kozi-Vaccine-Pipeline/1.0 (ask@kozi-ai.com)"
        })

    def fetch_protein_by_id(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch a single protein by UniProt ID.
        Returns standardized protein data or None if not found.
        """
        try:
            url = f"{self.base_url}/uniprotkb/{uniprot_id}.json"

            logger.info(f"Fetching protein {uniprot_id} from UniProt")
            response = self.session.get(url, timeout=30)

            if response.status_code == 404:
                logger.warning(f"Protein {uniprot_id} not found in UniProt")
                return None

            response.raise_for_status()
            data = response.json()

            # Extract key information
            protein_info = {
                "protein_id": uniprot_id,
                "protein_name": self._extract_protein_name(data),
                "sequence": data.get("sequence", {}).get("value", ""),
                "organism": self._extract_organism(data),
                "length": data.get("sequence", {}).get("length", 0),
                "keywords": self._extract_keywords(data),
                "subcellular_location": self._extract_location(data),
                "source": "uniprot"
            }

            logger.info(
                f"Successfully fetched {uniprot_id}: {protein_info['protein_name']}")
            return protein_info

        except requests.RequestException as e:
            logger.error(f"Failed to fetch protein {uniprot_id}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error parsing protein {uniprot_id}: {e}")
            return None

    def parse_fasta_input(self, fasta_text: str) -> List[Dict[str, Any]]:
        """
        Parse FASTA format text into protein data.
        Handles both single and multi-sequence FASTA.
        """
        try:
            proteins = []
            fasta_io = StringIO(fasta_text)

            for record in SeqIO.parse(fasta_io, "fasta"):
                protein_info = {
                    "protein_id": record.id,
                    "protein_name": record.description,
                    "sequence": str(record.seq),
                    "organism": "User provided",
                    "length": len(record.seq),
                    "keywords": [],
                    "subcellular_location": "Unknown",
                    "source": "user_input"
                }
                proteins.append(protein_info)

            logger.info(f"Parsed {len(proteins)} proteins from FASTA input")
            return proteins

        except Exception as e:
            logger.error(f"Failed to parse FASTA input: {e}")
            return []

    def _extract_protein_name(self, data: Dict) -> str:
        """Extract the main protein name from UniProt JSON."""
        try:
            protein_description = data.get("proteinDescription", {})
            recommended_name = protein_description.get("recommendedName", {})

            if "fullName" in recommended_name:
                return recommended_name["fullName"]["value"]

            # Fallback to alternative names
            alternative_names = protein_description.get("alternativeNames", [])
            if alternative_names:
                return alternative_names[0]["fullName"]["value"]

            # Last resort - use the entry name
            return data.get("uniProtkbId", "Unknown protein")

        except Exception:
            return "Unknown protein"

    def _extract_organism(self, data: Dict) -> str:
        """Extract organism name from UniProt JSON."""
        try:
            organism = data.get("organism", {})
            return organism.get("scientificName", "Unknown organism")
        except Exception:
            return "Unknown organism"

    def _extract_keywords(self, data: Dict) -> List[str]:
        """Extract keywords from UniProt JSON."""
        try:
            keywords = data.get("keywords", [])
            return [kw.get("name", "") for kw in keywords if kw.get("name")]
        except Exception:
            return []

    def _extract_location(self, data: Dict) -> str:
        """Extract subcellular location from UniProt JSON."""
        try:
            comments = data.get("comments", [])
            for comment in comments:
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    locations = comment.get("subcellularLocations", [])
                    if locations:
                        location = locations[0].get("location", {})
                        return location.get("value", "Unknown")
            return "Unknown"
        except Exception:
            return "Unknown"


# Global instance
uniprot = UniProtClient()
