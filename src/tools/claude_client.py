import os
import json
import logging
from typing import Dict, Any
from dotenv import load_dotenv

# Load .env file explicitly
load_dotenv()

logger = logging.getLogger(__name__)

class ClaudeClient:
    def __init__(self):
        # Load environment again just to be sure
        load_dotenv(override=True)

        self.api_key = os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            # Try to load from .env file one more time
            env_path = os.path.join(os.getcwd(), ".env")
            if os.path.exists(env_path):
                load_dotenv(env_path, override=True)
                self.api_key = os.getenv("ANTHROPIC_API_KEY")

        if not self.api_key:
            print("ANTHROPIC_API_KEY not found!")
            print("Current working directory:", os.getcwd())
            print("Looking for .env file at:", os.path.join(os.getcwd(), ".env"))
            print("File exists:", os.path.exists(os.path.join(os.getcwd(), ".env")))
            raise ValueError("ANTHROPIC_API_KEY must be set in .env file")

        self.client = None
        self.model = "claude-opus-4-1"

        # Initialize Anthropic client
        try:
            import anthropic
            self.client = anthropic.Anthropic(api_key=self.api_key)
            logger.info("Claude API client initialized successfully")
            print(f"Claude client initialized with key: {self.api_key[:20]}...")
        except ImportError as e:
            logger.error("anthropic package not installed. Run: pip install anthropic")
            raise ImportError("anthropic package not available") from e
        except Exception as e:
            logger.error(f"Failed to initialize Claude client: {e}")
            raise

    def make_antigen_decision(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:

        if not self.client:
            logger.error("Claude client not initialized")
            return self._get_fallback_decision(protein_data)

        prompt = f"""You are evaluating a protein for vaccine candidate screening.

Protein Information:
- Name: {protein_data.get('protein_name', 'Unknown')}
- Length: {protein_data.get('sequence_length', 'Unknown')} amino acids
- Cellular location: {protein_data.get('psortb_localization', 'Unknown')}
- Antigenicity score: {protein_data.get('vaxijen_score', 'Not calculated')} (>0.5 is good)
- Source: {protein_data.get('source', 'Unknown')}
- Flags: {protein_data.get('flags', [])}

Decision criteria:
- Surface proteins (extracellular, outer membrane) are preferred
- High antigenicity scores (>0.5) indicate good immune response
- Reasonable length (50-2000 amino acids) for vaccine development
- No major safety flags

Respond in JSON only:
{{
    "decision": "advance" or "deprioritize" or "discard",
    "reasoning": "one sentence explanation",
    "confidence": "high" or "medium" or "low",
    "flags": ["flag1", "flag2"] or []
}}"""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=256,
                messages=[{"role": "user", "content": prompt}]
            )

            decision_text = response.content[0].text.strip()

            # Clean up the response (remove any markdown formatting)
            if "```json" in decision_text:
                decision_text = decision_text.split("```json")[1].split("```")[0].strip()
            elif "```" in decision_text:
                # Handle cases where there's ``` but no json
                parts = decision_text.split("```")
                for part in parts:
                    if "{" in part and "}" in part:
                        decision_text = part.strip()
                        break

            decision = json.loads(decision_text)

            logger.info(f"Claude decision for {protein_data.get('protein_name')}: {decision['decision']}")
            return decision

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse Claude response as JSON: {e}")
            logger.error(f"Raw response: {decision_text}")
            return self._get_fallback_decision(protein_data)
        except Exception as e:
            logger.error(f"Claude API call failed: {e}")
            return self._get_fallback_decision(protein_data)

    def _get_fallback_decision(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Fallback decision logic when Claude API fails.
        Uses rule-based approach based on screening results.
        """
        logger.info("Using fallback decision logic")
        print("Using rule-based decision (Claude unavailable)")

        # Extract key metrics
        antigenicity = protein_data.get('vaxijen_score', 0)
        location = protein_data.get('psortb_localization', '').lower()
        flags = protein_data.get('flags', [])
        sequence_length = protein_data.get('sequence_length', 0)

        # Rule-based decision
        score = 0

        # Antigenicity scoring
        if antigenicity > 0.7:
            score += 3
        elif antigenicity > 0.5:
            score += 2
        elif antigenicity > 0.3:
            score += 1

        # Location scoring
        if location in ['outer_membrane', 'extracellular']:
            score += 3
        elif location in ['periplasmic']:
            score += 2
        elif location in ['inner_membrane']:
            score += 1

        # Length scoring
        if 100 <= sequence_length <= 1000:
            score += 1
        elif sequence_length < 50 or sequence_length > 2000:
            score -= 1

        # Flag penalties
        if 'too_short' in flags or 'very_long' in flags:
            score -= 2
        if 'screening_failed' in flags:
            score -= 3

        # Make decision
        if score >= 5:
            decision = "advance"
            confidence = "high"
            reasoning = f"High antigenicity ({antigenicity:.2f}) and good localization ({location})"
        elif score >= 3:
            decision = "advance"
            confidence = "medium"
            reasoning = f"Moderate scores: antigenicity {antigenicity:.2f}, location {location}"
        elif score >= 1:
            decision = "deprioritize"
            confidence = "low"
            reasoning = f"Low scores but keeping for potential: antigenicity {antigenicity:.2f}"
        else:
            decision = "discard"
            confidence = "low"
            reasoning = f"Poor vaccine candidate: low antigenicity ({antigenicity:.2f})"

        return {
            "decision": decision,
            "reasoning": reasoning,
            "confidence": confidence,
            "flags": ["fallback_decision"]
        }

    def test_api(self) -> bool:
        """Test Claude API connection."""
        if not self.client:
            return False

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=10,
                messages=[{"role": "user", "content": "Respond with just: API working"}]
            )

            result = "api working" in response.content[0].text.lower()
            logger.info(f"Claude API test: {'Success' if result else 'Failed'}")
            return result

        except Exception as e:
            logger.error(f"Claude API test failed: {e}")
            return False

# Global instance
claude = ClaudeClient()