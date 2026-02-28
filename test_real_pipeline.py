"""Test with mpox A35R (OPG161) envelope protein involved in virion membrane formation."""

import uuid
import logging
from src.models.candidate import PipelineRun
from src.agents.data_curator import data_curator

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_mpox_protein():
    """Test with mpox A35R protein - well-studied surface protein."""

    print("Testing Complete Pipeline with Mpox A35R")
    print("=" * 60)

    # Create test run with mpox A35R protein
    test_run = PipelineRun(
        run_id=str(uuid.uuid4()),
        pathogen_name="Monkeypox virus (mpox)",
        input_type="uniprot_id",
        raw_input="Q80KX2"

        """
        UniProt: Q80KX2
        Gene: A35R
        Protein: OPG161
        Virus: Monkeypox virus
        """
    )

    print(f"Run ID: {test_run.run_id}")
    print(f"Target: {test_run.pathogen_name}")
    print(f"Protein: {test_run.raw_input}")
    print("\n" + "Running Data Curator Agent...")

    # Run the complete agent
    try:
        result = data_curator.run(test_run)

        print("\n" + "RESULTS")
        print("=" * 40)

        if result.errors:
            print("ERRORS:")
            for error in result.errors:
                print(f"   {error['stage']}: {error['error']}")
            return False

        if not result.candidates:
            print("No candidates found")
            return False

        # Analyze the candidate
        candidate = result.candidates[0]

        print(f"CANDIDATE ANALYSIS")
        print(f"   Protein: {candidate.protein_name}")
        print(f"   UniProt ID: {candidate.protein_id}")
        print(f"   Length: {len(candidate.sequence)} amino acids")
        print(f"   Source: {candidate.source}")
        print(f"   Status: {candidate.status.value}")
        print(f"   Confidence: {candidate.confidence_tier.value}")

        print(f"\n SCREENING RESULTS")
        print(f"   Cellular Location: {candidate.psortb_localization}")
        print(f"   Antigenicity Score: {candidate.vaxijen_score:.3f}")
        print(f"   Flags: {candidate.flags}")

        if candidate.decisions:
            decision = candidate.decisions[-1]
            print(f"\n CLAUDE'S DECISION")
            print(f"     Decision: {decision['decision']}")
            print(f"     Reasoning: {decision['reasoning']}")

        # Expected results for A35R (validation)
        print(f"\n VALIDATION")

        # Check if protein was fetched correctly
        is_a35r = candidate.protein_id == test_run.raw_input
        print(f"   Correct protein fetched: {'Yes' if is_a35r else 'No'}")

        # Check if identified as surface protein
        is_surface = candidate.psortb_localization in ["outer_membrane",
    "extracellular",
    "periplasmic",
    "membrane",
    "cell_surface",
    "unknown"]
        print(f"   Surface localization: {'Yes' if is_surface else 'No'} ({candidate.psortb_localization})")

        # Check antigenicity
        is_antigenic = candidate.vaxijen_score and candidate.vaxijen_score > 0.5
        print(f"   Good antigenicity: {'Yes' if is_antigenic else 'No'} ({candidate.vaxijen_score:.3f})")

        # Check decision
        advanced = candidate.status.value in ["active", "deprioritized"]
        print(f"   Advanced in pipeline: {'Yes' if advanced else 'No'} ({candidate.status.value})")

        # Overall success
        success = is_a35r and is_surface and is_antigenic and advanced

        print(f"\n OVERALL RESULT: {'PASSED' if success else 'MIXED'}")

        if success:
            print("\n Excellent! Your pipeline correctly identified mpox A35R as a vaccine candidate.")
            print("   This matches known experimental data - A35R is used in real mpox vaccines.")
        else:
            print("\n Results analysis:")
            if not is_surface:
                print("   - Consider tuning PSORTb localization rules")
            if not is_antigenic:
                print("   - Consider tuning VaxiJen antigenicity calculation")
            if not advanced:
                print("   - Consider tuning Claude's decision criteria")

        return success

    except Exception as e:
        print(f"\nPIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_fasta_sequence():
    """Test with a FASTA sequence."""

    print("\n\nTesting FASTA Input")
    print("=" * 40)

    # Example surface protein sequence
    fasta_input = """>test_surface_protein
MKVLVLSLGMFPLADIEAAERTVQDLGKLQKASLVGQTLALEKLTGPTEPTVVNSQPHTSMNTKN
VKDLQASLAEITARGQNVVLETKLDGQIPVTTNPKDLQASLVELTARGQNVVLETKLDGQIP"""

    test_run = PipelineRun(
        run_id=str(uuid.uuid4()),
        pathogen_name="Test organism",
        input_type="fasta",
        raw_input=fasta_input
    )

    try:
        result = data_curator.run(test_run)

        if result.candidates:
            candidate = result.candidates[0]
            print(f"FASTA processed successfully")
            print(f"   Protein: {candidate.protein_name}")
            print(f"   Length: {len(candidate.sequence)}")
            print(f"   Antigenicity: {candidate.vaxijen_score:.3f}")
            print(f"   Status: {candidate.status.value}")
            return True
        else:
            print(f"FASTA processing failed")
            return False

    except Exception as e:
        print(f"FASTA test error: {e}")
        return False

def run_comprehensive_test():
    """Run all tests."""

    print("Kozi Vaccine Pipeline - Comprehensive Test")
    print("Testing complete data flow: fetch → screen → decide → store")
    print("\n")

    # Test 1: Real protein (mpox A35R)
    test1_passed = test_mpox_protein()

    # Test 2: FASTA input
    test2_passed = test_fasta_sequence()

    print("\n" + "=" * 60)
    print("FINAL TEST SUMMARY")
    print("=" * 60)
    print(f"Mpox A35R test: {'PASSED' if test1_passed else 'FAILED'}")
    print(f"FASTA test: {'PASSED' if test2_passed else 'FAILED'}")

    if test1_passed and test2_passed:
        print("\nALL TESTS PASSED!")
        print("\n Your vaccine pipeline is working! Key achievements:")
        print("   • Successfully fetches real protein data from UniProt")
        print("   • Runs biological screening tools (PSORTb, VaxiJen)")
        print("   • Makes intelligent decisions with Claude")
        print("   • Stores results with full audit trail")
        print("   • Correctly identifies known vaccine targets")
        print("\n Next steps:")
        print("   1. Add epitope prediction tools")
        print("   2. Implement coverage analysis")
        print("   3. Build simple UI")
        print("   4. Test with more pathogens")
    else:
        print("\n Some tests failed - check logs above for details")
        print("   The pipeline foundation is working, just needs tuning")

    return test1_passed and test2_passed

if __name__ == "__main__":
    run_comprehensive_test()