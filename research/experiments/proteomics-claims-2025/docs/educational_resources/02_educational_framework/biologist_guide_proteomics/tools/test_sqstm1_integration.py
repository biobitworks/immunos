#!/usr/bin/env python3
"""
Test script to demonstrate SQSTM1-enhanced UPS analysis
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from tools.uniprot_analysis import UPSAnalyzer
import pandas as pd

def test_sqstm1_ups_integration():
    """Test the enhanced UPS analyzer with SQSTM1-specific capabilities"""

    print("🧪 Testing SQSTM1-Enhanced UPS Analysis")
    print("=" * 50)

    # Initialize the enhanced UPS analyzer
    analyzer = UPSAnalyzer()

    # Test proteins including SQSTM1 and related proteins
    test_proteins = [
        'SQSTM1',  # p62 - should be highly classified
        'UBQLN1',  # Ubiquilin - UPS shuttle factor
        'NBR1',    # Neighbor of BRCA1 - autophagy receptor
        'OPTN',    # Optineurin - autophagy receptor
        'LC3B',    # LC3 - autophagy marker
        'KEAP1',   # Kelch-like ECH-associated protein 1
        'NFE2L2',  # NRF2 - transcription factor
        'MTOR',    # mTOR - nutrient sensing
        'PSMB5',   # Proteasome subunit
        'UBE2D1',  # E2 ubiquitin conjugating enzyme
        'HUWE1',   # E3 ubiquitin ligase
        'USP7',    # Deubiquitinase
        'TRAF6',   # TNF receptor associated factor 6
        'ATG5',    # Autophagy related 5
        'BECN1',   # Beclin 1
        'GAPDH'    # Control - non-UPS protein
    ]

    print(f"Analyzing {len(test_proteins)} test proteins...")
    print()

    results = []

    for protein in test_proteins:
        print(f"🔍 Analyzing {protein}...")

        # Classify the protein
        classification = analyzer.classify_ups_component(protein)

        # Calculate UPS score based on evidence
        ups_score = len(classification['evidence']) / 10.0  # Normalize evidence count

        # Combine results
        result = {
            'protein': protein,
            'is_ups_component': classification['is_ups_component'],
            'ups_category': classification['ups_category'],
            'ups_subcategory': classification['ups_subcategory'],
            'confidence': classification['confidence'],
            'evidence_count': len(classification['evidence']),
            'ups_score': ups_score
        }

        results.append(result)

        # Print detailed results for SQSTM1
        if protein == 'SQSTM1':
            print(f"  📊 SQSTM1 Detailed Analysis:")
            print(f"     UPS Component: {classification['is_ups_component']}")
            print(f"     Category: {classification['ups_category']}")
            print(f"     Subcategory: {classification['ups_subcategory']}")
            print(f"     Confidence: {classification['confidence']}")
            print(f"     Evidence items: {len(classification['evidence'])}")
            print(f"     UPS score: {ups_score:.2f}")

            # Show first few evidence items
            if classification['evidence']:
                print(f"     Sample evidence:")
                for i, evidence in enumerate(classification['evidence'][:3]):
                    print(f"       {i+1}. {evidence}")
                if len(classification['evidence']) > 3:
                    print(f"       ... and {len(classification['evidence'])-3} more")

        print(f"     ✅ UPS: {classification['is_ups_component']} | "
              f"Confidence: {classification['confidence']} | "
              f"Score: {ups_score:.2f}")
        print()

    # Create summary table
    df = pd.DataFrame(results)
    df = df.sort_values(['is_ups_component', 'ups_score'], ascending=[False, False])

    print("📈 Summary Results:")
    print("=" * 80)
    print(df.to_string(index=False))

    # SQSTM1-specific analysis
    sqstm1_result = df[df['protein'] == 'SQSTM1'].iloc[0]

    print(f"\n🎯 SQSTM1 Performance:")
    print(f"   UPS Classification: {'✅ PASS' if sqstm1_result['is_ups_component'] else '❌ FAIL'}")
    print(f"   Confidence Level: {sqstm1_result['confidence']}")
    print(f"   UPS Score: {sqstm1_result['ups_score']:.2f}")
    print(f"   Evidence Items: {sqstm1_result['evidence_count']}")

    # Count UPS classifications
    ups_count = df['is_ups_component'].sum()
    total_count = len(df)

    print(f"\n📊 Overall Performance:")
    print(f"   UPS Proteins Identified: {ups_count}/{total_count}")
    print(f"   Classification Rate: {ups_count/total_count*100:.1f}%")

    # Show high-confidence classifications
    high_conf = df[df['confidence'] == 'high']
    print(f"   High Confidence: {len(high_conf)} proteins")
    if not high_conf.empty:
        print(f"   High-conf proteins: {', '.join(high_conf['protein'].tolist())}")

    print(f"\n✨ SQSTM1 integration test completed!")
    return df

if __name__ == "__main__":
    test_sqstm1_ups_integration()