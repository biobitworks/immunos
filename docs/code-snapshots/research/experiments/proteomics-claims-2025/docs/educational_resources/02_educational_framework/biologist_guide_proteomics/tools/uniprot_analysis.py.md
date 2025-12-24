---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/uniprot_analysis.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/uniprot_analysis.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
UniProt Database Integration and Analysis Tools

This module provides comprehensive tools for:
1. Accessing UniProt database programmatically
2. Validating protein annotations against UniProtKB
3. Handling multiple UniProt IDs per gene
4. Identifying UPS (Ubiquitin-Proteasome System) interactions
5. Cross-referencing with external databases

Author: Proteomics Analysis Framework
"""

import requests
import pandas as pd
import numpy as np
import time
import re
from typing import Dict, List, Tuple, Optional, Union
import json
from urllib.parse import urlencode
import warnings


class UniProtAPI:
    """
    Interface for UniProt REST API with rate limiting and error handling
    """

    def __init__(self, rate_limit: float = 1.0):
        """
        Initialize UniProt API client

        Parameters:
        -----------
        rate_limit : float
            Minimum time between requests (seconds)
        """
        self.base_url = "https://rest.uniprot.org"
        self.rate_limit = rate_limit
        self.last_request_time = 0
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ProteomicsAnalysisFramework/1.0'
        })

    def _rate_limit_wait(self):
        """Ensure proper rate limiting between requests"""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        if time_since_last < self.rate_limit:
            time.sleep(self.rate_limit - time_since_last)
        self.last_request_time = time.time()

    def _make_request(self, endpoint: str, params: Dict = None) -> Optional[Dict]:
        """
        Make rate-limited request to UniProt API

        Parameters:
        -----------
        endpoint : str
            API endpoint
        params : dict
            Query parameters

        Returns:
        --------
        dict or None
            JSON response or None if failed
        """
        self._rate_limit_wait()

        try:
            url = f"{self.base_url}/{endpoint}"
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()

            if 'application/json' in response.headers.get('content-type', ''):
                return response.json()
            else:
                return {'text': response.text}

        except requests.exceptions.RequestException as e:
            print(f"API request failed: {e}")
            return None

    def validate_uniprot_id(self, uniprot_id: str) -> Dict:
        """
        Validate UniProt ID against UniProtKB

        Parameters:
        -----------
        uniprot_id : str
            UniProt accession to validate

        Returns:
        --------
        dict
            Validation results
        """
        # Basic format validation
        uniprot_pattern = r'^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]|[0-9]{5})$'

        result = {
            'uniprot_id': uniprot_id,
            'format_valid': bool(re.match(uniprot_pattern, uniprot_id)),
            'exists': False,
            'is_reviewed': False,
            'gene_name': None,
            'protein_name': None,
            'organism': None,
            'error': None
        }

        if not result['format_valid']:
            result['error'] = 'Invalid UniProt ID format'
            return result

        # Check existence in UniProtKB
        endpoint = f"uniprotkb/{uniprot_id}"
        params = {'format': 'json'}

        response = self._make_request(endpoint, params)

        if response:
            result['exists'] = True
            try:
                # Extract key information
                entry = response
                result['is_reviewed'] = entry.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)'

                # Gene name
                if 'genes' in entry and entry['genes']:
                    result['gene_name'] = entry['genes'][0].get('geneName', {}).get('value')

                # Protein name
                if 'proteinDescription' in entry:
                    result['protein_name'] = entry['proteinDescription'].get('recommendedName', {}).get('fullName', {}).get('value')

                # Organism
                if 'organism' in entry:
                    result['organism'] = entry['organism'].get('scientificName')

            except (KeyError, TypeError) as e:
                result['error'] = f'Error parsing UniProt response: {e}'
        else:
            result['error'] = 'UniProt ID not found or API error'

        return result

    def get_protein_info(self, uniprot_id: str, include_features: bool = True) -> Dict:
        """
        Retrieve comprehensive protein information from UniProt

        Parameters:
        -----------
        uniprot_id : str
            UniProt accession
        include_features : bool
            Whether to include sequence features

        Returns:
        --------
        dict
            Comprehensive protein information
        """
        endpoint = f"uniprotkb/{uniprot_id}"
        params = {'format': 'json'}

        response = self._make_request(endpoint, params)

        if not response:
            return {'error': 'Failed to retrieve protein information'}

        try:
            entry = response
            info = {
                'uniprot_id': uniprot_id,
                'gene_name': None,
                'protein_name': None,
                'alternative_names': [],
                'organism': None,
                'taxonomy_id': None,
                'sequence_length': None,
                'molecular_weight': None,
                'subcellular_location': [],
                'go_terms': [],
                'keywords': [],
                'pathways': [],
                'interactions': [],
                'features': [],
                'is_reviewed': False
            }

            # Basic information
            info['is_reviewed'] = entry.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)'

            # Gene information
            if 'genes' in entry and entry['genes']:
                gene = entry['genes'][0]
                info['gene_name'] = gene.get('geneName', {}).get('value')

            # Protein names
            if 'proteinDescription' in entry:
                desc = entry['proteinDescription']
                if 'recommendedName' in desc:
                    info['protein_name'] = desc['recommendedName'].get('fullName', {}).get('value')

                if 'alternativeNames' in desc:
                    for alt_name in desc['alternativeNames']:
                        if 'fullName' in alt_name:
                            info['alternative_names'].append(alt_name['fullName'].get('value'))

            # Organism
            if 'organism' in entry:
                org = entry['organism']
                info['organism'] = org.get('scientificName')
                info['taxonomy_id'] = org.get('taxonId')

            # Sequence information
            if 'sequence' in entry:
                seq = entry['sequence']
                info['sequence_length'] = seq.get('length')
                info['molecular_weight'] = seq.get('molWeight')

            # Subcellular location
            if 'comments' in entry:
                for comment in entry['comments']:
                    if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                        for location in comment.get('subcellularLocations', []):
                            if 'location' in location:
                                info['subcellular_location'].append(location['location'].get('value'))

            # GO terms
            if 'uniProtKBCrossReferences' in entry:
                for xref in entry['uniProtKBCrossReferences']:
                    if xref.get('database') == 'GO':
                        go_id = xref.get('id')
                        go_term = None
                        if 'properties' in xref:
                            for prop in xref['properties']:
                                if prop.get('key') == 'GoTerm':
                                    go_term = prop.get('value')
                                    break
                        info['go_terms'].append({'id': go_id, 'term': go_term})

            # Keywords
            if 'keywords' in entry:
                for keyword in entry['keywords']:
                    info['keywords'].append(keyword.get('value'))

            # Features (if requested)
            if include_features and 'features' in entry:
                for feature in entry['features']:
                    info['features'].append({
                        'type': feature.get('type'),
                        'description': feature.get('description'),
                        'location': feature.get('location')
                    })

            return info

        except (KeyError, TypeError) as e:
            return {'error': f'Error parsing protein information: {e}'}

    def batch_query(self, uniprot_ids: List[str], chunk_size: int = 25) -> Dict:
        """
        Query multiple UniProt IDs efficiently

        Parameters:
        -----------
        uniprot_ids : list
            List of UniProt accessions
        chunk_size : int
            Number of IDs to query per batch

        Returns:
        --------
        dict
            Results for all queried IDs
        """
        results = {}

        # Process in chunks to avoid overwhelming the API
        for i in range(0, len(uniprot_ids), chunk_size):
            chunk = uniprot_ids[i:i + chunk_size]
            print(f"Processing batch {i//chunk_size + 1}/{(len(uniprot_ids)-1)//chunk_size + 1}")

            for uniprot_id in chunk:
                results[uniprot_id] = self.get_protein_info(uniprot_id)

        return results


class UPSAnalyzer:
    """
    Analyzer for Ubiquitin-Proteasome System (UPS) interactions and annotations
    """

    def __init__(self):
        """Initialize UPS analyzer with known UPS components"""
        self.ups_components = {
            'proteasome_20s': [
                'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
                'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7'
            ],
            'proteasome_19s': [
                'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
                'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6',
                'PSMD7', 'PSMD8', 'PSMD9', 'PSMD10', 'PSMD11', 'PSMD12',
                'PSMD13', 'PSMD14'
            ],
            'ubiquitin_enzymes': [
                'UBA1', 'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2',
                'UBE2D3', 'UBE2E1', 'UBE2E2', 'UBE2E3', 'UBE2F',
                'UBE2G1', 'UBE2G2', 'UBE2H', 'UBE2I', 'UBE2J1', 'UBE2J2',
                'UBE2K', 'UBE2L3', 'UBE2L6', 'UBE2M', 'UBE2N', 'UBE2O',
                'UBE2Q1', 'UBE2Q2', 'UBE2R2', 'UBE2S', 'UBE2T', 'UBE2U',
                'UBE2V1', 'UBE2V2', 'UBE2W', 'UBE2Z'
            ],
            'deubiquitinases': [
                'USP1', 'USP2', 'USP3', 'USP4', 'USP5', 'USP6', 'USP7',
                'USP8', 'USP9X', 'USP9Y', 'USP10', 'USP11', 'USP12',
                'USP13', 'USP14', 'USP15', 'USP16', 'USP17L2', 'USP18',
                'USP19', 'USP20', 'USP21', 'USP22', 'USP24', 'USP25',
                'USP26', 'USP27X', 'USP28', 'USP29', 'USP30', 'USP31',
                'USP32', 'USP33', 'USP34', 'USP35', 'USP36', 'USP37',
                'USP38', 'USP39', 'USP40', 'USP41', 'USP42', 'USP43',
                'USP44', 'USP45', 'USP46', 'USP47', 'USP48', 'USP49',
                'USP50', 'USP51', 'USP52', 'USP53', 'USP54'
            ],
            'shuttle_factors': [
                'SQSTM1', 'NBR1', 'OPTN', 'TAX1BP1', 'TOLLIP', 'NDP52'
            ]
        }

        self.ups_go_terms = {
            # Proteasome Complex Terms
            'GO:0000502': 'proteasome complex',
            'GO:0005839': '26S proteasome',
            'GO:0019773': 'proteasome core complex, alpha-subunit complex',
            'GO:0019774': 'proteasome core complex, beta-subunit complex',
            'GO:0008540': 'proteasome regulatory particle, base subcomplex',
            'GO:0008541': 'proteasome regulatory particle, lid subcomplex',
            'GO:0031597': 'cytosolic proteasome complex',
            'GO:0022625': 'cytosolic large ribosomal subunit',

            # Proteasome Activity Terms
            'GO:0008233': 'peptidase activity',
            'GO:0004298': 'threonine-type endopeptidase activity',
            'GO:0061133': 'endopeptidase activator activity',
            'GO:0004175': 'endopeptidase activity',
            'GO:0070628': 'proteasome binding',

            # Ubiquitin-Proteasome Pathway Terms
            'GO:0006511': 'ubiquitin-dependent protein catabolic process',
            'GO:0043161': 'proteasome-mediated ubiquitin-dependent protein catabolic process',
            'GO:0071630': 'nuclear protein quality control by the ubiquitin-proteasome system',
            'GO:0030433': 'ubiquitin-dependent ERAD pathway',
            'GO:0038061': 'NIK/NF-kappaB signaling',
            'GO:0043248': 'proteasome assembly',

            # Ubiquitination Terms
            'GO:0016567': 'protein ubiquitination',
            'GO:0070647': 'protein modification by small protein conjugation or removal',
            'GO:0032436': 'positive regulation of proteasomal ubiquitin-dependent protein catabolic process',
            'GO:0032435': 'negative regulation of proteasomal ubiquitin-dependent protein catabolic process',

            # E1 Enzyme Terms
            'GO:0019778': 'small protein activating enzyme activity',
            'GO:0004839': 'ubiquitin activating enzyme activity',

            # E2 Enzyme Terms
            'GO:0004842': 'ubiquitin-protein transferase activity',
            'GO:0061631': 'ubiquitin conjugating enzyme activity',

            # E3 Ligase Terms
            'GO:0061630': 'ubiquitin protein ligase activity',
            'GO:0031625': 'ubiquitin protein ligase binding',
            'GO:0008270': 'zinc ion binding',  # Many RING E3s
            'GO:0046872': 'metal ion binding',  # E3 ligases often require metal ions

            # Deubiquitinase Terms
            'GO:0016579': 'protein deubiquitination',
            'GO:0004843': 'thiol-dependent ubiquitin-specific protease activity',
            'GO:0101005': 'ubiquitinyl hydrolase activity',
            'GO:0036459': 'thiol-dependent ubiquitinyl hydrolase activity',

            # Quality Control Terms
            'GO:0030968': 'endoplasmic reticulum unfolded protein response',
            'GO:0034976': 'response to endoplasmic reticulum stress',
            'GO:0070417': 'cellular response to cold',
            'GO:0006950': 'response to stress',

            # Protein Folding and Chaperone Terms (UPS-associated)
            'GO:0006457': 'protein folding',
            'GO:0051082': 'unfolded protein binding',
            'GO:0042026': 'protein refolding',

            # Autophagy-UPS Crosstalk Terms
            'GO:0016236': 'macroautophagy',
            'GO:0000045': 'autophagosome assembly',
            'GO:0034727': 'piecemeal microautophagy of nucleus',

            # Cell Cycle and UPS Terms
            'GO:0051436': 'negative regulation of ubiquitin-protein ligase activity involved in mitotic cell cycle',
            'GO:0051437': 'positive regulation of ubiquitin-protein ligase activity involved in mitotic cell cycle',

            # Signaling Pathway Terms
            'GO:0038095': 'Fc-epsilon receptor signaling pathway',
            'GO:0043123': 'positive regulation of I-kappaB kinase/NF-kappaB signaling',

            # Substrate-specific UPS Terms
            'GO:0070936': 'protein K48-linked ubiquitination',
            'GO:0070534': 'protein K63-linked ubiquitination',
            'GO:0033522': 'histone H2A ubiquitination',
            'GO:0033523': 'histone H2B ubiquitination',

            # SQSTM1/p62-specific Terms (Shuttle Factor Functions)
            'GO:0006914': 'autophagy',
            'GO:0061912': 'selective autophagy',
            'GO:0034727': 'piecemeal microautophagy of nucleus',
            'GO:0000045': 'autophagosome assembly',
            'GO:0000422': 'autophagy of mitochondrion',
            'GO:0035973': 'aggrephagy',
            'GO:0071230': 'cellular response to amino acid stimulus',
            'GO:0016236': 'macroautophagy',
            'GO:0019779': 'Atg8 ligase activity',
            'GO:0045732': 'positive regulation of protein catabolic process',

            # SQSTM1 Interaction and Binding Terms
            'GO:0031625': 'ubiquitin protein ligase binding',
            'GO:0019899': 'enzyme binding',
            'GO:0008270': 'zinc ion binding',
            'GO:0043130': 'ubiquitin binding',
            'GO:0051082': 'unfolded protein binding',
            'GO:0030544': 'Hsp70 protein binding',

            # SQSTM1 Signaling and Regulation Terms
            'GO:0038095': 'Fc-epsilon receptor signaling pathway',
            'GO:0043123': 'positive regulation of I-kappaB kinase/NF-kappaB signaling',
            'GO:0006355': 'regulation of transcription, DNA-templated',
            'GO:0030968': 'endoplasmic reticulum unfolded protein response',
            'GO:0034976': 'response to endoplasmic reticulum stress',
            'GO:0071456': 'cellular response to hypoxia',

            # SQSTM1 Localization Terms
            'GO:0005813': 'centrosome',
            'GO:0005634': 'nucleus',
            'GO:0005737': 'cytoplasm',
            'GO:0005783': 'endoplasmic reticulum',
            'GO:0005776': 'autophagosome',
            'GO:0000932': 'P-body',

            # SQSTM1 Aggregate and Inclusion Body Terms
            'GO:0031396': 'regulation of protein ubiquitination',
            'GO:0048525': 'negative regulation of viral process',
            'GO:0070845': 'polyubiquitinated protein binding',
            'GO:0072594': 'establishment of protein localization to organelle'
        }

        self.ups_keywords = [
            # Core UPS Components
            'Proteasome', 'Ubiquitin', 'Protease', 'Peptidase',
            'Protein degradation', 'Ubl conjugation pathway',
            'Threonine protease', 'Deubiquitinase',

            # Proteasome-specific
            '26S proteasome', '20S proteasome', '19S proteasome',
            'Proteasome subunit', 'Proteasome regulatory particle',
            'Proteasome activator', 'Proteasome inhibitor',

            # Ubiquitin System
            'E1 enzyme', 'E2 enzyme', 'E3 ligase', 'E3 ubiquitin ligase',
            'Ubiquitin activating', 'Ubiquitin conjugating', 'Ubiquitin ligase',
            'RING finger', 'HECT domain', 'U-box domain', 'RBR domain',
            'Cullin', 'SCF complex', 'APC/C complex',

            # Deubiquitinases
            'DUB', 'USP', 'UCH', 'OTU', 'JAMM', 'MJD',
            'Ubiquitin specific protease', 'Ubiquitin carboxyl-terminal hydrolase',
            'Ovarian tumor protease', 'Machado-Joseph disease protease',

            # Quality Control
            'Protein quality control', 'ERAD', 'Misfolded protein',
            'Protein refolding', 'Chaperone', 'Heat shock protein',
            'Endoplasmic reticulum stress', 'Unfolded protein response',

            # Degradation Signals
            'Degron', 'N-degron', 'C-degron', 'Destruction box',
            'KEN box', 'PEST sequence', 'Ubiquitin signal',

            # Autophagy-UPS Crosstalk
            'Shuttle factor', 'Autophagy receptor', 'Selective autophagy',
            'Aggrephagy', 'Proteaphagy', 'P62', 'SQSTM1', 'NBR1', 'OPTN',

            # SQSTM1/p62-specific Terms
            'Sequestosome', 'p62 protein', 'SQSTM1 protein',
            'LC3 interacting region', 'LIR motif', 'UBA domain',
            'PB1 domain', 'ZZ domain', 'TRAF6 binding',
            'Keap1 binding', 'Polyubiquitin binding',
            'Autophagosome formation', 'Protein aggregate clearance',
            'Oxidative stress response', 'NRF2 pathway',
            'mTORC1 signaling', 'Nutrient sensing',

            # Pathway Regulation
            'Cell cycle', 'DNA damage response', 'Stress response',
            'NF-kappa B', 'p53 pathway', 'Apoptosis',
            'Signal transduction', 'Transcription regulation',

            # Disease-related
            'Neurodegeneration', 'Protein aggregation', 'Inclusion body',
            'Lewy body', 'Neurofibrillary tangle', 'Huntington',
            'Alzheimer', 'Parkinson', 'ALS', 'Prion'
        ]

    def classify_ups_component(self, gene_name: str, protein_info: Dict = None) -> Dict:
        """
        Classify protein as UPS component and determine its role

        Parameters:
        -----------
        gene_name : str
            Gene symbol
        protein_info : dict
            UniProt protein information (optional)

        Returns:
        --------
        dict
            UPS classification results
        """
        result = {
            'gene_name': gene_name,
            'is_ups_component': False,
            'ups_category': None,
            'ups_subcategory': None,
            'confidence': 'low',
            'evidence': []
        }

        # Direct component check
        for category, components in self.ups_components.items():
            if gene_name in components:
                result['is_ups_component'] = True
                result['ups_category'] = 'core_ups'
                result['ups_subcategory'] = category
                result['confidence'] = 'high'
                result['evidence'].append(f'Direct match in {category}')
                break

        # Check UniProt annotations if provided
        if protein_info and not result['is_ups_component']:
            ups_evidence = []

            # Check GO terms
            go_terms = protein_info.get('go_terms', [])
            for go_term in go_terms:
                go_id = go_term.get('id')
                if go_id in self.ups_go_terms:
                    ups_evidence.append(f"GO term: {go_id} ({self.ups_go_terms[go_id]})")

            # Check keywords
            keywords = protein_info.get('keywords', [])
            for keyword in keywords:
                if any(ups_kw.lower() in keyword.lower() for ups_kw in self.ups_keywords):
                    ups_evidence.append(f"Keyword: {keyword}")

            # Check protein name for UPS-related terms
            protein_name = protein_info.get('protein_name', '')
            if protein_name:
                for ups_kw in self.ups_keywords:
                    if ups_kw.lower() in protein_name.lower():
                        ups_evidence.append(f"Protein name contains: {ups_kw}")

            if ups_evidence:
                result['is_ups_component'] = True
                result['ups_category'] = 'ups_associated'
                result['confidence'] = 'medium' if len(ups_evidence) > 1 else 'low'
                result['evidence'] = ups_evidence

        return result

    def analyze_ups_interactions(self, protein_list: List[str],
                                interaction_threshold: float = 0.7) -> Dict:
        """
        Analyze interactions between proteins and known UPS components

        Parameters:
        -----------
        protein_list : list
            List of gene symbols to analyze
        interaction_threshold : float
            Minimum confidence score for interactions

        Returns:
        --------
        dict
            UPS interaction analysis results
        """
        # This would typically query STRING, IntAct, or other databases
        # For now, return framework for implementation

        results = {
            'input_proteins': protein_list,
            'ups_classifications': {},
            'interaction_network': {},
            'summary': {
                'total_proteins': len(protein_list),
                'ups_components': 0,
                'ups_associated': 0,
                'high_confidence_interactions': 0
            }
        }

        # Classify each protein
        for gene in protein_list:
            classification = self.classify_ups_component(gene)
            results['ups_classifications'][gene] = classification

            if classification['is_ups_component']:
                if classification['ups_category'] == 'core_ups':
                    results['summary']['ups_components'] += 1
                else:
                    results['summary']['ups_associated'] += 1

        return results


class ProteinValidator:
    """
    Comprehensive protein annotation validator
    """

    def __init__(self):
        """Initialize validator with UniProt API client"""
        self.uniprot_api = UniProtAPI()
        self.ups_analyzer = UPSAnalyzer()

    def validate_protein_dataset(self, protein_df: pd.DataFrame,
                                uniprot_col: str = 'uniprot_id',
                                gene_col: str = 'gene_symbol') -> pd.DataFrame:
        """
        Validate entire protein dataset against UniProt

        Parameters:
        -----------
        protein_df : pd.DataFrame
            Protein annotation dataframe
        uniprot_col : str
            Column name containing UniProt IDs
        gene_col : str
            Column name containing gene symbols

        Returns:
        --------
        pd.DataFrame
            Validation results
        """
        results = []

        print(f"Validating {len(protein_df)} proteins against UniProt...")

        for idx, row in protein_df.iterrows():
            uniprot_id = row[uniprot_col]
            gene_symbol = row[gene_col]

            print(f"Validating {gene_symbol} ({uniprot_id})...")

            # Handle multiple UniProt IDs
            primary_id, all_ids = self._handle_multiple_ids(uniprot_id)

            # Validate primary ID
            validation = self.uniprot_api.validate_uniprot_id(primary_id)
            protein_info = self.uniprot_api.get_protein_info(primary_id)

            # UPS classification
            ups_classification = self.ups_analyzer.classify_ups_component(
                gene_symbol, protein_info
            )

            # Cross-validation checks
            validation_issues = self._cross_validate_annotations(
                row, validation, protein_info
            )

            result = {
                'original_uniprot_id': uniprot_id,
                'primary_uniprot_id': primary_id,
                'has_multiple_ids': len(all_ids) > 1,
                'all_uniprot_ids': all_ids,
                'gene_symbol': gene_symbol,
                'uniprot_exists': validation['exists'],
                'uniprot_reviewed': validation['is_reviewed'],
                'gene_name_match': validation.get('gene_name') == gene_symbol,
                'uniprot_gene_name': validation.get('gene_name'),
                'uniprot_protein_name': validation.get('protein_name'),
                'is_ups_component': ups_classification['is_ups_component'],
                'ups_category': ups_classification['ups_category'],
                'ups_confidence': ups_classification['confidence'],
                'validation_issues': validation_issues,
                'validation_status': self._determine_validation_status(
                    validation, validation_issues
                )
            }

            results.append(result)

        return pd.DataFrame(results)

    def _handle_multiple_ids(self, uniprot_string: str) -> Tuple[str, List[str]]:
        """
        Handle multiple UniProt IDs separated by semicolons

        Parameters:
        -----------
        uniprot_string : str
            UniProt ID string (may contain multiple IDs)

        Returns:
        --------
        tuple
            (primary_id, all_ids)
        """
        if pd.isna(uniprot_string) or not uniprot_string:
            return None, []

        # Split on semicolon and clean
        ids = [id_str.strip() for id_str in str(uniprot_string).split(';')]
        ids = [id_str for id_str in ids if id_str]  # Remove empty strings

        if not ids:
            return None, []

        # First ID is considered primary
        return ids[0], ids

    def _cross_validate_annotations(self, row: pd.Series,
                                   validation: Dict,
                                   protein_info: Dict) -> List[str]:
        """
        Cross-validate annotations between dataset and UniProt

        Parameters:
        -----------
        row : pd.Series
            Original protein annotation row
        validation : dict
            UniProt validation results
        protein_info : dict
            Detailed protein information from UniProt

        Returns:
        --------
        list
            List of validation issues found
        """
        issues = []

        # Gene symbol mismatch
        if (validation.get('gene_name') and
            validation['gene_name'] != row.get('gene_symbol')):
            issues.append(
                f"Gene symbol mismatch: dataset={row.get('gene_symbol')}, "
                f"UniProt={validation['gene_name']}"
            )

        # Molecular weight check (allow 10% tolerance)
        if ('molecular_weight' in row and
            protein_info.get('molecular_weight')):
            dataset_mw = row['molecular_weight']
            uniprot_mw = protein_info['molecular_weight'] / 1000  # Convert to kDa

            if abs(dataset_mw - uniprot_mw) / uniprot_mw > 0.1:
                issues.append(
                    f"Molecular weight mismatch: dataset={dataset_mw}kDa, "
                    f"UniProt={uniprot_mw:.1f}kDa"
                )

        # Protein name similarity check
        if ('protein_name' in row and
            validation.get('protein_name')):
            dataset_name = str(row['protein_name']).lower()
            uniprot_name = validation['protein_name'].lower()

            # Simple similarity check
            if (dataset_name not in uniprot_name and
                uniprot_name not in dataset_name):
                issues.append(
                    f"Protein name mismatch: dataset='{row['protein_name']}', "
                    f"UniProt='{validation['protein_name']}'"
                )

        return issues

    def _determine_validation_status(self, validation: Dict,
                                   issues: List[str]) -> str:
        """
        Determine overall validation status

        Parameters:
        -----------
        validation : dict
            UniProt validation results
        issues : list
            List of validation issues

        Returns:
        --------
        str
            Validation status: 'valid', 'warning', 'invalid'
        """
        if not validation['exists']:
            return 'invalid'

        if len(issues) == 0:
            return 'valid'

        # Check for critical issues
        critical_issues = [
            issue for issue in issues
            if 'gene symbol mismatch' in issue.lower()
        ]

        if critical_issues:
            return 'invalid'

        return 'warning'

    def generate_validation_report(self, validation_df: pd.DataFrame,
                                 output_file: str = None) -> Dict:
        """
        Generate comprehensive validation report

        Parameters:
        -----------
        validation_df : pd.DataFrame
            Validation results dataframe
        output_file : str
            Optional file to save detailed report

        Returns:
        --------
        dict
            Summary statistics and recommendations
        """
        report = {
            'summary': {
                'total_proteins': len(validation_df),
                'valid_proteins': len(validation_df[validation_df['validation_status'] == 'valid']),
                'warning_proteins': len(validation_df[validation_df['validation_status'] == 'warning']),
                'invalid_proteins': len(validation_df[validation_df['validation_status'] == 'invalid']),
                'multiple_id_proteins': len(validation_df[validation_df['has_multiple_ids']]),
                'ups_components': len(validation_df[validation_df['is_ups_component']]),
                'reviewed_entries': len(validation_df[validation_df['uniprot_reviewed']])
            },
            'recommendations': [],
            'critical_issues': [],
            'ups_analysis': {
                'core_ups': len(validation_df[validation_df['ups_category'] == 'core_ups']),
                'ups_associated': len(validation_df[validation_df['ups_category'] == 'ups_associated'])
            }
        }

        # Generate recommendations
        if report['summary']['invalid_proteins'] > 0:
            report['recommendations'].append(
                f"Review {report['summary']['invalid_proteins']} invalid protein entries"
            )

        if report['summary']['multiple_id_proteins'] > 0:
            report['recommendations'].append(
                f"Consider handling {report['summary']['multiple_id_proteins']} proteins with multiple IDs"
            )

        if report['summary']['reviewed_entries'] < report['summary']['total_proteins'] * 0.8:
            report['recommendations'].append(
                "Low percentage of reviewed UniProt entries - consider using Swiss-Prot IDs when possible"
            )

        # Identify critical issues
        critical_df = validation_df[validation_df['validation_status'] == 'invalid']
        for _, row in critical_df.iterrows():
            report['critical_issues'].append({
                'gene_symbol': row['gene_symbol'],
                'uniprot_id': row['original_uniprot_id'],
                'issues': row['validation_issues']
            })

        # Save detailed report if requested
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)

        return report


def main():
    """
    Example usage of UniProt analysis tools
    """
    # Load sample data
    protein_df = pd.read_csv('/Users/byron/project_plan/biologist_guide_proteomics/data/protein_annotations.csv')

    # Initialize validator
    validator = ProteinValidator()

    # Validate dataset
    print("Starting protein validation...")
    validation_results = validator.validate_protein_dataset(protein_df)

    # Generate report
    print("Generating validation report...")
    report = validator.generate_validation_report(validation_results)

    # Display summary
    print("\n=== VALIDATION SUMMARY ===")
    for key, value in report['summary'].items():
        print(f"{key}: {value}")

    print("\n=== UPS ANALYSIS ===")
    for key, value in report['ups_analysis'].items():
        print(f"{key}: {value}")

    if report['recommendations']:
        print("\n=== RECOMMENDATIONS ===")
        for rec in report['recommendations']:
            print(f"- {rec}")

    # Save results
    validation_results.to_csv('protein_validation_results.csv', index=False)
    print("\nResults saved to protein_validation_results.csv")


if __name__ == "__main__":
    main()
```
