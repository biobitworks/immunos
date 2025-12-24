---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/ups_interaction_analysis.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/ups_interaction_analysis.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
UPS Interaction Analysis Tool

Comprehensive analysis of protein interactions with the Ubiquitin-Proteasome System (UPS).
Integrates multiple databases and analysis methods to identify UPS components and interactions.

Author: Proteomics Analysis Framework
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import time
from typing import Dict, List, Tuple, Optional
import json
from urllib.parse import urlencode
import warnings
warnings.filterwarnings('ignore')

from uniprot_analysis import UPSAnalyzer, UniProtAPI


class UPSInteractionDatabase:
    """
    Comprehensive UPS interaction database with multiple data sources
    """

    def __init__(self):
        """Initialize UPS interaction database with known interactions"""
        self.string_api_url = "https://string-db.org/api"
        self.biogrid_api_url = "https://webservice.thebiogrid.org"

        # Known UPS protein complexes and interactions
        self.ups_complexes = self._initialize_ups_complexes()
        self.ups_pathways = self._initialize_ups_pathways()

    def _initialize_ups_complexes(self) -> Dict:
        """Initialize known UPS protein complexes"""
        return {
            '26S_proteasome': {
                '20S_core': {
                    'alpha_subunits': ['PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7'],
                    'beta_subunits': ['PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7']
                },
                '19S_regulatory': {
                    'base_complex': ['PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6'],
                    'lid_complex': ['PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6',
                                   'PSMD7', 'PSMD8', 'PSMD9', 'PSMD10', 'PSMD11', 'PSMD12',
                                   'PSMD13', 'PSMD14']
                }
            },
            'SCF_complex': {
                'core_components': ['SKP1', 'CUL1', 'RBX1'],
                'f_box_proteins': ['FBXW7', 'FBXO4', 'FBXL3', 'FBXW11', 'FBXO31']
            },
            'APC_complex': {
                'core_components': ['APC', 'APC2', 'APC3', 'APC4', 'APC5', 'APC6', 'APC7',
                                   'APC8', 'APC10', 'APC11', 'APC12', 'APC13', 'APC15', 'APC16'],
                'activators': ['CDC20', 'CDH1']
            },
            'CRL4_complex': {
                'core_components': ['CUL4A', 'CUL4B', 'DDB1', 'RBX1'],
                'substrate_receptors': ['DDB2', 'XPC', 'CSA', 'WDR23']
            }
        }

    def _initialize_ups_pathways(self) -> Dict:
        """Initialize UPS pathway information"""
        return {
            'proteasomal_degradation': {
                'pathway_id': 'hsa03050',
                'description': '26S proteasome-mediated protein degradation',
                'key_proteins': ['UBA1', 'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3']
            },
            'ubiquitin_conjugation': {
                'pathway_id': 'hsa04120',
                'description': 'Ubiquitin mediated proteolysis',
                'key_proteins': ['UBA1', 'UBA2', 'UBA3', 'UBA6']
            },
            'ERAD_pathway': {
                'pathway_id': 'hsa04141',
                'description': 'Protein processing in endoplasmic reticulum',
                'key_proteins': ['SEL1L', 'HRD1', 'DERL1', 'DERL2', 'VIMP']
            },
            'DNA_damage_response': {
                'pathway_id': 'hsa03460',
                'description': 'Fanconi anemia pathway',
                'key_proteins': ['BRCA1', 'BRCA2', 'TP53', 'MDM2', 'RNF8', 'RNF168']
            }
        }

    def query_string_interactions(self, proteins: List[str],
                                species: int = 9606,
                                score_threshold: float = 0.7) -> pd.DataFrame:
        """
        Query STRING database for protein interactions

        Parameters:
        -----------
        proteins : list
            List of gene symbols
        species : int
            NCBI taxonomy ID (9606 = Homo sapiens)
        score_threshold : float
            Minimum confidence score (0-1)

        Returns:
        --------
        pd.DataFrame
            Interaction data from STRING
        """
        interactions = []

        print(f"Querying STRING database for {len(proteins)} proteins...")

        # Query in batches to avoid overwhelming the API
        batch_size = 50
        for i in range(0, len(proteins), batch_size):
            batch = proteins[i:i + batch_size]

            try:
                # Build query URL
                params = {
                    'identifiers': '%0d'.join(batch),
                    'species': species,
                    'required_score': int(score_threshold * 1000),
                    'network_type': 'functional'
                }

                url = f"{self.string_api_url}/tsv/network"

                # Make request with rate limiting
                time.sleep(1)  # Respect API limits
                response = requests.get(url, params=params, timeout=30)
                response.raise_for_status()

                # Parse response
                lines = response.text.strip().split('\n')
                if len(lines) > 1:  # Skip header
                    for line in lines[1:]:
                        parts = line.split('\t')
                        if len(parts) >= 16:
                            interactions.append({
                                'protein_a': parts[2],
                                'protein_b': parts[3],
                                'combined_score': float(parts[5]) / 1000,
                                'string_id_a': parts[0],
                                'string_id_b': parts[1]
                            })

            except Exception as e:
                print(f"Error querying STRING for batch {i//batch_size + 1}: {e}")
                continue

        if interactions:
            return pd.DataFrame(interactions)
        else:
            return pd.DataFrame(columns=['protein_a', 'protein_b', 'combined_score',
                                       'string_id_a', 'string_id_b'])

    def identify_ups_interactions(self, proteins: List[str]) -> Dict:
        """
        Identify UPS interactions for a list of proteins

        Parameters:
        -----------
        proteins : list
            List of gene symbols to analyze

        Returns:
        --------
        dict
            Comprehensive UPS interaction analysis
        """
        print("Performing comprehensive UPS interaction analysis...")

        # Initialize UPS analyzer
        ups_analyzer = UPSAnalyzer()

        # Classify all proteins for UPS involvement
        ups_classifications = {}
        ups_proteins = []

        for protein in proteins:
            classification = ups_analyzer.classify_ups_component(protein)
            ups_classifications[protein] = classification

            if classification['is_ups_component']:
                ups_proteins.append(protein)

        # Get interaction data from STRING
        string_interactions = self.query_string_interactions(proteins)

        # Analyze UPS-specific interactions
        ups_interactions = self._analyze_ups_specific_interactions(
            ups_proteins, string_interactions, ups_classifications
        )

        # Build UPS network
        ups_network = self._build_ups_network(ups_interactions, ups_classifications)

        # Calculate network metrics
        network_metrics = self._calculate_network_metrics(ups_network)

        # Identify functional modules
        functional_modules = self._identify_functional_modules(
            ups_network, ups_classifications
        )

        return {
            'input_proteins': proteins,
            'ups_classifications': ups_classifications,
            'ups_proteins': ups_proteins,
            'string_interactions': string_interactions,
            'ups_interactions': ups_interactions,
            'network_metrics': network_metrics,
            'functional_modules': functional_modules,
            'summary': self._generate_interaction_summary(
                proteins, ups_proteins, ups_interactions, network_metrics
            )
        }

    def _analyze_ups_specific_interactions(self, ups_proteins: List[str],
                                         string_interactions: pd.DataFrame,
                                         classifications: Dict) -> pd.DataFrame:
        """Identify interactions specifically involving UPS components"""

        if string_interactions.empty:
            return pd.DataFrame()

        # Filter for interactions involving at least one UPS protein
        ups_mask = (
            string_interactions['protein_a'].isin(ups_proteins) |
            string_interactions['protein_b'].isin(ups_proteins)
        )

        ups_specific = string_interactions[ups_mask].copy()

        if ups_specific.empty:
            return ups_specific

        # Add UPS classification information
        ups_specific['protein_a_ups_category'] = ups_specific['protein_a'].map(
            lambda x: classifications.get(x, {}).get('ups_category', 'non_ups')
        )
        ups_specific['protein_b_ups_category'] = ups_specific['protein_b'].map(
            lambda x: classifications.get(x, {}).get('ups_category', 'non_ups')
        )

        # Classify interaction types
        interaction_types = []
        for _, row in ups_specific.iterrows():
            cat_a = row['protein_a_ups_category']
            cat_b = row['protein_b_ups_category']

            if cat_a == 'core_ups' and cat_b == 'core_ups':
                interaction_types.append('core_ups_internal')
            elif (cat_a == 'core_ups' and cat_b == 'ups_associated') or \
                 (cat_a == 'ups_associated' and cat_b == 'core_ups'):
                interaction_types.append('core_ups_associated')
            elif (cat_a in ['core_ups', 'ups_associated'] and cat_b == 'non_ups') or \
                 (cat_a == 'non_ups' and cat_b in ['core_ups', 'ups_associated']):
                interaction_types.append('ups_substrate')
            else:
                interaction_types.append('other')

        ups_specific['interaction_type'] = interaction_types

        return ups_specific.sort_values('combined_score', ascending=False)

    def _build_ups_network(self, ups_interactions: pd.DataFrame,
                          classifications: Dict) -> nx.Graph:
        """Build UPS interaction network"""

        G = nx.Graph()

        # Add all proteins as nodes
        for protein, info in classifications.items():
            G.add_node(protein,
                      ups_component=info['is_ups_component'],
                      ups_category=info.get('ups_category', 'non_ups'),
                      ups_confidence=info.get('ups_confidence', 'none'))

        # Add interactions as edges
        if not ups_interactions.empty:
            for _, row in ups_interactions.iterrows():
                G.add_edge(row['protein_a'], row['protein_b'],
                          score=row['combined_score'],
                          interaction_type=row['interaction_type'])

        return G

    def _calculate_network_metrics(self, network: nx.Graph) -> Dict:
        """Calculate comprehensive network metrics"""

        if len(network.nodes()) == 0:
            return {
                'nodes': 0, 'edges': 0, 'density': 0, 'clustering': 0,
                'connected_components': 0, 'largest_component_size': 0
            }

        metrics = {
            'nodes': len(network.nodes()),
            'edges': len(network.edges()),
            'density': nx.density(network),
            'clustering': nx.average_clustering(network),
            'connected_components': nx.number_connected_components(network),
            'largest_component_size': len(max(nx.connected_components(network),
                                            key=len, default=[]))
        }

        # Calculate centrality measures for UPS proteins
        ups_nodes = [node for node, data in network.nodes(data=True)
                    if data.get('ups_component', False)]

        if ups_nodes and len(network.edges()) > 0:
            degree_centrality = nx.degree_centrality(network)
            betweenness_centrality = nx.betweenness_centrality(network)

            metrics['ups_hub_proteins'] = {
                'degree': sorted([(node, degree_centrality[node]) for node in ups_nodes],
                               key=lambda x: x[1], reverse=True)[:5],
                'betweenness': sorted([(node, betweenness_centrality[node]) for node in ups_nodes],
                                    key=lambda x: x[1], reverse=True)[:5]
            }

        return metrics

    def _identify_functional_modules(self, network: nx.Graph,
                                   classifications: Dict) -> Dict:
        """Identify functional modules in the UPS network"""

        modules = {
            'proteasome_module': [],
            'e3_ligase_module': [],
            'deubiquitinase_module': [],
            'substrate_module': [],
            'mixed_modules': []
        }

        # Get connected components
        components = list(nx.connected_components(network))

        for component in components:
            component_list = list(component)

            # Classify component based on UPS categories
            ups_categories = [classifications.get(protein, {}).get('ups_category', 'non_ups')
                            for protein in component_list]

            core_ups_count = ups_categories.count('core_ups')
            ups_associated_count = ups_categories.count('ups_associated')

            if core_ups_count > len(component_list) * 0.7:
                modules['proteasome_module'].extend(component_list)
            elif ups_associated_count > len(component_list) * 0.5:
                modules['e3_ligase_module'].extend(component_list)
            elif core_ups_count + ups_associated_count > len(component_list) * 0.3:
                modules['mixed_modules'].append(component_list)
            else:
                modules['substrate_module'].extend(component_list)

        return modules

    def _generate_interaction_summary(self, all_proteins: List[str],
                                    ups_proteins: List[str],
                                    ups_interactions: pd.DataFrame,
                                    network_metrics: Dict) -> Dict:
        """Generate comprehensive interaction analysis summary"""

        summary = {
            'total_proteins': len(all_proteins),
            'ups_proteins': len(ups_proteins),
            'ups_coverage': len(ups_proteins) / len(all_proteins) * 100,
            'total_interactions': len(ups_interactions),
            'network_density': network_metrics.get('density', 0),
            'largest_component': network_metrics.get('largest_component_size', 0)
        }

        if not ups_interactions.empty:
            interaction_types = ups_interactions['interaction_type'].value_counts()
            summary['interaction_breakdown'] = interaction_types.to_dict()

            # High-confidence interactions
            high_conf = ups_interactions[ups_interactions['combined_score'] > 0.8]
            summary['high_confidence_interactions'] = len(high_conf)

        return summary


class UPSVisualization:
    """
    Visualization tools for UPS interaction analysis
    """

    def __init__(self):
        """Initialize visualization tools"""
        plt.style.use('seaborn-v0_8-darkgrid')

    def plot_ups_network(self, network: nx.Graph, output_file: str = None,
                        figsize: Tuple[int, int] = (12, 10)) -> None:
        """
        Create comprehensive UPS network visualization

        Parameters:
        -----------
        network : nx.Graph
            UPS interaction network
        output_file : str
            Optional file to save plot
        figsize : tuple
            Figure size
        """

        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('UPS Interaction Network Analysis', fontsize=16, fontweight='bold')

        # 1. Full network layout
        ax = axes[0, 0]
        pos = nx.spring_layout(network, k=3, iterations=50)

        # Color nodes by UPS category
        node_colors = []
        for node in network.nodes():
            category = network.nodes[node].get('ups_category', 'non_ups')
            if category == 'core_ups':
                node_colors.append('red')
            elif category == 'ups_associated':
                node_colors.append('orange')
            else:
                node_colors.append('lightgray')

        nx.draw(network, pos, ax=ax, node_color=node_colors,
               node_size=50, edge_color='gray', alpha=0.7, width=0.5)
        ax.set_title('Full UPS Network')

        # Add legend
        legend_elements = [
            plt.scatter([], [], c='red', s=50, label='Core UPS'),
            plt.scatter([], [], c='orange', s=50, label='UPS Associated'),
            plt.scatter([], [], c='lightgray', s=50, label='Other')
        ]
        ax.legend(handles=legend_elements, loc='upper right')

        # 2. Degree distribution
        ax = axes[0, 1]
        degrees = [network.degree(node) for node in network.nodes()]
        ax.hist(degrees, bins=20, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Node Degree')
        ax.set_ylabel('Frequency')
        ax.set_title('Degree Distribution')

        # 3. UPS component subnetwork
        ax = axes[1, 0]
        ups_nodes = [node for node, data in network.nodes(data=True)
                    if data.get('ups_component', False)]

        if ups_nodes:
            ups_subgraph = network.subgraph(ups_nodes)
            ups_pos = nx.spring_layout(ups_subgraph, k=2)

            ups_node_colors = []
            for node in ups_subgraph.nodes():
                category = ups_subgraph.nodes[node].get('ups_category', 'non_ups')
                if category == 'core_ups':
                    ups_node_colors.append('darkred')
                else:
                    ups_node_colors.append('darkorange')

            nx.draw(ups_subgraph, ups_pos, ax=ax, node_color=ups_node_colors,
                   node_size=100, edge_color='darkgray', alpha=0.8, width=1)

            # Add labels for important nodes
            high_degree_nodes = [node for node in ups_subgraph.nodes()
                               if ups_subgraph.degree(node) >= 3]
            labels = {node: node for node in high_degree_nodes}
            nx.draw_networkx_labels(ups_subgraph, ups_pos, labels, ax=ax, font_size=8)

        ax.set_title('UPS Components Subnetwork')

        # 4. Interaction confidence distribution
        ax = axes[1, 1]
        if network.edges():
            edge_scores = [data.get('score', 0) for _, _, data in network.edges(data=True)]
            ax.hist(edge_scores, bins=20, edgecolor='black', alpha=0.7, color='skyblue')
            ax.axvline(x=0.7, color='red', linestyle='--', label='High confidence (0.7)')
            ax.set_xlabel('Interaction Confidence Score')
            ax.set_ylabel('Frequency')
            ax.set_title('Interaction Confidence Distribution')
            ax.legend()

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')

        plt.show()

    def plot_ups_analysis_summary(self, analysis_results: Dict,
                                output_file: str = None) -> None:
        """
        Create summary visualization of UPS analysis results

        Parameters:
        -----------
        analysis_results : dict
            Results from UPS interaction analysis
        output_file : str
            Optional file to save plot
        """

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('UPS Interaction Analysis Summary', fontsize=16, fontweight='bold')

        summary = analysis_results['summary']
        classifications = analysis_results['ups_classifications']

        # 1. UPS Coverage
        ax = axes[0, 0]
        ups_data = {
            'UPS Proteins': summary['ups_proteins'],
            'Other Proteins': summary['total_proteins'] - summary['ups_proteins']
        }
        colors = ['lightcoral', 'lightgray']
        wedges, texts, autotexts = ax.pie(ups_data.values(), labels=ups_data.keys(),
                                         autopct='%1.1f%%', colors=colors)
        ax.set_title(f"UPS Coverage ({summary['ups_coverage']:.1f}%)")

        # 2. UPS Categories
        ax = axes[0, 1]
        categories = {}
        for protein, info in classifications.items():
            if info['is_ups_component']:
                cat = info.get('ups_category', 'unknown')
                categories[cat] = categories.get(cat, 0) + 1

        if categories:
            ax.bar(categories.keys(), categories.values(), color=['darkred', 'darkorange'])
            ax.set_xlabel('UPS Category')
            ax.set_ylabel('Number of Proteins')
            ax.set_title('UPS Component Categories')
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # 3. Confidence Distribution
        ax = axes[0, 2]
        confidences = {}
        for protein, info in classifications.items():
            if info['is_ups_component']:
                conf = info.get('ups_confidence', 'none')
                confidences[conf] = confidences.get(conf, 0) + 1

        if confidences:
            conf_colors = {'high': 'darkgreen', 'medium': 'orange', 'low': 'lightcoral'}
            colors = [conf_colors.get(conf, 'gray') for conf in confidences.keys()]
            ax.bar(confidences.keys(), confidences.values(), color=colors)
            ax.set_xlabel('Confidence Level')
            ax.set_ylabel('Number of Proteins')
            ax.set_title('UPS Classification Confidence')

        # 4. Interaction Types
        ax = axes[1, 0]
        if 'interaction_breakdown' in summary:
            interaction_data = summary['interaction_breakdown']
            ax.pie(interaction_data.values(), labels=interaction_data.keys(), autopct='%1.1f%%')
            ax.set_title('Interaction Type Distribution')

        # 5. Network Metrics
        ax = axes[1, 1]
        metrics = analysis_results['network_metrics']
        metric_names = ['Nodes', 'Edges', 'Density', 'Clustering']
        metric_values = [
            metrics.get('nodes', 0),
            metrics.get('edges', 0),
            metrics.get('density', 0) * 100,  # Convert to percentage
            metrics.get('clustering', 0) * 100  # Convert to percentage
        ]

        bars = ax.bar(metric_names, metric_values, color='steelblue')
        ax.set_ylabel('Value')
        ax.set_title('Network Metrics')

        # Add value labels on bars
        for bar, value in zip(bars, metric_values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(metric_values)*0.01,
                   f'{value:.1f}', ha='center', va='bottom')

        # 6. Top Hub Proteins
        ax = axes[1, 2]
        if 'ups_hub_proteins' in metrics and metrics['ups_hub_proteins']['degree']:
            hub_data = metrics['ups_hub_proteins']['degree'][:5]
            proteins = [item[0] for item in hub_data]
            centralities = [item[1] for item in hub_data]

            ax.barh(proteins, centralities, color='lightgreen')
            ax.set_xlabel('Degree Centrality')
            ax.set_title('Top UPS Hub Proteins')

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')

        plt.show()


def main():
    """
    Example usage of UPS interaction analysis tools
    """
    # Load sample data
    protein_df = pd.read_csv('/Users/byron/project_plan/biologist_guide_proteomics/data/protein_annotations.csv')
    proteins = protein_df['gene_symbol'].tolist()

    # Initialize UPS interaction database
    ups_db = UPSInteractionDatabase()

    # Perform comprehensive UPS analysis
    print("Starting UPS interaction analysis...")
    analysis_results = ups_db.identify_ups_interactions(proteins)

    # Create visualizations
    visualizer = UPSVisualization()

    # Plot network if interactions found
    if analysis_results['network_metrics']['edges'] > 0:
        print("Creating network visualization...")
        network = ups_db._build_ups_network(
            analysis_results['ups_interactions'],
            analysis_results['ups_classifications']
        )
        visualizer.plot_ups_network(network, 'ups_network.png')

    # Plot analysis summary
    print("Creating summary visualization...")
    visualizer.plot_ups_analysis_summary(analysis_results, 'ups_analysis_summary.png')

    # Print summary
    print("\n=== UPS INTERACTION ANALYSIS SUMMARY ===")
    summary = analysis_results['summary']
    for key, value in summary.items():
        if isinstance(value, dict):
            print(f"{key}:")
            for subkey, subvalue in value.items():
                print(f"  {subkey}: {subvalue}")
        else:
            print(f"{key}: {value}")

    # Save detailed results
    with open('ups_interaction_analysis.json', 'w') as f:
        # Convert NetworkX graph to serializable format
        serializable_results = analysis_results.copy()
        serializable_results.pop('network_metrics', None)  # Remove non-serializable network
        json.dump(serializable_results, f, indent=2, default=str)

    print("\nDetailed results saved to ups_interaction_analysis.json")
    print("Analysis complete!")


if __name__ == "__main__":
    main()
```
