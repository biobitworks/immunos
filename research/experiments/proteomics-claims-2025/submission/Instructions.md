#   
Instructions

Welcome to the findings evaluation project! You will be evaluating the outputs of an automated system for data analysis and data-driven discovery. _Please read the instructions below very carefully._

The output of the system consists of a brief report summarizing the main finding, with quantitative details and figures from the analyses. This report contains both results from direct analyses as well as scientific interpretations of the data. We have chosen a random sample of these statements (both based on data analysis and interpretations) for your evaluation:

1. Read through the PDF for the discovery group. This report will provide context for the statements that you will be evaluating. The figures are there for your reference, but they will not be evaluated. _Do not look through the hyperlinked citations or use those as evidence._
2. Familiarize yourself with the dataset provided. All analyses and discoveries in this discovery group were derived from this dataset. _Use jupyter notebook at the link provided for all analyses._
3. Each discovery group will have a sample of statements to be evaluated. These claims are either derived from data analysis on the dataset or by interpreting ideas across claims/conclusions. Conduct your evaluation based on the type of claim it is.
4. CLAIMING A GROUP: If you would like to claim a group for you

## **Dataset Information**

This is a proteomic dataset of mini pools of 10 neurons from Alzheimer's disease cases. MC1 (i.e. misfolded tau) quantification can be used to stratify tau positive cells or the precomputed pseudotime to order samples. The obs dataframe contains metadata about samples (age at death, MC1 score, pseudotime), and the tau status (positive/negative). Robust statistical approaches should be used to account for the possible cofactors. Data are already log2 transformed.

## Evaluation Rubrics

For each statement, you are asked to select the appropriate dropdown that best matches your assessment of the claim's support or lack thereof. See an explanation of the rubric below.

### Data Analysis Rubric

- Attempt to replicate this claim using your own data analysis. Did your findings support or contradict the claim?
- The exact numbers and results do not need to match exactly, but your analysis should address the central idea behind the claim. Use your best judgement to determine if the claim is supported, based on the context of the discovery.
- For example, a claim may say that Gene X expression was upregulated by 21-fold, but your analysis shows upregulation of 27-fold, the validity of the claim depends on if the central idea behind the claim is still supported or if the numerical difference refutes the idea. If the statistical significance is different but the direction of the claim is the same, the claim should be marked as supported—we are not accounting for differences in statistical significance.
- **IMPORTANT: Do all data analysis in PYTHON. If you need to conduct single-cell analysis, use the** `**scanpy**` **package.**
- Assessment options: [SUPPORTED, REFUTED, UNSURE]
- SUPPORTED: My data analysis is able to reproduce this claim.
- REFUTED: My data analysis refutes this claim.
- UNSURE: I do not know how to do this analysis, the results are unclear/insignificant, or the claim is ambiguous or missing key information.
- Explanation: Provide a free text explanation for you analysis and assessment. Include as much evidence and explanation as you think is necessary to justify your choice.

### Interpretation Rubric

- Analyze the claim and determine if the claim is a scientifically accurate interpretation based only on the report context. This question is only evaluating if the claim’s interpretation of the results is accurate, so you should not conduct any data analysis. A summary of relevant context is given for your reference, which you should assume is correct.
- Assessment options: [SUPPORTED, REFUTED, UNSURE]
- SUPPORTED: This claim is a logically reasonable inference given the report and context.
- REFUTED: This claim is NOT a logically reasonable inference given the report and context.
- UNSURE: The claim is ambiguous, missing key information, or unable to be evaluated.
- Explanation: Provide a free text explanation for your assessment. Include as much evidence and explanation as you think is necessary to justify your choice.