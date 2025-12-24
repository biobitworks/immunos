---
source: /Users/byron/projects/experiments/benford_pi.py
relative: experiments/benford_pi.py
generated_at: 2025-12-23 10:28
---

```python
import sys
import numpy as np
from mpmath import mp
from scipy.stats import chisquare
import jax
import jax.numpy as jnp
from thrml import CategoricalNode
from thrml.block_management import Block
from thrml.block_sampling import BlockGibbsSpec, sample_states, SamplingSchedule
from thrml.models.discrete_ebm import CategoricalEBM
from thrml.sampling_programs import CategoricalSamplingProgram

# User input for number of digits
X = int(input("Enter the number of digits after the decimal point for pi: "))

# Compute pi to X digits using mpmath
mp.dps = X + 20  # Extra digits for accuracy
pi_str = str(mp.pi)[2:2 + X]  # Extract decimal digits (skip '3.')
digits = [int(d) for d in pi_str]

# Parameters for extracting leading digits
l = 4  # Digits per number (adjust for more/less samples; e.g., l=1 for single digits)
step = 1  # Step size for overlapping windows

observed_leading = []
for i in range(0, len(digits) - l + 1, step):
    num = 0
    for j in range(l):
        num = num * 10 + digits[i + j]
    if num >= 10 ** (l - 1):  # Only consider numbers with l significant digits (no leading zero)
        leading = num // 10 ** (l - 1)
        observed_leading.append(leading)

observed = np.array(observed_leading)
num_obs = len(observed)
print(f"Generated {num_obs} leading digit observations from pi's decimal expansion.")

if num_obs == 0:
    print("No valid observations. Try smaller l or larger X.")
    sys.exit(1)

# Benford's Law expected probabilities for digits 1-9
benford_probs = np.log10(1 + 1 / np.arange(1, 10))

# Observed frequencies (proportions)
obs_counts, _ = np.histogram(observed, bins=np.arange(1, 11))
obs_probs = obs_counts.astype(float) / num_obs

# Perform Chi-square goodness-of-fit test
expected_counts = benford_probs * num_obs
chi_stat, p_value = chisquare(obs_counts, expected_counts)

print(f"\nChi-square statistic: {chi_stat:.4f}")
print(f"P-value: {p_value:.4f}")
if p_value > 0.05:
    print("Conclusion: The leading digits follow Benford's Law (fail to reject null hypothesis at 5% significance).")
else:
    print("Conclusion: The leading digits do not follow Benford's Law (reject null hypothesis at 5% significance).")

# THRML: Model Benford as a categorical EBM and sample to estimate frequencies
n_cats = 9  # States 0-8 corresponding to digits 1-9
# Biases for EBM: energy[k] = -log(p[k]) so that p(k) \propto exp(-energy[k])
biases = -jnp.log(jnp.array(benford_probs))

# Map observed leading digits to states (1->0, 2->1, ..., 9->8)
observed_states = observed - 1

# Define a single categorical node
node = CategoricalNode(num_states=n_cats)

# Create unary EBM with biases (no interactions)
model = CategoricalEBM([node], biases=biases)

# Define blocks: single free block for the node
free_blocks = [Block([node])]
spec = BlockGibbsSpec(free_blocks, [])  # No clamped blocks

# Sampling program for categorical EBM
program = CategoricalSamplingProgram(model, free_blocks, [])

# Sampling parameters
n_sim_samples = 10000  # Number of samples from Benford model (adjust as needed)
steps_per_sample = 1  # Sufficient for independent categorical samples
schedule = SamplingSchedule(
    n_warmup=0,
    n_samples=n_sim_samples,
    steps_per_sample=steps_per_sample
)

# Initialize state and keys
key = jax.random.PRNGKey(42)
key, subkey = jax.split(key)
init_state = jax.random.randint(subkey, (1,), minval=0, maxval=n_cats, dtype=jnp.uint8)  # Single batch

keys = jax.random.split(key, 1)

# Sample states (vmap over batches if scaling up)
def sample_batch(init, k):
    return sample_states(k, program, schedule, init, [], free_blocks)

samples = jax.jit(jax.vmap(sample_batch))(init_state[None, :], keys)  # Reshape for vmap

# Extract samples (shape: (1, n_samples, 1) -> flatten)
sim_samples = samples[0, :, 0]  # For single batch and single node

# Compute simulated frequencies
sim_samples_np = np.asarray(sim_samples)
sim_counts, _ = np.histogram(sim_samples_np + 1, bins=np.arange(1, 11))  # Back to 1-9
sim_probs = sim_counts / n_sim_samples
print(f"\nSimulated frequencies from THRML Benford EBM (top 3): {dict(zip(range(1,10), sim_probs[:9]))}")

# Optional: Chi-square using simulated expected (should match theoretical closely)
sim_expected_counts = sim_probs * num_obs
chi_sim, p_sim = chisquare(obs_counts, sim_expected_counts)
print(f"Chi-square with THRML-simulated expected: {chi_sim:.4f}, P-value: {p_sim:.4f}")

```
