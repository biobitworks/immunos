---
source: /Users/byron/projects/experiments/benford_pi_simple.py
relative: experiments/benford_pi_simple.py
generated_at: 2025-12-23 10:28
---

```python
import sys
import numpy as np
from mpmath import mp
from scipy.stats import chisquare
import jax
import jax.numpy as jnp

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
print(f"\nGenerated {num_obs} leading digit observations from pi's decimal expansion.")

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

print(f"\nBenford's Law Analysis:")
print(f"Chi-square statistic: {chi_stat:.4f}")
print(f"P-value: {p_value:.4f}")
if p_value > 0.05:
    print("Conclusion: The leading digits follow Benford's Law (fail to reject null hypothesis at 5% significance).")
else:
    print("Conclusion: The leading digits do not follow Benford's Law (reject null hypothesis at 5% significance).")

# Display distribution comparison
print(f"\nDistribution Comparison (Leading Digit, Observed, Expected):")
for i in range(9):
    print(f"  {i+1}: {obs_probs[i]:.4f} vs {benford_probs[i]:.4f} (diff: {abs(obs_probs[i] - benford_probs[i]):.4f})")

# Simple JAX-based sampling from Benford distribution (multinomial)
print(f"\nSampling {10000} values from theoretical Benford distribution using JAX...")
key = jax.random.PRNGKey(42)
# Normalize benford_probs to sum to 1
benford_probs_normalized = benford_probs / benford_probs.sum()
samples = jax.random.choice(key, jnp.arange(1, 10), shape=(10000,), p=jnp.array(benford_probs_normalized))

# Compute simulated frequencies
sim_samples_np = np.asarray(samples)
sim_counts, _ = np.histogram(sim_samples_np, bins=np.arange(1, 11))
sim_probs = sim_counts / 10000

print(f"\nSimulated frequencies from Benford distribution:")
for i in range(9):
    print(f"  Digit {i+1}: {sim_probs[i]:.4f} (expected: {benford_probs[i]:.4f})")

# Chi-square with simulated expected
sim_expected_counts = sim_probs * num_obs
chi_sim, p_sim = chisquare(obs_counts, sim_expected_counts)
print(f"\nChi-square with JAX-simulated expected: {chi_sim:.4f}, P-value: {p_sim:.4f}")

```
