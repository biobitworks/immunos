---
source: /Users/byron/projects/experiments/benford_weights.py
relative: experiments/benford_weights.py
generated_at: 2025-12-23 10:28
---

```python
import numpy as np
import jax
import jax.numpy as jnp
import time
from thrml import CategoricalNode
from thrml.block_management import Block
from thrml.block_sampling import BlockGibbsSpec, sample_states, SamplingSchedule
from thrml.models.discrete_ebm import CategoricalEBM
from thrml.sampling_programs import CategoricalSamplingProgram

# Constants
MAX_EPOCHS = 10_000_000
NUM_WEIGHTS = 1_000_000
LEARNING_RATE = 0.1
TEMPERATURE = 2.0
ENTROPY_WEIGHT = 0.05
EPS = 1e-8
INDIVIDUAL_LOSS_THRESHOLD = 0.005

# Benford distribution
def benford_distribution():
    return jnp.log10(1.0 + 1.0 / jnp.arange(1, 10))

# Softmax with temperature
def softmax(logits, T):
    logits = jnp.array(logits)
    max_logit = jnp.max(logits)
    exp_logits = jnp.exp((logits - max_logit) / T)
    return exp_logits / jnp.sum(exp_logits)

# KL divergence
def kl_divergence(p, q):
    p, q = jnp.array(p), jnp.array(q)
    return jnp.sum(p * jnp.log((p + EPS) / (q + EPS)))

# Entropy
def entropy(p):
    p = jnp.array(p)
    return -jnp.sum(p * jnp.log(p + EPS))

# Leading digit distribution
def get_digit_distribution(weights):
    digit_counts = jnp.zeros(9)
    for w in weights:
        if w <= 0:
            continue
        scaled = w / 10.0 ** jnp.floor(jnp.log10(w))
        first_digit = int(scaled)
        if 1 <= first_digit <= 9:
            digit_counts = digit_counts.at[first_digit - 1].add(1)
    total = jnp.sum(digit_counts)
    return digit_counts / total if total > 0 else jnp.zeros(9)

# THRML setup
def initialize_thrml_model():
    node = CategoricalNode(num_states=9)  # States 0-8 for digits 1-9
    biases = jnp.zeros(9)  # Initialize logits as biases
    return node, biases

def sample_weights(node, biases, n_samples, key):
    model = CategoricalEBM([node], biases=biases)
    free_blocks = [Block([node])]
    spec = BlockGibbsSpec(free_blocks, [])
    program = CategoricalSamplingProgram(model, free_blocks, [])
    schedule = SamplingSchedule(n_warmup=0, n_samples=n_samples, steps_per_sample=1)
    init_state = jnp.zeros((1, 1), dtype=jnp.uint8)
    keys = jax.random.split(key, 1)
    samples = sample_states(keys[0], program, schedule, init_state, [], free_blocks)
    # Map samples (0-8) to digits (1-9) and generate weights
    digits = samples[0, :, 0] + 1  # Shape: (n_samples,)
    rng = np.random.default_rng()
    magnitudes = rng.uniform(0.1, 100.0, n_samples)
    weights = digits * magnitudes  # Approximate weights with correct leading digits
    return weights

def main():
    start_time = time.time()

    # Initialize THRML model
    node, biases = initialize_thrml_model()
    target = benford_distribution()
    loss = 0.0
    key = jax.random.PRNGKey(42)

    for epoch in range(MAX_EPOCHS):
        # Sample weights
        key, subkey = jax.split(key)
        weights = sample_weights(node, biases, NUM_WEIGHTS, subkey)

        # Compute observed distribution
        observed_dist = get_digit_distribution(weights)

        # Compute loss
        probs = softmax(biases, TEMPERATURE)
        kl = kl_divergence(target, probs)
        H = entropy(probs)
        loss = kl - ENTROPY_WEIGHT * H

        # Check individual losses
        individual_losses = jnp.abs(observed_dist - target)
        all_below_threshold = jnp.all(individual_losses < INDIVIDUAL_LOSS_THRESHOLD)

        # Update biases (logits)
        if not all_below_threshold:
            grad_kl = observed_dist - target
            grad_entropy = -jnp.log(probs + EPS) - 1.0
            grad = grad_kl - ENTROPY_WEIGHT * grad_entropy
            biases = biases - LEARNING_RATE * grad

        if all_below_threshold:
            print(f"Early stopping at epoch {epoch} with KL loss {loss:.4f}")
            break

        if epoch % 1000 == 0:
            print(f"Epoch {epoch}, KL loss: {loss:.4f}")

    # Final output
    final_probs = softmax(biases, TEMPERATURE)
    final_dist = get_digit_distribution(weights)
    print("Final leading digit distribution:")
    for i in range(9):
        print(f"Digit {i + 1}: {final_dist[i]:.4f} (expected: {target[i]:.4f})")
    print(f"Final KL loss: {loss:.4f}")
    print(f"Execution time: {time.time() - start_time:.4f} seconds")

if __name__ == "__main__":
    main()

```
