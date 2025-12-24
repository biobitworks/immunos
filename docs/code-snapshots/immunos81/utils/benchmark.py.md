---
source: /Users/byron/projects/immunos81/utils/benchmark.py
relative: immunos81/utils/benchmark.py
generated_at: 2025-12-23 10:28
---

```python
"""
Performance Benchmarking Tools for Immunos-81

Compare sequential vs parallel implementations to measure speedup.
"""

import time
from typing import List, Dict, Tuple
from ..core.antigen import Antigen
from ..engine.trainer import Trainer
from ..engine.parallel_trainer import ParallelTrainer
from ..engine.recognizer import Recognizer, RecognitionStrategy
from ..engine.parallel_recognizer import ParallelRecognizer


class Benchmark:
    """
    Performance benchmarking for Immunos-81.

    Compares sequential vs parallel implementations to measure:
    - Training time reduction
    - Recognition time reduction
    - Overall speedup factor
    """

    @staticmethod
    def benchmark_training(train_data: List[Antigen],
                          n_workers: int = None,
                          n_runs: int = 3) -> Dict[str, any]:
        """
        Benchmark training performance.

        Args:
            train_data: Training dataset
            n_workers: Number of workers for parallel version
            n_runs: Number of runs to average over

        Returns:
            Dictionary with benchmark results
        """
        print(f"Benchmarking Training Performance")
        print(f"Dataset size: {len(train_data)} instances")
        print(f"Number of runs: {n_runs}\n")

        library = train_data[0].library

        # Sequential training
        print("Running sequential training...")
        seq_times = []
        for run in range(n_runs):
            trainer = Trainer(library)
            start = time.time()
            trainer.train(train_data)
            elapsed = time.time() - start
            seq_times.append(elapsed)
            print(f"  Run {run+1}: {elapsed:.2f}s")

        seq_mean = sum(seq_times) / len(seq_times)
        print(f"Sequential mean: {seq_mean:.2f}s\n")

        # Parallel training
        print("Running parallel training...")
        par_times = []
        for run in range(n_runs):
            trainer = ParallelTrainer(library, n_workers=n_workers)
            start = time.time()
            trainer.train(train_data)
            elapsed = time.time() - start
            par_times.append(elapsed)
            print(f"  Run {run+1}: {elapsed:.2f}s")

        par_mean = sum(par_times) / len(par_times)
        print(f"Parallel mean: {par_mean:.2f}s\n")

        # Calculate speedup
        speedup = seq_mean / par_mean
        improvement = (1 - par_mean / seq_mean) * 100

        print(f"{'='*50}")
        print(f"TRAINING BENCHMARK RESULTS")
        print(f"{'='*50}")
        print(f"Sequential:  {seq_mean:.2f}s ± {max(seq_times) - min(seq_times):.2f}s")
        print(f"Parallel:    {par_mean:.2f}s ± {max(par_times) - min(par_times):.2f}s")
        print(f"Speedup:     {speedup:.2f}x")
        print(f"Improvement: {improvement:.1f}% faster")
        print(f"{'='*50}\n")

        return {
            "dataset_size": len(train_data),
            "n_runs": n_runs,
            "n_workers": n_workers,
            "sequential_times": seq_times,
            "parallel_times": par_times,
            "sequential_mean": seq_mean,
            "parallel_mean": par_mean,
            "speedup": speedup,
            "improvement_pct": improvement
        }

    @staticmethod
    def benchmark_recognition(train_data: List[Antigen],
                             test_data: List[Antigen],
                             n_workers: int = None,
                             n_runs: int = 3) -> Dict[str, any]:
        """
        Benchmark recognition performance.

        Args:
            train_data: Training dataset
            test_data: Test dataset
            n_workers: Number of workers for parallel version
            n_runs: Number of runs to average over

        Returns:
            Dictionary with benchmark results
        """
        print(f"Benchmarking Recognition Performance")
        print(f"Training size: {len(train_data)} instances")
        print(f"Test size: {len(test_data)} instances")
        print(f"Number of runs: {n_runs}\n")

        library = train_data[0].library

        # Train once (use parallel for speed)
        print("Training model...")
        trainer = ParallelTrainer(library, n_workers=n_workers)
        tcells = trainer.train(train_data)
        print()

        # Sequential recognition
        print("Running sequential recognition...")
        seq_times = []
        recognizer_seq = Recognizer(tcells, strategy=RecognitionStrategy.SHA)
        for run in range(n_runs):
            start = time.time()
            recognizer_seq.batch_recognize(test_data)
            elapsed = time.time() - start
            seq_times.append(elapsed)
            print(f"  Run {run+1}: {elapsed:.2f}s")

        seq_mean = sum(seq_times) / len(seq_times)
        print(f"Sequential mean: {seq_mean:.2f}s\n")

        # Parallel recognition
        print("Running parallel recognition...")
        par_times = []
        recognizer_par = ParallelRecognizer(tcells, strategy=RecognitionStrategy.SHA,
                                           n_workers=n_workers)
        for run in range(n_runs):
            start = time.time()
            recognizer_par.batch_recognize(test_data)
            elapsed = time.time() - start
            par_times.append(elapsed)
            print(f"  Run {run+1}: {elapsed:.2f}s")

        par_mean = sum(par_times) / len(par_times)
        print(f"Parallel mean: {par_mean:.2f}s\n")

        # Calculate speedup
        speedup = seq_mean / par_mean
        improvement = (1 - par_mean / seq_mean) * 100

        print(f"{'='*50}")
        print(f"RECOGNITION BENCHMARK RESULTS")
        print(f"{'='*50}")
        print(f"Sequential:  {seq_mean:.2f}s ± {max(seq_times) - min(seq_times):.2f}s")
        print(f"Parallel:    {par_mean:.2f}s ± {max(par_times) - min(par_times):.2f}s")
        print(f"Speedup:     {speedup:.2f}x")
        print(f"Improvement: {improvement:.1f}% faster")
        print(f"{'='*50}\n")

        return {
            "train_size": len(train_data),
            "test_size": len(test_data),
            "n_runs": n_runs,
            "n_workers": n_workers,
            "sequential_times": seq_times,
            "parallel_times": par_times,
            "sequential_mean": seq_mean,
            "parallel_mean": par_mean,
            "speedup": speedup,
            "improvement_pct": improvement
        }

    @staticmethod
    def benchmark_end_to_end(train_data: List[Antigen],
                            test_data: List[Antigen],
                            n_workers: int = None) -> Dict[str, any]:
        """
        Benchmark complete training + recognition pipeline.

        Args:
            train_data: Training dataset
            test_data: Test dataset
            n_workers: Number of workers for parallel version

        Returns:
            Dictionary with benchmark results
        """
        print(f"Benchmarking End-to-End Performance")
        print(f"Training size: {len(train_data)}, Test size: {len(test_data)}\n")

        library = train_data[0].library

        # Sequential pipeline
        print("Running sequential pipeline...")
        start_seq = time.time()
        trainer_seq = Trainer(library)
        tcells_seq = trainer_seq.train(train_data)
        recognizer_seq = Recognizer(tcells_seq, strategy=RecognitionStrategy.SHA)
        results_seq = recognizer_seq.batch_recognize(test_data)
        metrics_seq = recognizer_seq.evaluate(test_data)
        time_seq = time.time() - start_seq
        print(f"Sequential total: {time_seq:.2f}s")
        print(f"Accuracy: {metrics_seq['accuracy_overall']:.1%}\n")

        # Parallel pipeline
        print("Running parallel pipeline...")
        start_par = time.time()
        trainer_par = ParallelTrainer(library, n_workers=n_workers)
        tcells_par = trainer_par.train(train_data)
        recognizer_par = ParallelRecognizer(tcells_par, strategy=RecognitionStrategy.SHA,
                                           n_workers=n_workers)
        results_par = recognizer_par.batch_recognize(test_data)
        metrics_par = recognizer_par.evaluate(test_data)
        time_par = time.time() - start_par
        print(f"Parallel total: {time_par:.2f}s")
        print(f"Accuracy: {metrics_par['accuracy_overall']:.1%}\n")

        # Calculate speedup
        speedup = time_seq / time_par
        improvement = (1 - time_par / time_seq) * 100

        print(f"{'='*50}")
        print(f"END-TO-END BENCHMARK RESULTS")
        print(f"{'='*50}")
        print(f"Sequential:  {time_seq:.2f}s")
        print(f"Parallel:    {time_par:.2f}s")
        print(f"Speedup:     {speedup:.2f}x")
        print(f"Improvement: {improvement:.1f}% faster")
        print(f"Accuracy match: {abs(metrics_seq['accuracy_overall'] - metrics_par['accuracy_overall']) < 0.001}")
        print(f"{'='*50}\n")

        return {
            "train_size": len(train_data),
            "test_size": len(test_data),
            "n_workers": n_workers,
            "sequential_time": time_seq,
            "parallel_time": time_par,
            "speedup": speedup,
            "improvement_pct": improvement,
            "sequential_accuracy": metrics_seq['accuracy_overall'],
            "parallel_accuracy": metrics_par['accuracy_overall']
        }

```
