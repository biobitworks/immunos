---
source: /Users/byron/projects/immunos81/examples/dashboard_launcher.py
relative: immunos81/examples/dashboard_launcher.py
generated_at: 2025-12-23 10:28
---

```python
"""
Interactive Dashboard Launcher

Launches the Immunos-81 web dashboard for real-time visualization
and analysis.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def main():
    print("="*70)
    print("  IMMUNOS-81: Interactive Web Dashboard")
    print("="*70)
    print()
    print("Starting dashboard server...")
    print()
    print("Dashboard features:")
    print("  ✓ Train models with custom parameters")
    print("  ✓ Compare SHA vs RHA strategies")
    print("  ✓ Interactive confusion matrices")
    print("  ✓ Decision confidence analysis")
    print("  ✓ Variable importance ranking")
    print("  ✓ Performance benchmarking")
    print()
    print("-"*70)
    print()
    print("🌐 Dashboard URL: http://127.0.0.1:8050")
    print()
    print("Press Ctrl+C to stop the server")
    print()
    print("="*70)
    print()

    # Import and run dashboard
    from immunos81.dashboard.app import app
    app.run_server(debug=False, host='127.0.0.1', port=8050)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n" + "="*70)
        print("Dashboard stopped.")
        print("="*70)

```
