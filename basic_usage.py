"""
Basic Usage Example for EnzymeKinetics-Skill

This example demonstrates the basic workflow for analyzing enzyme kinetic data
using the EnzymeKineticsAnalyzer class.
"""

import numpy as np
import sys
sys.path.insert(0, '../src')

from enzymekinetics import EnzymeKineticsAnalyzer

# Example 1: Glucose Oxidase Data
print("=" * 60)
print("Example 1: Glucose Oxidase Kinetic Analysis")
print("=" * 60)

# Experimental data (substrate concentrations in mM, velocities in μmol/min/mg)
substrate = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0])
velocity = np.array([0.15, 0.28, 0.52, 0.89, 1.35, 1.58, 1.68, 1.72, 1.75, 1.78])

# Initialize analyzer
analyzer = EnzymeKineticsAnalyzer(confidence_level=0.95)

# Fit Michaelis-Menten model
results = analyzer.fit_michaelis_menten(substrate, velocity)

# Print results
print("\nFitting Results:")
print("-" * 40)
print(f"Vmax = {results['Vmax']:.3f} ± {results['Vmax_std']:.3f} μmol/min/mg")
print(f"Km = {results['Km']:.3f} ± {results['Km_std']:.3f} mM")
print(f"R² = {results['r_squared']:.4f}")

# Perform bootstrap analysis
print("\nBootstrap Analysis (n=1000):")
print("-" * 40)
bootstrap_results = analyzer.bootstrap_analysis(substrate, velocity, n_bootstrap=1000)
print(f"Vmax 95% CI: [{bootstrap_results['Vmax_ci_lower']:.3f}, {bootstrap_results['Vmax_ci_upper']:.3f}]")
print(f"Km 95% CI: [{bootstrap_results['Km_ci_lower']:.3f}, {bootstrap_results['Km_ci_upper']:.3f}]")
print(f"Convergence rate: {bootstrap_results['convergence_rate']:.1f}%")

# Detect outliers
print("\nOutlier Detection:")
print("-" * 40)
outliers = analyzer.detect_outliers(substrate, velocity, 
                                     (results['Vmax'], results['Km']))
if np.any(outliers):
    print(f"Detected {np.sum(outliers)} outlier(s) at indices: {np.where(outliers)[0]}")
else:
    print("No outliers detected (all |r| < 3)")

# Compare methods
print("\nMethod Comparison:")
print("-" * 40)
comparison = analyzer.compare_methods(substrate, velocity)
for method, params in comparison.items():
    print(f"{method:20s}: Vmax={params['Vmax']:.3f}, Km={params['Km']:.3f}")

# Generate figure
print("\nGenerating figure...")
analyzer.plot_michaelis_menten(substrate, velocity, results, 
                                save_path='../figures/example1_michaelis_menten.png')
print("Figure saved to: figures/example1_michaelis_menten.png")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
