"""
Unit tests for EnzymeKinetics-Skill

Run with: pytest tests/test_enzymekinetics.py -v
"""

import unittest
import numpy as np
import sys
sys.path.insert(0, '../src')

from enzymekinetics import EnzymeKineticsAnalyzer


class TestEnzymeKineticsAnalyzer(unittest.TestCase):
    """Test cases for EnzymeKineticsAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = EnzymeKineticsAnalyzer(confidence_level=0.95)
        
        # Generate synthetic data with known parameters
        self.true_vmax = 1.0
        self.true_km = 0.5
        self.substrate = np.array([0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0])
        self.velocity = (self.true_vmax * self.substrate) / (self.true_km + self.substrate)
        # Add small noise
        np.random.seed(42)
        self.velocity_noisy = self.velocity + np.random.normal(0, 0.02, len(self.velocity))
    
    def test_michaelis_menten_equation(self):
        """Test Michaelis-Menten equation calculation."""
        S = np.array([0.5, 1.0, 2.0])
        Vmax, Km = 1.0, 0.5
        expected = np.array([0.5, 0.667, 0.8])
        result = self.analyzer.michaelis_menten(S, Vmax, Km)
        np.testing.assert_array_almost_equal(result, expected, decimal=3)
    
    def test_fit_michaelis_menten(self):
        """Test Michaelis-Menten fitting."""
        results = self.analyzer.fit_michaelis_menten(self.substrate, self.velocity_noisy)
        
        # Check that fitted parameters are close to true values
        self.assertAlmostEqual(results['Vmax'], self.true_vmax, delta=0.1)
        self.assertAlmostEqual(results['Km'], self.true_km, delta=0.1)
        self.assertGreater(results['r_squared'], 0.99)
    
    def test_bootstrap_analysis(self):
        """Test bootstrap confidence interval estimation."""
        results = self.analyzer.bootstrap_analysis(self.substrate, self.velocity_noisy, 
                                                    n_bootstrap=100)
        
        # Check that confidence intervals contain true values
        self.assertLess(results['Vmax_ci_lower'], self.true_vmax)
        self.assertGreater(results['Vmax_ci_upper'], self.true_vmax)
        self.assertLess(results['Km_ci_lower'], self.true_km)
        self.assertGreater(results['Km_ci_upper'], self.true_km)
        
        # Check convergence rate
        self.assertGreater(results['convergence_rate'], 90)
    
    def test_outlier_detection(self):
        """Test outlier detection methods."""
        # Create data with an outlier
        velocity_with_outlier = self.velocity_noisy.copy()
        velocity_with_outlier[3] = velocity_with_outlier[3] * 2  # Add outlier
        
        results = self.analyzer.fit_michaelis_menten(self.substrate, velocity_with_outlier)
        outliers = self.analyzer.detect_outliers(self.substrate, velocity_with_outlier,
                                                  (results['Vmax'], results['Km']))
        
        # Check that outlier is detected
        self.assertTrue(outliers[3])
    
    def test_compare_methods(self):
        """Test method comparison."""
        comparison = self.analyzer.compare_methods(self.substrate, self.velocity_noisy)
        
        # Check that all methods are present
        expected_methods = ['Nonlinear', 'Lineweaver-Burk', 'Eadie-Hofstee', 'Hanes-Woolf']
        for method in expected_methods:
            self.assertIn(method, comparison)
        
        # Check that all methods return valid parameters
        for method, params in comparison.items():
            self.assertGreater(params['Vmax'], 0)
            self.assertGreater(params['Km'], 0)
    
    def test_weighted_regression(self):
        """Test weighted regression."""
        results_unweighted = self.analyzer.fit_michaelis_menten(
            self.substrate, self.velocity_noisy, weights=None)
        results_weighted = self.analyzer.fit_michaelis_menten(
            self.substrate, self.velocity_noisy, weights='1/v')
        
        # Both should produce reasonable results
        self.assertGreater(results_unweighted['r_squared'], 0.95)
        self.assertGreater(results_weighted['r_squared'], 0.95)


class TestDataValidation(unittest.TestCase):
    """Test data validation functions."""
    
    def setUp(self):
        self.analyzer = EnzymeKineticsAnalyzer()
    
    def test_invalid_input_lengths(self):
        """Test that mismatched input lengths raise error."""
        substrate = np.array([1, 2, 3])
        velocity = np.array([0.1, 0.2])
        
        with self.assertRaises(ValueError):
            self.analyzer.fit_michaelis_menten(substrate, velocity)
    
    def test_negative_values(self):
        """Test that negative values are handled appropriately."""
        substrate = np.array([-1, 2, 3])
        velocity = np.array([0.1, 0.2, 0.3])
        
        with self.assertRaises(ValueError):
            self.analyzer.fit_michaelis_menten(substrate, velocity)


if __name__ == '__main__':
    unittest.main()
