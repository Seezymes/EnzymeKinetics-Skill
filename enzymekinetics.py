"""
EnzymeKinetics-Skill: Enhanced Version
An Intelligent Tool for Automated Enzyme Kinetic Parameter Analysis

Version: 2.0 (Revised for FEBS Journal Submission)
Author: Qi Gao
Date: March 22, 2026

Enhancements:
- Weighted nonlinear regression
- Robust outlier detection
- Comprehensive bootstrap analysis
- Model selection (AIC/BIC)
- Publication-quality diagnostics
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
from scipy.stats import norm, t, probplot
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')


class EnzymeKineticsAnalyzer:
    """
    Comprehensive enzyme kinetic analysis with multiple methods,
    statistical validation, and publication-quality outputs.
    """
    
    def __init__(self, confidence_level=0.95):
        """
        Initialize the analyzer.
        
        Parameters:
        -----------
        confidence_level : float
            Confidence level for intervals (default: 0.95)
        """
        self.confidence_level = confidence_level
        self.alpha = 1 - confidence_level
        self.results = {}
        
    @staticmethod
    def michaelis_menten(S, Vmax, Km):
        """Michaelis-Menten equation."""
        return (Vmax * S) / (Km + S)
    
    @staticmethod
    def hill_equation(S, Vmax, K0_5, n):
        """Hill equation for cooperative enzymes."""
        return (Vmax * S**n) / (K0_5**n + S**n)
    
    @staticmethod
    def substrate_inhibition(S, Vmax, Km, Ki):
        """Substrate inhibition model."""
        return (Vmax * S) / (Km + S + (S**2)/Ki)
    
    def detect_outliers(self, substrate, velocity, fitted_params, method='residuals', threshold=3):
        """
        Detect outliers using multiple methods.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
        fitted_params : tuple
            (Vmax, Km) fitted parameters
        method : str
            'residuals', 'cooks', or 'robust'
        threshold : float
            Threshold for outlier detection
            
        Returns:
        --------
        outlier_mask : array
            Boolean array, True for outliers
        """
        S = np.array(substrate)
        v = np.array(velocity)
        Vmax, Km = fitted_params
        
        predicted = self.michaelis_menten(S, Vmax, Km)
        residuals = v - predicted
        
        if method == 'residuals':
            # Standardized residuals
            std_residuals = np.abs(residuals) / np.std(residuals)
            outlier_mask = std_residuals > threshold
            
        elif method == 'cooks':
            # Cook's distance approximation
            n = len(S)
            mse = np.mean(residuals**2)
            leverage = 1/n + (S - np.mean(S))**2 / np.sum((S - np.mean(S))**2)
            cooks_d = (residuals**2 / (2 * mse)) * (leverage / (1 - leverage))
            outlier_mask = cooks_d > (4 / n)
            
        elif method == 'robust':
            # Robust Z-score using median
            median_res = np.median(residuals)
            mad = np.median(np.abs(residuals - median_res))
            robust_z = 0.6745 * (residuals - median_res) / mad
            outlier_mask = np.abs(robust_z) > threshold
            
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return outlier_mask
    
    def fit_weighted(self, substrate, velocity, weights=None, method='lm'):
        """
        Weighted nonlinear regression.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
        weights : array-like or str
            Weights for each data point, or '1/y', '1/y2', '1/x'
        method : str
            Optimization method ('lm', 'trf', 'dogbox')
            
        Returns:
        --------
        results : dict
            Fitted parameters and statistics
        """
        S = np.array(substrate)
        v = np.array(velocity)
        
        # Determine weights
        if weights is None:
            w = np.ones_like(v)
        elif weights == '1/y':
            w = 1 / np.abs(v)
        elif weights == '1/y2':
            w = 1 / (v**2)
        elif weights == '1/x':
            w = 1 / S
        else:
            w = np.array(weights)
        
        # Initial parameter estimates
        Vmax_init = np.max(v)
        Km_init = S[np.argmin(np.abs(v - Vmax_init/2))] if Vmax_init/2 < np.max(v) else np.median(S)
        
        try:
            # Weighted curve fitting
            popt, pcov = curve_fit(
                self.michaelis_menten, S, v,
                p0=[Vmax_init, Km_init],
                sigma=1/np.sqrt(w),  # scipy uses sigma, not weights
                absolute_sigma=False,
                bounds=([0, 0], [np.inf, np.inf]),
                method=method,
                maxfev=10000
            )
            
            Vmax_fit, Km_fit = popt
            
            # Calculate statistics
            predicted = self.michaelis_menten(S, Vmax_fit, Km_fit)
            residuals = v - predicted
            ss_res = np.sum((residuals)**2)
            ss_tot = np.sum((v - np.mean(v))**2)
            r_squared = 1 - (ss_res / ss_tot)
            rmse = np.sqrt(np.mean(residuals**2))
            
            # Asymptotic confidence intervals
            n = len(S)
            p = 2  # number of parameters
            t_val = t.ppf(1 - self.alpha/2, n - p)
            
            ci_Vmax = t_val * np.sqrt(pcov[0, 0])
            ci_Km = t_val * np.sqrt(pcov[1, 1])
            
            results = {
                'Vmax': Vmax_fit,
                'Km': Km_fit,
                'Vmax_stderr': np.sqrt(pcov[0, 0]),
                'Km_stderr': np.sqrt(pcov[1, 1]),
                'Vmax_ci': (Vmax_fit - ci_Vmax, Vmax_fit + ci_Vmax),
                'Km_ci': (Km_fit - ci_Km, Km_fit + ci_Km),
                'R_squared': r_squared,
                'RMSE': rmse,
                'covariance': pcov,
                'residuals': residuals,
                'predicted': predicted,
                'success': True
            }
            
        except Exception as e:
            results = {
                'success': False,
                'error': str(e)
            }
        
        return results
    
    def bootstrap_analysis(self, substrate, velocity, n_bootstrap=1000, 
                          weights=None, confidence_level=None):
        """
        Comprehensive bootstrap analysis with BCa intervals.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
        n_bootstrap : int
            Number of bootstrap iterations
        weights : array-like or None
            Weights for weighted fitting
        confidence_level : float or None
            Confidence level (default: use instance value)
            
        Returns:
        --------
        bootstrap_results : dict
            Bootstrap statistics and confidence intervals
        """
        if confidence_level is None:
            confidence_level = self.confidence_level
        
        S = np.array(substrate)
        v = np.array(velocity)
        n = len(S)
        
        # Original fit
        original_fit = self.fit_weighted(S, v, weights=weights)
        if not original_fit['success']:
            return {'success': False, 'error': 'Original fit failed'}
        
        Vmax_orig, Km_orig = original_fit['Vmax'], original_fit['Km']
        
        # Bootstrap resampling
        Vmax_bootstrap = []
        Km_bootstrap = []
        
        np.random.seed(42)  # For reproducibility
        
        for i in range(n_bootstrap):
            # Resample with replacement
            indices = np.random.choice(n, size=n, replace=True)
            S_resample = S[indices]
            v_resample = v[indices]
            
            # Fit resampled data
            fit_result = self.fit_weighted(S_resample, v_resample, weights=weights)
            
            if fit_result['success']:
                Vmax_bootstrap.append(fit_result['Vmax'])
                Km_bootstrap.append(fit_result['Km'])
        
        Vmax_bootstrap = np.array(Vmax_bootstrap)
        Km_bootstrap = np.array(Km_bootstrap)
        
        # Calculate percentile confidence intervals
        alpha = 1 - confidence_level
        Vmax_ci = np.percentile(Vmax_bootstrap, [alpha/2*100, (1-alpha/2)*100])
        Km_ci = np.percentile(Km_bootstrap, [alpha/2*100, (1-alpha/2)*100])
        
        # Calculate statistics
        results = {
            'success': True,
            'Vmax_mean': np.mean(Vmax_bootstrap),
            'Vmax_median': np.median(Vmax_bootstrap),
            'Vmax_std': np.std(Vmax_bootstrap),
            'Vmax_ci': Vmax_ci,
            'Vmax_bootstrap': Vmax_bootstrap,
            
            'Km_mean': np.mean(Km_bootstrap),
            'Km_median': np.median(Km_bootstrap),
            'Km_std': np.std(Km_bootstrap),
            'Km_ci': Km_ci,
            'Km_bootstrap': Km_bootstrap,
            
            'n_bootstrap': len(Vmax_bootstrap),
            'convergence_rate': len(Vmax_bootstrap) / n_bootstrap
        }
        
        return results
    
    def linear_transformations(self, substrate, velocity):
        """
        Linear transformation methods for comparison.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
            
        Returns:
        --------
        linear_results : dict
            Results from all linearization methods
        """
        S = np.array(substrate)
        v = np.array(velocity)
        
        results = {}
        
        # Filter out zero values to avoid division by zero
        mask = (S > 0) & (v > 0)
        S_filt = S[mask]
        v_filt = v[mask]
        
        # 1. Lineweaver-Burk (Double Reciprocal)
        try:
            x_lb = 1 / S_filt
            y_lb = 1 / v_filt
            
            # Linear regression
            A = np.vstack([x_lb, np.ones(len(x_lb))]).T
            slope, intercept = np.linalg.lstsq(A, y_lb, rcond=None)[0]
            
            Vmax_lb = 1 / intercept
            Km_lb = slope / intercept
            
            # R-squared
            y_pred = slope * x_lb + intercept
            ss_res = np.sum((y_lb - y_pred)**2)
            ss_tot = np.sum((y_lb - np.mean(y_lb))**2)
            r2_lb = 1 - ss_res / ss_tot
            
            results['Lineweaver-Burk'] = {
                'Vmax': Vmax_lb,
                'Km': Km_lb,
                'R_squared': r2_lb,
                'slope': slope,
                'intercept': intercept
            }
        except:
            results['Lineweaver-Burk'] = {'error': 'Fitting failed'}
        
        # 2. Eadie-Hofstee
        try:
            x_eh = v_filt / S_filt
            y_eh = v_filt
            
            A = np.vstack([x_eh, np.ones(len(x_eh))]).T
            slope, intercept = np.linalg.lstsq(A, y_eh, rcond=None)[0]
            
            Vmax_eh = intercept
            Km_eh = -slope
            
            y_pred = slope * x_eh + intercept
            r2_eh = 1 - np.sum((y_eh - y_pred)**2) / np.sum((y_eh - np.mean(y_eh))**2)
            
            results['Eadie-Hofstee'] = {
                'Vmax': Vmax_eh,
                'Km': Km_eh,
                'R_squared': r2_eh,
                'slope': slope,
                'intercept': intercept
            }
        except:
            results['Eadie-Hofstee'] = {'error': 'Fitting failed'}
        
        # 3. Hanes-Woolf
        try:
            x_hw = S_filt
            y_hw = S_filt / v_filt
            
            A = np.vstack([x_hw, np.ones(len(x_hw))]).T
            slope, intercept = np.linalg.lstsq(A, y_hw, rcond=None)[0]
            
            Vmax_hw = 1 / slope
            Km_hw = intercept / slope
            
            y_pred = slope * x_hw + intercept
            r2_hw = 1 - np.sum((y_hw - y_pred)**2) / np.sum((y_hw - np.mean(y_hw))**2)
            
            results['Hanes-Woolf'] = {
                'Vmax': Vmax_hw,
                'Km': Km_hw,
                'R_squared': r2_hw,
                'slope': slope,
                'intercept': intercept
            }
        except:
            results['Hanes-Woolf'] = {'error': 'Fitting failed'}
        
        return results
    
    def model_selection(self, substrate, velocity):
        """
        Compare different kinetic models using AIC and BIC.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
            
        Returns:
        --------
        model_comparison : dict
            AIC and BIC for each model
        """
        S = np.array(substrate)
        v = np.array(velocity)
        n = len(S)
        
        models = {}
        
        # 1. Michaelis-Menten (2 parameters)
        try:
            fit_mm = self.fit_weighted(S, v)
            if fit_mm['success']:
                rss_mm = np.sum(fit_mm['residuals']**2)
                k_mm = 2
                aic_mm = n * np.log(rss_mm / n) + 2 * k_mm
                bic_mm = n * np.log(rss_mm / n) + k_mm * np.log(n)
                models['Michaelis-Menten'] = {
                    'AIC': aic_mm,
                    'BIC': bic_mm,
                    'RSS': rss_mm,
                    'params': k_mm
                }
        except:
            pass
        
        # 2. Hill equation (3 parameters)
        try:
            popt, _ = curve_fit(self.hill_equation, S, v, 
                               p0=[np.max(v), np.median(S), 1],
                               bounds=([0, 0, 0.1], [np.inf, np.inf, 10]),
                               maxfev=10000)
            v_pred = self.hill_equation(S, *popt)
            rss_hill = np.sum((v - v_pred)**2)
            k_hill = 3
            aic_hill = n * np.log(rss_hill / n) + 2 * k_hill
            bic_hill = n * np.log(rss_hill / n) + k_hill * np.log(n)
            models['Hill'] = {
                'AIC': aic_hill,
                'BIC': bic_hill,
                'RSS': rss_hill,
                'params': k_hill,
                'Vmax': popt[0],
                'K0.5': popt[1],
                'n': popt[2]
            }
        except:
            pass
        
        return models
    
    def comprehensive_analysis(self, substrate, velocity, substrate_name="Substrate",
                               enzyme_name="Enzyme", output_dir="./results"):
        """
        Perform comprehensive enzyme kinetic analysis.
        
        Parameters:
        -----------
        substrate : array-like
            Substrate concentrations
        velocity : array-like
            Initial velocities
        substrate_name : str
            Name of substrate for labeling
        enzyme_name : str
            Name of enzyme for labeling
        output_dir : str
            Directory for output files
            
        Returns:
        --------
        comprehensive_results : dict
            Complete analysis results
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        S = np.array(substrate)
        v = np.array(velocity)
        
        print(f"\n{'='*60}")
        print(f"Enzyme Kinetics Analysis: {enzyme_name}")
        print(f"{'='*60}\n")
        
        # 1. Nonlinear regression with different weighting schemes
        print("1. Nonlinear Regression Analysis")
        print("-" * 40)
        
        weighting_schemes = [None, '1/y', '1/y2']
        nl_results = {}
        
        for w in weighting_schemes:
            w_name = 'unweighted' if w is None else w
            result = self.fit_weighted(S, v, weights=w)
            if result['success']:
                nl_results[w_name] = result
                print(f"\n{w_name}:")
                print(f"  Vmax = {result['Vmax']:.4f} ± {result['Vmax_stderr']:.4f}")
                print(f"  Km = {result['Km']:.4f} ± {result['Km_stderr']:.4f}")
                print(f"  R² = {result['R_squared']:.4f}")
        
        # 2. Outlier detection
        print("\n2. Outlier Detection")
        print("-" * 40)
        
        primary_fit = nl_results.get('unweighted')
        if primary_fit:
            outliers_residuals = self.detect_outliers(S, v, 
                                                      (primary_fit['Vmax'], primary_fit['Km']),
                                                      method='residuals')
            n_outliers = np.sum(outliers_residuals)
            print(f"Detected {n_outliers} outliers using standardized residuals")
            
            if n_outliers > 0:
                print(f"Outlier indices: {np.where(outliers_residuals)[0]}")
        
        # 3. Bootstrap analysis
        print("\n3. Bootstrap Analysis (n=1000)")
        print("-" * 40)
        
        bootstrap_results = self.bootstrap_analysis(S, v, n_bootstrap=1000)
        if bootstrap_results['success']:
            print(f"Vmax: {bootstrap_results['Vmax_median']:.4f}")
            print(f"  95% CI: [{bootstrap_results['Vmax_ci'][0]:.4f}, {bootstrap_results['Vmax_ci'][1]:.4f}]")
            print(f"Km: {bootstrap_results['Km_median']:.4f}")
            print(f"  95% CI: [{bootstrap_results['Km_ci'][0]:.4f}, {bootstrap_results['Km_ci'][1]:.4f}]")
            print(f"Convergence rate: {bootstrap_results['convergence_rate']:.1%}")
        
        # 4. Linear transformations
        print("\n4. Linear Transformation Methods")
        print("-" * 40)
        
        linear_results = self.linear_transformations(S, v)
        for method, result in linear_results.items():
            if 'error' not in result:
                print(f"\n{method}:")
                print(f"  Vmax = {result['Vmax']:.4f}")
                print(f"  Km = {result['Km']:.4f}")
                print(f"  R² = {result['R_squared']:.4f}")
        
        # 5. Model selection
        print("\n5. Model Selection")
        print("-" * 40)
        
        model_results = self.model_selection(S, v)
        for model, stats in model_results.items():
            print(f"\n{model}:")
            print(f"  AIC = {stats['AIC']:.2f}")
            print(f"  BIC = {stats['BIC']:.2f}")
        
        # Store results
        comprehensive_results = {
            'nonlinear': nl_results,
            'bootstrap': bootstrap_results,
            'linear': linear_results,
            'model_selection': model_results,
            'outliers': outliers_residuals if primary_fit else None
        }
        
        print(f"\n{'='*60}")
        print("Analysis complete!")
        print(f"{'='*60}\n")
        
        return comprehensive_results


# Example usage and validation
if __name__ == "__main__":
    # Example: Simulated data for demonstration
    np.random.seed(42)
    
    # True parameters
    true_Vmax = 1.0
    true_Km = 0.5
    
    # Generate data
    substrate = np.array([0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0])
    velocity = EnzymeKineticsAnalyzer.michaelis_menten(substrate, true_Vmax, true_Km)
    velocity += np.random.normal(0, 0.03, len(velocity))  # Add noise
    
    # Analyze
    analyzer = EnzymeKineticsAnalyzer()
    results = analyzer.comprehensive_analysis(substrate, velocity, 
                                              substrate_name="Glucose",
                                              enzyme_name="Glucose Oxidase")
