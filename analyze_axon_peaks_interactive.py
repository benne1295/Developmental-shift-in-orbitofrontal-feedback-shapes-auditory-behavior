"""
Enhanced Interactive Axon Peak Analysis
- Drag threshold line directly in graph
- Adjustable minimum distance between peaks
- Navigate freely between axons (Previous/Next)
- Auto-save progress and resume from last position
- Extract cortical layer (L1/L5) from filename
"""

# ============================================================================
# CONFIGURATION - CHANGE THESE PATHS FOR YOUR DATA
# ============================================================================

INPUT_DIRECTORY = '/Users/bennepraegel/Downloads/adult_filtered'
OUTPUT_DIRECTORY = '/Users/bennepraegel/Downloads/adult_filtered'

# Auto-run mode: Set to True to automatically process all axons with auto-detected parameters
# You can still review and edit them afterwards using the GUI
AUTO_RUN_ALL = True  # Change to False for manual mode

# Force re-process: Set to True to re-calculate parameters for ALL axons (ignores saved progress)
# Useful if you want to start fresh with new optimization
FORCE_REPROCESS = True  # Change to True to ignore previous progress and re-optimize all

# Output files will be named:
# - consolidated_axon_peaks_interactive.xlsx (main results)
# - consolidated_axon_peaks_interactive_summary_by_layer.xlsx (L1 vs L5 summary)
# - analysis_progress.json (auto-saved progress for resuming)

# ============================================================================

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.lines import Line2D
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from typing import Tuple, List, Dict, Optional
import glob


class InteractiveAxonAnalyzer:
    """Interactive peak detection with draggable threshold and navigation."""
    
    def __init__(self, all_data: List[Dict], progress_file: str):
        """
        Initialize analyzer with all axons from all files.
        
        Parameters:
        -----------
        all_data : list of dicts containing axon information
        progress_file : path to save/load progress
        """
        self.all_data = all_data
        self.progress_file = progress_file
        self.current_idx = 0
        self.results = {}
        
        # Load previous progress if exists
        self.load_progress()
        
        # GUI elements
        self.fig = None
        self.ax = None
        self.threshold_line = None
        self.data_line = None
        self.smoothed_line = None
        self.peak_scatter = None
        self.slider_min_dist = None
        self.slider_smoothing = None
        self.slider_intercept = None
        self.slider_slope = None
        self.btn_previous = None
        self.btn_next = None
        self.btn_skip = None
        self.btn_save_quit = None
        self.btn_auto_detect = None
        self.info_text = None
        
        # Current state
        self.dragging = False
        self.drag_mode = None  # 'intercept' or 'slope'
        self.loading_axon = False  # Flag to prevent slider callbacks during load
        
    def load_progress(self):
        """Load previous progress from file."""
        if os.path.exists(self.progress_file):
            try:
                with open(self.progress_file, 'r') as f:
                    data = json.load(f)
                    # Always start from axon 1 (index 0) when opening GUI
                    self.current_idx = 0
                    self.results = data.get('results', {})
                print(f"\n✓ Loaded previous progress: {len(self.results)} axons already analyzed")
            except Exception as e:
                print(f"Warning: Could not load progress file: {e}")
                self.current_idx = 0
                self.results = {}
        else:
            print("\nNo previous progress found. Starting from beginning.")
            self.current_idx = 0
            self.results = {}
    
    def save_progress(self):
        """Save current progress to file."""
        progress_data = {
            'current_idx': self.current_idx,
            'results': self.results
        }
        try:
            with open(self.progress_file, 'w') as f:
                json.dump(progress_data, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save progress: {e}")
    
    def get_current_axon(self) -> Dict:
        """Get current axon data."""
        return self.all_data[self.current_idx]
    
    def get_result_key(self, axon_data: Dict) -> str:
        """Generate unique key for axon results."""
        return f"{axon_data['file_name']}__axon_{axon_data['axon_num']}"
    
    def auto_detect_optimal_smoothing(self, distance_um: np.ndarray, intensity: np.ndarray) -> float:
        """
        Auto-detect optimal smoothing sigma for this axon.
        Tests different smoothing values and picks one that reduces noise while preserving peaks.
        """
        try:
            # Ensure data is numeric float type
            intensity = np.asarray(intensity, dtype=np.float64)
            distance_um = np.asarray(distance_um, dtype=np.float64)
            
            # Calculate signal-to-noise ratio for different smoothing values
            # Include a lighter smoothing (σ=1.0) and use it as default
            test_sigmas = [0.0, 1.0, 2.0, 5.0, 8.0, 10.0]
            best_sigma = 1.0
            best_score = -np.inf
            
            # Use a simple heuristic: we want to reduce noise but keep clear peaks
            original_std = np.std(intensity)
            
            for sigma in test_sigmas:
                if sigma > 0:
                    smoothed = gaussian_filter1d(intensity, sigma=sigma)
                else:
                    smoothed = intensity
                
                # Calculate metrics
                smoothed_std = np.std(smoothed)
                
                # Try to find peaks with this smoothing
                threshold = np.mean(smoothed) + 0.5 * np.std(smoothed)
                threshold_line = threshold * np.ones_like(smoothed)
                corrected = smoothed - threshold_line
                
                avg_spacing = np.mean(np.diff(distance_um))
                # Use a very small base distance (1.0 µm) when exploring raw peak distribution
                min_samples = max(1, int(1.0 / avg_spacing))
                
                try:
                    peak_indices, properties = find_peaks(
                        corrected,
                        height=0,
                        distance=min_samples,
                        prominence=0.3 * np.std(smoothed)
                    )
                    n_peaks = len(peak_indices)
                except:
                    n_peaks = 0
                
                # Score: prefer moderate smoothing that still finds reasonable number of peaks
                # Penalize too much smoothing (low std) and too little (high std, few peaks)
                if n_peaks > 0:
                    score = n_peaks * (1.0 - abs(smoothed_std - 0.5 * original_std) / original_std)
                else:
                    score = -1000
                
                if score > best_score:
                    best_score = score
                    best_sigma = sigma
            
            # Default to 1.0 if optimization fails
            return best_sigma if best_sigma > 0 else 1.0
        except Exception as e:
            print(f"Warning: Could not auto-detect smoothing: {e}")
            return 1.0
    
    def auto_detect_baseline_slope(self, distance_um: np.ndarray, intensity: np.ndarray) -> float:
        """
        Auto-detect baseline slope by analyzing intensity trend at lower percentiles.
        Returns suggested slope for threshold line.
        """
        try:
            # Ensure data is numeric float type
            intensity = np.asarray(intensity, dtype=np.float64)
            distance_um = np.asarray(distance_um, dtype=np.float64)
            
            # Use 10th percentile as baseline proxy (excludes peaks)
            window_size = max(10, len(intensity) // 20)
            baseline_points = []
            baseline_positions = []
            
            for i in range(0, len(intensity), window_size):
                window = intensity[i:i+window_size]
                if len(window) > 0:
                    # Use numpy's nanpercentile to handle any NaN values
                    percentile_val = np.nanpercentile(window, 10)
                    if not np.isnan(percentile_val):
                        baseline_points.append(percentile_val)
                        baseline_positions.append(distance_um[i + len(window)//2])
            
            if len(baseline_points) >= 2:
                # Fit line to baseline points
                slope, _ = np.polyfit(baseline_positions, baseline_points, 1)
                return float(slope)
        except Exception as e:
            print(f"Warning: Could not auto-detect slope: {e}")
        
        return 0.0
    
    def auto_detect_min_distance(self, distance_um: np.ndarray, intensity: np.ndarray, 
                                  threshold_intercept: float, threshold_slope: float,
                                  smoothing_sigma: float = 5.0) -> float:
        """
        Auto-detect optimal minimum distance between peaks.
        Analyzes peak density to suggest distance that reduces false positives.
        """
        try:
            # Ensure data is numeric float type
            intensity = np.asarray(intensity, dtype=np.float64)
            distance_um = np.asarray(distance_um, dtype=np.float64)
            
            # Apply smoothing first
            if smoothing_sigma > 0:
                intensity_smoothed = gaussian_filter1d(intensity, sigma=smoothing_sigma)
            else:
                intensity_smoothed = intensity
            
            # First pass: find peaks with very low min distance to see raw distribution
            threshold_line = threshold_intercept + threshold_slope * distance_um
            corrected = intensity_smoothed - threshold_line
            
            # Find all possible peaks with minimal constraints
            avg_spacing = np.mean(np.diff(distance_um))
            min_samples_test = max(1, int(0.1 / avg_spacing))  # 0.1 µm minimum
            
            peak_indices, _ = find_peaks(
                corrected,
                height=0,
                distance=min_samples_test,
                prominence=0.2 * np.std(intensity_smoothed)
            )
            
            if len(peak_indices) < 2:
                return 1.0  # Default if no peaks found
            
            # Calculate inter-peak distances
            peak_distances_um = np.diff(distance_um[peak_indices])
            
            if len(peak_distances_um) == 0:
                return 1.0  # Default to 1.0 µm
            
            # Use median inter-peak distance as a good estimate
            # This naturally filters out clustered false positives
            median_dist = np.median(peak_distances_um)
            
            # Use 50% of median as minimum (allows some variation but filters tight clusters)
            # But ensure at least 1.0 µm as per requirement
            suggested_min = max(1.0, min(median_dist * 0.5, 10.0))
            
            return round(suggested_min, 2)
        except Exception as e:
            print(f"Warning: Could not auto-detect min distance: {e}")
            return 1.0
    
    def get_current_params(self) -> Dict:
        """Get parameters for current axon (from saved results or defaults)."""
        axon_data = self.get_current_axon()
        key = self.get_result_key(axon_data)
        
        if key in self.results:
            return self.results[key]
        else:
            # Auto-detect only baseline slope; keep other defaults simple
            distance_um = axon_data['distance']
            intensity = axon_data['intensity']
            
            # Auto-detect baseline slope
            auto_slope = self.auto_detect_baseline_slope(distance_um, intensity)
            
            # Initial threshold (intercept)
            threshold_intercept = np.max(intensity) + 0.5 * np.std(intensity)
            return {
                'threshold_intercept': threshold_intercept,
                'threshold_slope': auto_slope,
                # Default: allow peaks to be as close as spacing resolution allows
                'min_dist_um': 0.0,
                # Default smoothing sigma
                'smoothing_sigma': 5.0,
                'min_prominence_factor': 0.3,
                'auto_detected': True  # Flag to show initial slope was auto-detected
            }
    
    def find_peaks_with_params(self, distance_um: np.ndarray, intensity: np.ndarray,
                                threshold_intercept: float, threshold_slope: float,
                                min_dist_um: float, min_prominence_factor: float,
                                smoothing_sigma: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        """Find peaks with given parameters on smoothed distribution."""
        # Apply Gaussian smoothing if sigma > 0
        if smoothing_sigma > 0:
            intensity_smoothed = gaussian_filter1d(intensity, sigma=smoothing_sigma)
        else:
            intensity_smoothed = intensity
        
        threshold_line = threshold_intercept + threshold_slope * distance_um
        corrected_intensity = intensity_smoothed - threshold_line
        min_prominence = min_prominence_factor * np.std(intensity_smoothed)
        
        # Convert min distance from microns to samples
        if len(distance_um) > 1:
            avg_spacing = np.mean(np.diff(distance_um))
            if min_dist_um <= 0:
                min_samples = 1
            else:
                min_samples = max(1, int(min_dist_um / avg_spacing))
        else:
            min_samples = 1
        
        peak_indices, properties = find_peaks(
            corrected_intensity,
            height=0,
            distance=min_samples,
            prominence=min_prominence
        )
        
        if len(peak_indices) > 0:
            peak_locs = distance_um[peak_indices]
            # Use smoothed intensity values so peaks sit on the smoothed curve
            peak_vals = intensity_smoothed[peak_indices]
        else:
            peak_locs = np.array([])
            peak_vals = np.array([])
        
        return peak_locs, peak_vals
    
    def update_plot(self):
        """Update the plot with current parameters."""
        axon_data = self.get_current_axon()
        params = self.get_current_params()
        
        distance_um = axon_data['distance']
        intensity = axon_data['intensity']
        
        # Apply smoothing for visualization
        smoothing_sigma = params.get('smoothing_sigma', 0.0)
        if smoothing_sigma > 0:
            intensity_smoothed = gaussian_filter1d(intensity, sigma=smoothing_sigma)
        else:
            intensity_smoothed = intensity
        
        # Find peaks
        peak_locs, peak_vals = self.find_peaks_with_params(
            distance_um, intensity,
            params['threshold_intercept'],
            params['threshold_slope'],
            params['min_dist_um'],
            params['min_prominence_factor'],
            smoothing_sigma
        )

        # ------------------------------------------------------------------
        # NEW: Automatically skip axons that have NO peaks with current params
        # This ensures flat/no-bouton axons are not shown in the GUI at all.
        # ------------------------------------------------------------------
        if len(peak_locs) == 0:
            # Mark as skipped and move to next axon (if possible)
            key = self.get_result_key(axon_data)
            params['skipped'] = True
            self.results[key] = params.copy()
            # If there is another axon, go to it; otherwise just return
            if self.current_idx < len(self.all_data) - 1:
                self.go_to_axon(self.current_idx + 1)
            return
        
        # Update smoothed line (always visible)
        self.smoothed_line.set_data(distance_um, intensity_smoothed)
        
        # Update threshold line
        threshold_line = params['threshold_intercept'] + params['threshold_slope'] * distance_um
        self.threshold_line.set_data(distance_um, threshold_line)
        
        # Update peaks
        if len(peak_locs) > 0:
            self.peak_scatter.set_offsets(np.c_[peak_locs, peak_vals])
            self.peak_scatter.set_visible(True)
        else:
            self.peak_scatter.set_visible(False)
        
        # Update title
        axon_length = distance_um[-1] - distance_um[0]
        bouton_density = len(peak_locs) / axon_length if axon_length > 0 else 0
        
        # Calculate average peak distance
        if len(peak_locs) > 1:
            avg_peak_dist = np.mean(np.diff(peak_locs))
            avg_dist_str = f" | Avg Peak Dist: {avg_peak_dist:.2f} µm"
        else:
            avg_dist_str = ""
        
        # Calculate density per 50 µm
        density_50um = (len(peak_locs) / axon_length) * 50.0 if axon_length > 0 else 0
        
        self.ax.set_title(
            f"File: {axon_data['file_name']} | Axon {axon_data['axon_num']}/{axon_data['total_axons']} | "
            f"Layer: {axon_data['layer']} | Mouse: {axon_data['mouse_id']}\n"
            f"Peaks: {len(peak_locs)} | Density: {bouton_density:.4f}/µm ({density_50um:.2f}/50µm){avg_dist_str} | "
            f"Progress: {self.current_idx + 1}/{len(self.all_data)}",
            fontsize=9, pad=10
        )
        
        # Update info text
        status = "✓ SAVED" if self.get_result_key(axon_data) in self.results else "Not saved"
        auto_note = " (AUTO-RUN)" if params.get('auto_detected', False) else " (MANUAL)"
        smoothing = params.get('smoothing_sigma', 0.0)
        
        # Calculate metrics for display
        if len(peak_locs) > 1:
            avg_peak_dist = np.mean(np.diff(peak_locs))
            avg_dist_text = f"Avg Peak Dist: {avg_peak_dist:.2f} µm\n"
        else:
            avg_dist_text = ""
        
        # Show only concise status/summary here to avoid duplicating parameter values
        self.info_text.set_text(
            f"Status: {status}{auto_note}\n"
            f"{avg_dist_text}"
        )
        
        self.fig.canvas.draw_idle()
    
    def update_from_sliders(self, val=None):
        """Update parameters from slider values."""
        # Skip if we're loading an axon (prevents overwriting loaded params)
        if self.loading_axon:
            return
            
        params = self.get_current_params()
        params['threshold_intercept'] = self.slider_intercept.val
        params['threshold_slope'] = self.slider_slope.val
        params['min_dist_um'] = self.slider_min_dist.val
        params['smoothing_sigma'] = self.slider_smoothing.val
        
        # Mark as manually adjusted (no longer auto-detected)
        params['auto_detected'] = False
        
        # Save to results
        axon_data = self.get_current_axon()
        key = self.get_result_key(axon_data)
        self.results[key] = params.copy()
        
        self.update_plot()
        self.save_progress()
    
    def on_mouse_press(self, event):
        """Handle mouse press for dragging threshold line."""
        if event.inaxes != self.ax:
            return
        
        # Check if click is near threshold line
        if self.threshold_line.contains(event)[0]:
            self.dragging = True
            # Determine if dragging intercept (middle) or slope (edges)
            axon_data = self.get_current_axon()
            distance_um = axon_data['distance']
            mid_point = (distance_um[0] + distance_um[-1]) / 2
            
            if abs(event.xdata - mid_point) < (distance_um[-1] - distance_um[0]) * 0.3:
                self.drag_mode = 'intercept'
            else:
                self.drag_mode = 'slope'
    
    def on_mouse_release(self, event):
        """Handle mouse release."""
        self.dragging = False
        self.drag_mode = None
    
    def on_mouse_motion(self, event):
        """Handle mouse motion for dragging."""
        if not self.dragging or event.inaxes != self.ax:
            return
        
        params = self.get_current_params()
        axon_data = self.get_current_axon()
        distance_um = axon_data['distance']
        
        if self.drag_mode == 'intercept':
            # Move threshold up/down
            params['threshold_intercept'] = event.ydata
        elif self.drag_mode == 'slope':
            # Adjust slope based on drag
            mid_point_x = (distance_um[0] + distance_um[-1]) / 2
            mid_point_y = params['threshold_intercept'] + params['threshold_slope'] * mid_point_x
            
            # Calculate new slope
            dx = event.xdata - mid_point_x
            dy = event.ydata - mid_point_y
            if abs(dx) > 1:
                params['threshold_slope'] = dy / dx
        
        # Update sliders
        self.slider_intercept.set_val(params['threshold_intercept'])
        self.slider_slope.set_val(params['threshold_slope'])
        
        # Save and update
        key = self.get_result_key(axon_data)
        self.results[key] = params.copy()
        self.update_plot()
    
    def go_to_axon(self, idx: int):
        """Navigate to specific axon index."""
        if 0 <= idx < len(self.all_data):
            self.current_idx = idx
            self.load_axon()
            self.save_progress()
    
    def on_previous(self, event):
        """Go to previous axon."""
        if self.current_idx > 0:
            self.go_to_axon(self.current_idx - 1)
    
    def on_next(self, event):
        """Go to next axon."""
        if self.current_idx < len(self.all_data) - 1:
            self.go_to_axon(self.current_idx + 1)
    
    def on_skip(self, event):
        """Skip current axon (mark with defaults and move to next)."""
        params = self.get_current_params()
        axon_data = self.get_current_axon()
        key = self.get_result_key(axon_data)
        
        # Save with current params but mark as skipped
        params['skipped'] = True
        self.results[key] = params.copy()
        self.save_progress()
        
        # Move to next
        self.on_next(event)
    
    def on_save_quit(self, event):
        """Save and quit."""
        self.save_progress()
        print("\n✓ Progress saved. Closing...")
        plt.close(self.fig)
    
    def on_auto_detect(self, event):
        """Re-run auto-detection for current axon with individual optimization."""
        axon_data = self.get_current_axon()
        distance_um = axon_data['distance']
        intensity = axon_data['intensity']
        
        # Get current params
        params = self.get_current_params()
        
        # Auto-detect all parameters individually for this axon
        auto_slope = self.auto_detect_baseline_slope(distance_um, intensity)
        threshold_intercept = params['threshold_intercept']  # Keep current intercept
        
        # Auto-detect optimal smoothing
        auto_smoothing = self.auto_detect_optimal_smoothing(distance_um, intensity)
        
        # Auto-detect minimum distance with optimized smoothing
        auto_min_dist = self.auto_detect_min_distance(
            distance_um, intensity, threshold_intercept, auto_slope, auto_smoothing
        )
        
        # Update parameters
        params['threshold_slope'] = auto_slope
        params['min_dist_um'] = max(1.0, auto_min_dist)
        params['smoothing_sigma'] = auto_smoothing
        params['auto_detected'] = True
        
        # Update sliders
        self.slider_slope.set_val(auto_slope)
        self.slider_min_dist.set_val(params['min_dist_um'])
        self.slider_smoothing.set_val(auto_smoothing)
        
        # Save and update
        key = self.get_result_key(axon_data)
        self.results[key] = params.copy()
        self.update_plot()
        self.save_progress()
    
    def load_axon(self):
        """Load current axon data into GUI."""
        # Set flag to prevent slider callbacks from overwriting loaded params
        self.loading_axon = True
        
        axon_data = self.get_current_axon()
        params = self.get_current_params()
        
        # Debug: print parameters for this axon
        print(f"\nLoading Axon {self.current_idx + 1}: "
              f"slope={params['threshold_slope']:.3f}, "
              f"σ={params.get('smoothing_sigma', 5.0):.1f}, "
              f"minDist={params['min_dist_um']:.2f}")
        
        # Update data line
        distance_um = axon_data['distance']
        intensity = axon_data['intensity']
        self.data_line.set_data(distance_um, intensity)
        
        # Reset axis limits - Y axis always from 0 to 1.5x max intensity
        self.ax.set_xlim(distance_um[0], distance_um[-1])
        y_max = intensity.max() * 1.5  # Consistent headroom above peaks
        self.ax.set_ylim(0, y_max)
        
        # Update sliders with fixed range (0 to 2x max intensity)
        intercept_max = intensity.max() * 2.0
        self.slider_intercept.valmin = 0.0
        self.slider_intercept.valmax = intercept_max
        self.slider_intercept.set_val(params['threshold_intercept'])
        self.slider_slope.set_val(params['threshold_slope'])
        self.slider_min_dist.set_val(params['min_dist_um'])
        self.slider_smoothing.set_val(params.get('smoothing_sigma', 1.0))

        self.update_plot()
        
        # Re-enable slider callbacks
        self.loading_axon = False
    
    def auto_run_all(self):
        """Auto-run all axons with individually optimized parameters."""
        print("\n" + "="*70)
        print("AUTO-RUN MODE: Processing all axons with individual optimization...")
        if FORCE_REPROCESS:
            print("FORCE_REPROCESS enabled: Re-optimizing ALL axons from scratch")
        print("="*70)
        
        processed_count = 0
        skipped_count = 0
        
        for idx in range(len(self.all_data)):
            self.current_idx = idx
            axon_data = self.get_current_axon()
            key = self.get_result_key(axon_data)
            
            # Skip if already processed (unless FORCE_REPROCESS is enabled)
            if key in self.results and not FORCE_REPROCESS:
                skipped_count += 1
                continue
            
            distance_um = axon_data['distance']
            intensity = axon_data['intensity']
            
            # Auto-detect parameters individually for this axon
            auto_slope = self.auto_detect_baseline_slope(distance_um, intensity)
            threshold_intercept = np.mean(intensity) + 0.5 * np.std(intensity)
            
            # Try different smoothing values and pick the best
            best_smoothing = self.auto_detect_optimal_smoothing(distance_um, intensity)
            
            # Auto-detect minimum distance with optimal smoothing
            auto_min_dist = self.auto_detect_min_distance(
                distance_um, intensity, threshold_intercept, auto_slope, best_smoothing
            )
            
            # Save individually optimized parameters
            params = {
                'threshold_intercept': threshold_intercept,
                'threshold_slope': auto_slope,
                # Default minimum now 1.0 µm
                'min_dist_um': max(1.0, auto_min_dist),
                'smoothing_sigma': best_smoothing,
                'min_prominence_factor': 0.3,
                'auto_detected': True
            }
            
            self.results[key] = params.copy()
            processed_count += 1
            
            # Show some details every 10 axons
            if processed_count % 10 == 0:
                print(f"  Processed {processed_count} new axons (example: slope={auto_slope:.3f}, σ={best_smoothing:.1f})...")
        
        print(f"✓ Auto-processed {processed_count} new axons")
        if skipped_count > 0:
            print(f"  Skipped {skipped_count} already-processed axons from previous run")
        print(f"  Total axons with parameters: {len(self.results)}")
        print("  You can now review and edit them in the GUI")
        self.save_progress()
    
    def run(self):
        """Run the interactive GUI."""
        # Create figure
        self.fig = plt.figure(figsize=(14, 10))
        plt.subplots_adjust(bottom=0.25, right=0.85, top=0.98)
        
        # Main plot
        self.ax = plt.axes([0.1, 0.35, 0.7, 0.6])
        
        # Initial data
        axon_data = self.get_current_axon()
        distance_um = axon_data['distance']
        intensity = axon_data['intensity']
        params = self.get_current_params()
        
        # Plot elements - only show smoothed data
        # Original raw data (hidden, kept for reference)
        self.data_line, = self.ax.plot(distance_um, intensity, 'k-', linewidth=1.0, 
                                        label='Original', alpha=0.0, visible=False)
        
        # Smoothed line - this is what we show and analyze
        smoothing_sigma = params.get('smoothing_sigma', 0.0)
        if smoothing_sigma > 0:
            intensity_smoothed = gaussian_filter1d(intensity, sigma=smoothing_sigma)
        else:
            intensity_smoothed = intensity
        self.smoothed_line, = self.ax.plot(distance_um, intensity_smoothed, 'k-', 
                                           linewidth=1.5, label='Intensity (smoothed)', alpha=1.0)
        
        threshold_line_data = params['threshold_intercept'] + params['threshold_slope'] * distance_um
        self.threshold_line, = self.ax.plot(distance_um, threshold_line_data, 'b--', 
                                            linewidth=2.5, label='Threshold (draggable)', 
                                            picker=5, alpha=0.8)
        self.peak_scatter = self.ax.scatter([], [], c='red', s=80, zorder=5, 
                                           marker='o', edgecolors='darkred', 
                                           linewidths=1.5, label='Peaks')
        
        self.ax.set_xlabel('Distance (µm)', fontsize=11)
        self.ax.set_ylabel('Fluorescence Intensity', fontsize=11)
        self.ax.legend(loc='upper right', fontsize=9)
        self.ax.grid(True, alpha=0.3)
        
        # Info text box
        self.info_text = self.fig.text(0.87, 0.7, '', fontsize=9, 
                                       verticalalignment='top',
                                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Sliders
        ax_intercept = plt.axes([0.15, 0.27, 0.50, 0.02])
        ax_slope = plt.axes([0.15, 0.23, 0.50, 0.02])
        ax_min_dist = plt.axes([0.15, 0.19, 0.50, 0.02])
        ax_smoothing = plt.axes([0.15, 0.15, 0.50, 0.02])
        
        # Fixed range from 0 to 2x the max intensity (consistent scale across all axons)
        intercept_max = intensity.max() * 2.0
        self.slider_intercept = Slider(
            ax_intercept, 'Threshold\nIntercept',
            0.0, intercept_max * 3.0,
            valinit=params['threshold_intercept'],
            valstep=max(intercept_max * 3.0 / 500, 1e-6)
        )
        self.slider_slope = Slider(
            ax_slope, 'Threshold\nSlope',
            -2000.0, 2000.0,
            valinit=params['threshold_slope'],
            valstep=1.0
        )
        self.slider_min_dist = Slider(
            ax_min_dist, 'Min Peak\nDist (µm)',
            0.0, 20.0,
            valinit=params['min_dist_um'],
            valstep=0.1
        )
        self.slider_smoothing = Slider(
            ax_smoothing, 'Smoothing\n(σ)',
            0.0, 10.0,
            valinit=params.get('smoothing_sigma', 1.0),
            valstep=0.1
        )

        # Navigation buttons (bottom row)
        ax_prev = plt.axes([0.08, 0.08, 0.08, 0.035])
        ax_skip = plt.axes([0.19, 0.08, 0.08, 0.035])
        ax_auto = plt.axes([0.30, 0.08, 0.10, 0.035])
        ax_next = plt.axes([0.45, 0.08, 0.08, 0.035])
        ax_save = plt.axes([0.57, 0.08, 0.12, 0.035])

        self.btn_previous = Button(ax_prev, '← Previous', color='lightblue')
        self.btn_skip = Button(ax_skip, 'Skip', color='lightyellow')
        self.btn_auto_detect = Button(ax_auto, '🔍 Auto-Detect', color='lightgreen')
        self.btn_next = Button(ax_next, 'Next →', color='lightgreen')
        self.btn_save_quit = Button(ax_save, 'Save & Quit', color='lightcoral')

        self.btn_previous.on_clicked(self.on_previous)
        self.btn_skip.on_clicked(self.on_skip)
        self.btn_auto_detect.on_clicked(self.on_auto_detect)
        self.btn_next.on_clicked(self.on_next)
        self.btn_save_quit.on_clicked(self.on_save_quit)

        # Connect sliders
        self.slider_intercept.on_changed(self.update_from_sliders)
        self.slider_slope.on_changed(self.update_from_sliders)
        self.slider_min_dist.on_changed(self.update_from_sliders)
        self.slider_smoothing.on_changed(self.update_from_sliders)
        
        # Mouse events for dragging
        self.fig.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_motion)
        
        # Instructions - concise, slider/drag only
        instructions = (
            "CONTROLS: DRAG blue threshold line | SLIDERS adjust params | "
            "Black line = smoothed intensity | AUTO-SAVES"
        )
        self.fig.text(0.5, 0.03, instructions, fontsize=8, 
                     ha='center', verticalalignment='bottom',
                     bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
        
        # Initial load
        self.load_axon()
        
        plt.show()
        
        # After closing, save final progress
        self.save_progress()
        
        return self.results


def parse_filename_metadata(filename: str) -> Tuple[str, str, str, str]:
    """
    Extract metadata from filename including cortical layer (L1/L5).
    
    Returns: (mouse_id, brain_region, image_num, layer)
    """
    basename = os.path.splitext(os.path.basename(filename))[0]
    
    # Extract layer (L1 or L5)
    layer = 'Unknown'
    if '_L5_' in basename or '_L5' in basename:
        layer = 'L5'
    elif '_L1_' in basename or '_L1' in basename:
        layer = 'L1'
    
    # Parse other components
    parts = basename.replace('Values_', '').replace('Values__', '').replace('_', ' ').split()
    
    mouse_id = 'Unknown'
    brain_region = 'Unknown'
    image_num = 'Unknown'
    
    if len(parts) >= 2:
        mouse_id = parts[0]  # e.g., "Adult"
        if len(parts) >= 3:
            brain_region = f"{parts[1]}_{parts[2]}"  # e.g., "1_sl1"
        if len(parts) >= 4:
            image_num = parts[-1]  # Last part
    
    return mouse_id, brain_region, image_num, layer


def load_all_axons(input_dir: str) -> List[Dict]:
    """Load all axons from all files in directory."""
    # Find all data files
    file_patterns = [
        os.path.join(input_dir, '*.xlsx'),
        os.path.join(input_dir, '*.xls'),
        os.path.join(input_dir, '*.csv')
    ]
    
    data_files = []
    for pattern in file_patterns:
        data_files.extend(glob.glob(pattern))
    
    # Filter to only include files with "values" in the filename (case-insensitive)
    data_files = [f for f in data_files if 'values' in os.path.basename(f).lower()]
    data_files = sorted(set(data_files))
    
    if not data_files:
        print(f"No Excel or CSV files with 'values' in filename found in {input_dir}")
        return []
    
    print(f"Found {len(data_files)} data file(s) with 'values' in filename")
    
    all_axons = []
    
    for file_path in data_files:
        # Skip temporary Excel files (start with ~$)
        if os.path.basename(file_path).startswith('~$'):
            continue
            
        try:
            # Load file
            if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
                data = pd.read_excel(file_path, header=0)
            elif file_path.endswith('.csv'):
                data = pd.read_csv(file_path, header=0)
            
            # Convert to numeric
            for col in data.columns:
                data[col] = pd.to_numeric(data[col], errors='coerce')
            
            # Parse metadata
            mouse_id, brain_region, image_num, layer = parse_filename_metadata(file_path)
            
            # Build column pairs
            column_names = data.columns.tolist()
            axon_col_pairs = []
            
            if len(column_names) >= 2:
                axon_col_pairs.append((0, 1))
            
            for i in range(2, len(column_names), 2):
                if i + 1 < len(column_names):
                    axon_col_pairs.append((i, i + 1))
            
            # Extract each axon
            for ax_idx, (dist_col, intensity_col) in enumerate(axon_col_pairs):
                ax_num = ax_idx + 1
                
                dist_vals = data.iloc[:, dist_col].values
                intensity_vals = data.iloc[:, intensity_col].values
                
                # Clean NaNs
                valid_idx = ~np.isnan(dist_vals) & ~np.isnan(intensity_vals)
                dist_vals = dist_vals[valid_idx]
                intensity_vals = intensity_vals[valid_idx]
                
                if len(dist_vals) == 0:
                    continue
                
                # ------------------------------------------------------------------
                # NEW: Automatically exclude axons with flat signal / no peaks
                # We do a quick peak check here and completely skip axons that have
                # zero detectable peaks, so they are not considered anywhere.
                # ------------------------------------------------------------------
                try:
                    # Basic smoothing with small sigma
                    intensity_arr = np.asarray(intensity_vals, dtype=np.float64)
                    dist_arr = np.asarray(dist_vals, dtype=np.float64)
                    
                    if len(intensity_arr) > 1 and not np.all(intensity_arr == intensity_arr[0]):
                        # Light smoothing
                        intensity_smoothed = gaussian_filter1d(intensity_arr, sigma=1.0)
                        
                        # Simple threshold at mean + 0.5*std
                        thr = np.mean(intensity_smoothed) + 0.5 * np.std(intensity_smoothed)
                        corrected = intensity_smoothed - thr
                        
                        # Convert 1 µm to samples
                        avg_spacing = np.mean(np.diff(dist_arr)) if len(dist_arr) > 1 else 1.0
                        min_samples = max(1, int(1.0 / avg_spacing))
                        
                        # Find peaks
                        peak_indices, _ = find_peaks(
                            corrected,
                            height=0,
                            distance=min_samples,
                            prominence=0.3 * np.std(intensity_smoothed)
                        )
                        if len(peak_indices) == 0:
                            # No peaks at all → skip this axon entirely
                            continue
                    else:
                        # Completely flat or single-point signal → skip
                        continue
                except Exception as e:
                    # If anything goes wrong in this quick check, fall back to keeping the axon
                    pass
                
                all_axons.append({
                    'file_path': file_path,
                    'file_name': os.path.basename(file_path),
                    'mouse_id': mouse_id,
                    'brain_region': brain_region,
                    'image_num': image_num,
                    'layer': layer,
                    'axon_num': ax_num,
                    'total_axons': len(axon_col_pairs),
                    'distance_col': column_names[dist_col],
                    'intensity_col': column_names[intensity_col],
                    'distance': dist_vals,
                    'intensity': intensity_vals,
                    'axon_length_um': dist_vals[-1] - dist_vals[0]
                })
                
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
            continue
    
    print(f"Loaded {len(all_axons)} total axons from all files")
    return all_axons


def export_results(results: Dict, all_axons: List[Dict], output_file: str):
    """Export all results to Excel file."""
    rows = []
    
    for axon_data in all_axons:
        key = f"{axon_data['file_name']}__axon_{axon_data['axon_num']}"
        
        if key in results:
            params = results[key]
            
            # Find peaks with final parameters
            intensity = axon_data['intensity']
            distance_um = axon_data['distance']
            
            # Apply smoothing if specified
            smoothing_sigma = params.get('smoothing_sigma', 0.0)
            if smoothing_sigma > 0:
                intensity_smoothed = gaussian_filter1d(intensity, sigma=smoothing_sigma)
            else:
                intensity_smoothed = intensity
            
            threshold_line = params['threshold_intercept'] + params['threshold_slope'] * distance_um
            corrected = intensity_smoothed - threshold_line
            
            avg_spacing = np.mean(np.diff(distance_um))
            if params['min_dist_um'] <= 0:
                min_samples = 1
            else:
                min_samples = max(1, int(params['min_dist_um'] / avg_spacing))
            min_prominence = params['min_prominence_factor'] * np.std(intensity_smoothed)
            
            peak_indices, _ = find_peaks(
                corrected, height=0, distance=min_samples, prominence=min_prominence
            )
            
            if len(peak_indices) > 0:
                peak_locs = distance_um[peak_indices]
                # Use smoothed intensity for peak values (so they match the displayed curve)
                peak_vals = intensity_smoothed[peak_indices]
            else:
                peak_locs = np.array([])
                peak_vals = np.array([])
            
            peak_count = len(peak_locs)
            
            # If this axon has 0 detected boutons/peaks, exclude it from export and downstream analysis
            if peak_count == 0:
                continue
            
            # Calculate three versions of spike density
            # 1. Total spike count for the length of the axon
            spike_count_total = peak_count
            
            # 2. Spike density per micron
            if axon_data['axon_length_um'] > 0:
                spike_density_per_um = peak_count / axon_data['axon_length_um']
            else:
                spike_density_per_um = np.nan
            
            # 3. Spike density per 50 microns
            if axon_data['axon_length_um'] > 0:
                spike_density_per_50um = (peak_count / axon_data['axon_length_um']) * 50.0
            else:
                spike_density_per_50um = np.nan
            
            # Calculate average inter-peak distance
            if peak_count > 1:
                inter_peak_distances = np.diff(peak_locs)
                avg_peak_distance = np.mean(inter_peak_distances)
            else:
                avg_peak_distance = np.nan
            
            # Keep old names for backward compatibility
            bouton_density = spike_density_per_um
            density_per_50um = spike_density_per_50um
            
            peak_locs_str = ';'.join([f'{loc:.4f}' for loc in peak_locs])
            peak_vals_str = ';'.join([f'{val:.2f}' for val in peak_vals])
            
            rows.append({
                'MouseID': axon_data['mouse_id'],
                'BrainRegion': axon_data['brain_region'],
                'ImageNum': axon_data['image_num'],
                'Layer': axon_data['layer'],
                'FileName': axon_data['file_name'],
                'Axon': axon_data['axon_num'],
                'DistanceColumn': axon_data['distance_col'],
                'IntensityColumn': axon_data['intensity_col'],
                'ThresholdIntercept': params['threshold_intercept'],
                'ThresholdSlope': params['threshold_slope'],
                'MinDist_um': params['min_dist_um'],
                'SmoothingSigma': params.get('smoothing_sigma', 0.0),
                'AxonLength_um': axon_data['axon_length_um'],
                # Three versions of spike/bouton density
                'SpikeCount_Total': spike_count_total,
                'SpikeDensity_per_um': spike_density_per_um,
                'SpikeDensity_per_50um': spike_density_per_50um,
                'AvgPeakDistance_um': avg_peak_distance,
                'PeakLocations_um': peak_locs_str,
                'PeakIntensities': peak_vals_str,
                'Skipped': params.get('skipped', False),
                # Keep old names for backward compatibility
                'PeakCount': peak_count,
                'BoutonDensity_per_um': bouton_density,
                'PeakDensity_per_50um': density_per_50um,
            })
    
    if rows:
        df = pd.DataFrame(rows)
        df.to_excel(output_file, index=False)
        print(f"\n✓ Results exported to: {output_file}")
        print(f"  Total axons analyzed: {len(rows)}")
        print(f"  Total peaks detected: {df['PeakCount'].sum()}")
        return df
    else:
        print("\nNo results to export.")
        return None


def main():
    """Main function."""
    # Use directories from configuration at top of file
    input_dir = INPUT_DIRECTORY
    progress_file = os.path.join(OUTPUT_DIRECTORY, 'analysis_progress.json')
    output_file = os.path.join(OUTPUT_DIRECTORY, 'consolidated_axon_peaks_interactive.xlsx')
    
    print("="*70)
    print("Enhanced Interactive Axon Peak Analysis")
    print("="*70)
    
    # Load all axons
    all_axons = load_all_axons(input_dir)
    
    if not all_axons:
        print("No axons to analyze.")
        return
    
    # Load previous results to filter out already-analyzed axons
    previous_results = {}
    if os.path.exists(progress_file):
        try:
            with open(progress_file, 'r') as f:
                data = json.load(f)
                previous_results = data.get('results', {})
        except Exception as e:
            print(f"Warning: Could not load progress file: {e}")
    
    # Filter out axons that are already analyzed (unless FORCE_REPROCESS)
    if FORCE_REPROCESS:
        axons_to_analyze = all_axons
        print(f"\nFORCE_REPROCESS enabled: Will re-analyze all {len(all_axons)} axons")
    else:
        axons_to_analyze = []
        for axon in all_axons:
            key = f"{axon['file_name']}__axon_{axon['axon_num']}"
            if key not in previous_results:
                axons_to_analyze.append(axon)
        print(f"\nProgress status: {len(previous_results)}/{len(all_axons)} axons already analyzed")
        print(f"  Showing {len(axons_to_analyze)} remaining axons in GUI (starting from axon 1)")
    
    if not axons_to_analyze:
        print("\nAll axons have been analyzed! Set FORCE_REPROCESS=True to re-analyze.")
        # Still export the results even if nothing new to analyze
        if previous_results:
            print("\nExporting existing results...")
            export_results(previous_results, all_axons, output_file)
        return
    
    # Run analysis with filtered axon list
    # The analyzer will load previous_results automatically in load_progress()
    analyzer = InteractiveAxonAnalyzer(axons_to_analyze, progress_file)
    # Ensure previous results are preserved (they should already be loaded, but make sure)
    analyzer.results.update(previous_results)
    
    # Explicitly start from axon 1 (index 0) when opening GUI
    analyzer.current_idx = 0
    if len(axons_to_analyze) > 0:
        first_axon = axons_to_analyze[0]
        print(f"\n✓ GUI will start at axon 1: {first_axon['file_name']} - Axon {first_axon['axon_num']}")
    
    # Auto-run all axons if enabled
    if AUTO_RUN_ALL:
        if FORCE_REPROCESS:
            print("Will re-optimize all axons (FORCE_REPROCESS=True)")
        analyzer.auto_run_all()
        # Reset to axon 1 after auto-run
        analyzer.current_idx = 0
        if len(axons_to_analyze) > 0:
            print(f"✓ Resetting to axon 1 after auto-run")
    
    # Open GUI for review/editing (will start at axon 1)
    results = analyzer.run()
    
    # Export results (analyzer.results already includes all previous + new results)
    print("\nExporting results...")
    df = export_results(analyzer.results, all_axons, output_file)
    
    if df is not None:
        # Summary by layer
        print("\n" + "="*70)
        print("SUMMARY BY LAYER")
        print("="*70)
        layer_summary = df.groupby('Layer').agg({
            'PeakCount': ['sum', 'mean', 'std'],
            'BoutonDensity_per_um': ['mean', 'std'],
            'Axon': 'count'
        })
        print(layer_summary)
        
        # Save summary
        summary_file = output_file.replace('.xlsx', '_summary_by_layer.xlsx')
        layer_summary.to_excel(summary_file)
        print(f"\n✓ Layer summary saved to: {summary_file}")
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)


if __name__ == "__main__":
    main()

