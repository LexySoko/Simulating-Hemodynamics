# Simulating-Hemodynamics

This repository contains MATLAB code to process wide-field optical mapping (WFOM) data for calculating changes in total hemoglobin concentration. The HRF helps in understanding the neurovascular coupling by analyzing the changes in fluorescence and hemoglobin concentrations.

## Description

The script `process_hemodynamics.m` performs the following steps:

1. Fetches raw fluorescence data for violet light.
2. Computes the mode-based baseline for violet fluorescence.
3. Calculates the fractional fluorescence for violet light.
4. Computes changes in total hemoglobin concentration.
5. Uses least-square deconvolution to estimate the HRF.

## Usage

1. Ensure you have the necessary data and dependencies.
2. Run the `process_hemodynamics.m` script in MATLAB.

## Example

The script includes an example to plot the raw fluorescence and the computed change in total hemoglobin for a selected pixel, and to estimate the HRF using least-square deconvolution.

## Code

Please refer to `process_hemodynamics.m` for the implementation.



