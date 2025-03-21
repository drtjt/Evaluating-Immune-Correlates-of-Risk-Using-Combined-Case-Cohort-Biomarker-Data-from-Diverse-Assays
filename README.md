# Simulation Study for "Meta-Analysis Approaches for Evaluating Immune Correlates of Risk Using Combined Case-Cohort Biomarker Data from Diverse Assays"

This repository contains the code used for the simulation study in our paper:

> **Bridging Biomarker Measurements to Identify Cost-Effective Biomarkers**  
> Authors: Trevor J Thomson, Ying Huang, and Yunda Huang  
> Submitted to "Statistics in Medicine"

## 📌 Overview
This repository provides:
- Code to generate simulated data under a probit regression model, and obtain information to reproduce Tables 1 and 2
- Code to generate simulated data under a logistic regression model, and obtain information to reproduce Tables 1 and 3

## 🛠️ Requirements
The code was written in **R** and uses the following packages:
- `MASS`
- `nleqslv`
- `mvtnorm`
- `Rcpp`
- `RcppArmadillo`
- `stringr`
- `numDeriv`
- `cubature`
- `survey` 

The **R** script calls **C++** files to execute our proposed estimation strategy.

## 📧 Contact
For questions, feedback, or collaboration inquiries, feel free to reach out to:

👤 Trevor Thomson

📍 Fred Hutchinson Cancer Center

📬 Email: tthomson@fredhutch.org
