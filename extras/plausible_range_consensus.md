# Learning Plausible Measurement Ranges via Consensus Outlier Detection

## Context

This document summarizes a design discussion on how to learn **biologically plausible value ranges** for a single OMOP CDM measurement concept using **multiple univariate outlier detection methods**.

The goal is to move beyond:
- Clinical "normal ranges"
- Simple statistical thresholds

Toward:
> A **robust, data-driven estimate of plausible min/max values**

---

## Problem Statement

Given a dataset for a single measurement concept with:
- Measurement values
- Frequency of values

We want to:
1. Detect outliers using multiple statistical methods
2. Combine results using a **consensus (voting) approach**
3. Derive a **plausible value range** that is:
   - Robust to noise
   - Resistant to ETL/unit errors
   - Inclusive of rare but biologically plausible values

---

## Candidate Methods

The following univariate outlier detection methods are used:

| Method | Description |
|------|--------|
| Tukey fences | IQR-based, distribution-free |
| Quantile thresholds | Percentile cutoffs |
| Z-score | Mean/std-based |
| Modified Z-score | MAD-based (robust) |
| Generalized ESD | Statistical test for outliers |

---

## Key Design Principle

> Do NOT average bounds across methods.

Instead:
> Use **consensus voting on outlier flags**, then derive bounds from retained values.

---

## Proposed Approach

### Step 1: Compute bounds per method

Each method produces:

[min_m, max_m]


---

### Step 2: Flag observations

For each value `x_i`:

flag_m(i) = 1 if x_i < min_m OR x_i > max_m
flag_m(i) = 0 otherwise


---

### Step 3: Voting / consensus

Aggregate across methods:
score(i) = sum(flag_m(i))


Define thresholds:
- Strict: score ≥ 4
- Moderate: score ≥ 3  ← recommended
- Lenient: score ≥ 2

---

### Step 4: Define valid data
valid_points = { x_i where score(i) < threshold }
---

### Step 5: Learn plausible range

Preferred (robust):
min_plausible = quantile(valid_points, 0.001)
max_plausible = quantile(valid_points, 0.999)


Alternative (less robust):


min_plausible = min(valid_points)
max_plausible = max(valid_points)


---

## Incorporating Frequency

Frequency provides important signal about reliability.

### Option A: Weighted voting
weighted_score(i) = sum(flag_m(i)) * freq_i


---

### Option B: Density filtering

Remove values that are:
- Low frequency
- AND flagged by multiple methods

---

### Option C: Support-based filtering
support(x) = frequency(x) / total_observations


Keep values that:
- Have high support
- OR are not flagged

---

## Important Distinction

> Plausible ≠ Frequent

Examples:
- Potassium = 6.5 → rare but plausible
- Potassium = 50 → implausible

The method must:
- Remove structurally inconsistent values
- Retain rare but valid observations

---

## Enhancements

### 1. Log transformation

For skewed distributions:
x_log = log(x)


Apply detection → back-transform

---

### 2. Asymmetric bounds

Lower and upper tails behave differently:

- Learn bounds separately
- Use different thresholds per tail if needed

---

### 3. Bootstrap stability

- Recompute plausible range on resampled data
- Ensure bounds are stable

---

### 4. Stratification (optional)

Split by:
- Inpatient vs outpatient
- Age groups

Then recombine

---

## Failure Modes

Be cautious when:

- Data contains unit inconsistencies
- Distribution is multi-modal
- Dataset is heavily corrupted

In these cases:
- Preprocessing is required before applying this method

---

## Example (R-like pseudocode)

```r
# Step 1: compute bounds
bounds <- list(
  tukey = c(min_t, max_t),
  quantile = c(min_q, max_q),
  z = c(min_z, max_z),
  modz = c(min_mz, max_mz),
  esd = c(min_esd, max_esd)
)

# Step 2: flag values
flags <- sapply(bounds, function(b) {
  (x < b[1]) | (x > b[2])
})

# Step 3: consensus score
score <- rowSums(flags)

# Step 4: filter valid values
valid <- x[score < 3]

# Step 5: plausible range
min_plausible <- quantile(valid, 0.001)
max_plausible <- quantile(valid, 0.999)

Relation to Ensemble Methods

This approach is analogous to:

Ensemble anomaly detection
Majority voting classifiers

Key idea:

Combine diverse statistical assumptions to reduce bias and false positives

Future Extensions
Add multivariate methods (e.g., LOF, Isolation Forest)
Build concept-level anomaly graphs
Integrate with OHDSI DataQualityDashboard
Learn ranges across sites for network-wide validation
Summary
Use multiple univariate detectors
Convert outputs to binary flags
Apply consensus voting
Derive plausible range from retained values
Incorporate frequency to improve robustness

This provides a simple, explainable, and extensible framework for learning biologically plausible measurement ranges in OMOP CDM.