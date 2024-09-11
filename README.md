# Power Transformations to Improve Normality or Symmetry

This project is a replication of the method presented in the influential paper, **"A New Family of Power Transformations to Improve Normality or Symmetry"** by In-Kwon Yeo and Richard A. Johnson (2000). The main objective of this work is to introduce a robust transformation technique that can be applied across the entire real line, addressing limitations found in earlier approaches, such as the Box-Cox transformation, which are constrained to positive values.

## Overview

In this repository, we implement the R version of the new family of transformations proposed by Yeo and Johnson. These transformations are specifically designed to reduce skewness and help approximate normality in data, making them invaluable in statistical analysis and modeling where normality assumptions are critical. The method also provides a solution for handling data that includes both positive and negative values, extending the versatility of traditional transformation techniques.

## Features
- Implementation of power transformations that handle both positive and negative values.
- Detailed comparison with the Box-Cox transformation.
- Real-world data applications to demonstrate the transformation's effectiveness in improving normality and reducing skewness.
- Visualizations to help interpret the results of the transformations.

## Motivation
Statistical models often require data to be normally distributed. Traditional transformations like Box-Cox fail to handle datasets with negative values or provide sufficient skewness correction. This project highlights how the Yeo-Johnson transformation not only addresses these issues but also offers improved flexibility and applicability in practical scenarios.

## Paper Reference
- In-Kwon Yeo and Richard A. Johnson (2000). *A New Family of Power Transformations to Improve Normality or Symmetry*, Biometrika, 87(4), 954-959. [Link to the paper](https://www.jstor.org/stable/2673623)
