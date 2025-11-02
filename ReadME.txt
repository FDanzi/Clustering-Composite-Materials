=======================================================================
 Cracking Composite Design: Trace-Driven Scaling for Cluster-Based 
 Universal Design
=======================================================================

 Author:
     Francesco Danzi

 Reference:
     F. Danzi, "Cracking Composite Design: Trace-Driven Scaling for 
     Cluster-Based Universal Design," Composite Structures, 2025.
     DOI: [insert DOI once published]

-----------------------------------------------------------------------
 Overview
-----------------------------------------------------------------------
This repository contains MATLAB scripts accompanying the paper
"Cracking Composite Design: Trace-Driven Scaling for Cluster-Based
Universal Design."

The codes implement the π-space clustering, sensitivity analysis,
and validation procedures described in the manuscript. Together,
they reproduce Tables 1–3, Figure 1, and results reported in Appendix 1.

-----------------------------------------------------------------------
 Workflow
-----------------------------------------------------------------------
To reproduce the results or apply the framework to a new dataset,
run the scripts in the following order:

 1. Composite_Materials_Clusters.m
    - Clusters the material systems using the K-means algorithm.
    - Computes the cluster centroids and plots the clusters with their
      respective centroids.
    - Includes all materials listed in Table 1 and reproduces the results
      shown in Table 3, Figure 1, and Appendix 1.
    - The number of clusters can be modified in the input section
      (see the line indicated in the script).

 2. davies_bouldin.m
    - Computes the Davies–Bouldin Index (DBI) to evaluate the 
      cohesiveness and separation of clusters.

 3. my_silhouette.m
    - Computes the Silhouette Index to assess intra-cluster compactness
      and inter-cluster separation.

 4. Sensitivity Analysis (run independently)
    - Sensitivity_Trace_EigenvalueNorm.m
         Computes sensitivity with respect to the Eigenvalue Norm
         (Section 3.2 of the manuscript).
    - Sensitivity_Trace_FrobeniusNorm.m
         Computes sensitivity with respect to the Frobenius Norm.
    - Sensitivity_Trace_SpectralNorm.m
         Computes sensitivity with respect to the Spectral Norm.

-----------------------------------------------------------------------
 Output
-----------------------------------------------------------------------
The scripts generate:
    - Cluster plots and centroids (corresponding to Figure 1)
    - Davies–Bouldin and Silhouette indices (Table 2 metrics)
    - Sensitivity results (Section 3.2 of the manuscript)

-----------------------------------------------------------------------
 License
-----------------------------------------------------------------------
All scripts are released under the MIT License for research and 
educational use.

-----------------------------------------------------------------------
 Citation
-----------------------------------------------------------------------
If you use this code, please cite:

    F. Danzi, "Cracking Composite Design: Trace-Driven Scaling for 
    Cluster-Based Universal Design," Composite Structures, 2025.
    DOI: [insert DOI once published]

-----------------------------------------------------------------------
 Repository Link
-----------------------------------------------------------------------
GitHub Repository:
    https://github.com/FDanzi/Clustering-Composite-Materials
-----------------------------------------------------------------------
