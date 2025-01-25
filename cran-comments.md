Dear CRAN volunteers,

Thank you for reviewing this submission. It contains an update to the prior3D R
package. Specifically, the update contains assorted minor improvements, bug
fixes, and updates to the package documentation. It also addresses the NOTEs
currently produced during CRAN package checks related to dependence of
R >= 4.1.0.

Best,

Christos Adam

## CRAN check notes

* checking DESCRIPTION meta-information ... NOTE
  Missing dependency on R >= 4.1.0 because package code uses the pipe
  |> or function shorthand \(...) syntax added in R 4.1.0.
  File(s) using such syntax:
    ‘rfunctions.R’

**Added the following field in the DESCRIPTION:**

Depends:
    R (>= 4.1.0)

## R CMD check results

0 errors | 0 warnings | 0 notes
