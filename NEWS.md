# nQuack 1.0.5

* Expected-maximization code updates to export the probability of loci belonging to each mixture, denoted `pir`.
    * When, `glreturn = TRUE` each emstep function will now return a matrix where each row is a locus and each column is a mixture. The observations represent the probability that the locus observed belongs to each mixture. In the supplement for our publication, this is the Zi value in the E step.

# nQuack 1.0.4

* Documentation updates including:
  * Updated README.md with new features and usage instructions.
  * Improved examples and tutorials to help users use nQuack correctly!
  * Input data expanded to include VCFs and output from Qploidy2. 
  * Updated formatting on function documentation. 
  
* Preparing package for CRAN:
  * Adding checks for dependencies to ensure smooth installation and usage.
  * Moved full example data to Zenodo to decrease file size of the package and improve installation speed.
  
