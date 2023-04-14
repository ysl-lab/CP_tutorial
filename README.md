# CP_tutorial
This public GitHub repository provides a tutorial for setting-up BE-META simulations of cyclic peptides.<sup>1,2</sup> 

•	To access the step-by-step tutorials, please see the Word documents in the /Tutorials directory, which includes documents with step-by-step explanations for building cyclic peptides in Chimera, solvating and equilibrating the cyclic peptide in GROMACS<sup>3</sup> using the RSFF2 forcefield<sup>4</sup> , setting-up the BE-META simulation<sup>5,6</sup>, and analyzing the trajectories using dihedral PCA<sup>7,8</sup> and grid-based density-peak based clustering<sup>9</sup> in MATLAB. The tutorial documents will also reference scripts that our group has developed, which are also included in the GitHub repository. 

•	Scripts that reference other software (for example, we have a script called chimScriptMaker.py that will use Chimera to build cyclic peptides) have a line (with comments) at the start of the code that will point to the directory that the software is installed, and this line should be updated according to where the user has installed the software.

### References
1. Yu, H.; Lin, Y.S. Toward structure prediction of cyclic peptides. <i>Phys. Chem. Chem. Phys.</i> <b>2015</b>, 17, 4210-4219.
2. McHugh, S.M.; Rogers, J.R.; Yu, H.; Lin, Y.-S. Insights into how cyclic peptides switch conformations. <i>J. Chem. Theory Comput.</i> <b>2016</b>, 12, 2480.
3. Abraham, M.J.; Murtola, T.; Schulz, R.; Páll, S.; Smith, J.C.; Hess, B.; Lindahl, E. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX <b>2015</b>, 1-2, 19.
4. Zhou, C.-Y.; Jiang, F.; Wu, Y.-D. Residue-specific force field based on protein coil library. RSFF2: modification of AMBER ff99SB. <i>J. Phys. Chem. B</i> <b>2015</b>, 119, 1035.
5. Piana, S.; Laio, A. A bias-exchange approach to protein folding. <i>J. Phys. Chem. B</i> <b>2007</b>, 111, 4553.
6. Laio, A.; Parrinello, M. Escaping free-energy minima. <i>Proc. Natl. Acad. Sci. U.S.A.</i> <b>2002</b>, 99, 12562.
7. Mu, Y.; Nguyen, P. H.; Stock, G. Energy landscape of a small peptide revealed by dihedral angle principal component analysis. <i>Proteins</i> <b>2005</b>, 58, 45. 
8. Sittel, F.; Jain, A.; Stock, G. Principal component analysis of molecular dynamics: on the use of cartesian vs. internal coordinates. <i>J. Chem. Phys.</i> <b>2014</b>, 141, 014111.
9. Rodriguez, A.; Laio, A. Clustering by fast search and find of density peaks. <i>Science</i> <b>2014</b>, 344, 1492.

### Contact
Please contact yu-shan.lin@tufts.edu if you have any questions.
