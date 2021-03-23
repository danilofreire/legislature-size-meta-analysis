# TODO Appendix

1. Section A: Change the numbers to add new matched files. **To Catarina and Huzeyfe**

2. Section B: Check if automated code is getting the article numbers rights. Should be 31 now. **To Catarina and Huzeyfe**

3. Section B: Change the Flow Chart to the new numbers of selected, and change the numbers in each step. **To Catarina and Huzeyfe**

4. Section C: check if the numbers are correct. **To Catarina and Huzeyfe**

5. Section C: describe new approach using three-level models. **DONE**

6. Section G1: add citations supporting why REML is best for continuous data. **DONE**

7. Section G1: Explain three-level meta-analysis. **DONE**

8. Section G1: Explain Funnel Plots and publication bias. **DONE**

9. Sections G and I: Names of the models in the graphs, regressions, and plots are inconsistently used. **DONE**

10. Sections G and I: Check if all funnel plots were added. **DONE**

11. Section I: Remove all the funnel plots and funnel tests. Preferably: think what to do with this section. **DONE**

12. Section H: Remove the section and make the code not echo. **DONE**

13. Fix meta-regressions to add the multilevel estimation. **DONE**

14. Remove permutation tests. No permutation tests for the multilevel meta regressions. **DONE**

15. Change the interpretation of the meta-regressions. **DONE**

16. Do a final sweep to check if it is all ok. **DONE**

17. Add the model section with info about weingast. **DONE**

18. Check others such as Primo, Lee, and others. **DONE**

19. Catarina points:

```
Line 1107: Should be column = 3, not 1.

Paper 542: Institutional Design variable should be coded as Mixed; it is currently coded as Bicameral (This was my mistake, sorry guys). Must change to Mixed.

Lines 1283-1314 : I think the coefficient is wrong (I am actually not sure this coefficient shows up anywhere in the tables...).

Line 1479: SE is wrong. It is supposed to be 0.295/7.545

Lines 2307-2338: I think the coefficient is wrong (I am actually not sure this coefficient shows up anywhere in the tables...).

Line 2790: Should be table = 3, not 4.

Line 3088: Should be column = 2, not 1.

Line 3288: SE is wrong. It is supposed to be 0.1594/3.05

Line 3360: SE is wrong. It is supposed to be 0.039 (according to the table).

Line 3489: Add 0 before . (it is not being considered as numeric, I think).

--

Paper 165: Coefficient using N in Table 4, column 7, is in the dataset twice. Delete one of the occurences.

Paper A258: Coefficient in Table 4, column 4, is in the dataset twice. Delete one of the occurences.

Paper 230: I think the Unicameral Chamber Size data should maybe be included. Here is the info: 
table = 2,
column = 1,
coef = 0.0107,
se = 0.0033

--

Paper 288: This was the major mistake. In the dataset, we coded the indvar as N, when it is actually log(N). Please check the note below Table 3 of the attached .pdf just to make sure I read correctly.
We need to replace the indvar2 variable with the logN value and re-run everything.
```