
Major Concerns:
---------------

*1. page 6; I found this notation to be needlessly confusing. Is it possible to give more explicit examples or change the forms? R^\pi means \pi(r_{ij}) is my primary problem.*

We agree. The notation was not only confusing, but also incorrect, because the operator \pi_{ij} was formally defined on binary distributions (D \propto (p,p')), while later we wrote expressions of the form \pi_{ij}(r_{ij}), where r_{ij} denotes only a scalar probability. We made it more precise and hopefully clearer. We abandoned the definition of \pi_{ij} completely and now we define the transformed matrix R^\pi through pairwise distributions:
(r^\pi_{ij},r^\pi_{ji}) \propto (r_{ij}\frac{P_i'}{P_i}, r_{ji}\frac{P_j'}{P_j}).

*2. page 15; tol=0.05 is not clearly justified. From the references figure, it should be 2^-4 or possibly higher.*

We have extended Figure 2 so that the plot now covers values of log_2(tol) up to −2, making the drop in accuracy for larger values of tol clearly visible. The choice of tol = 0.05 (used in subsequent experiments) corresponds approximately to 2^−4, and lies in the range where performance is stable across data sets. We chose this value as a representative point in this region; similar values (e.g., 2^−4) yield comparable results.

*3. In Figures 2 & 3, match the order of the data sets so that it is simpler to compare between the figures.*

We have reordered the data sets in Figure 3.

*4. Page 16, something is wrong with this sentence: “As table 2 shows, this simplex contained, in more than half cases, a base cover method matching the performance of W-L-Wang’s method.” (Commas added for clarity) Looking at this for a long time, I think the problem is that the last paragraph on page 15 should be phrased in a matching way, now it says slightly less than half of the cases have no good Bayes covariant method. Match the effect direction or put more details in somewhere on page 16 so that it’s not so hard to understand what you mean.*

We have added new sentences both on page 15 and page 16 to make the intention clearer.

*5. Page 17, why do table 2 and figure 3 disagree?*

The results in Table 2 and Figure 3 differ because they evaluate different subsets of the data using different evaluation metrics:

- They are evaluated on different datasets. Table 2 compares parameterless coupling methods on all available three-class subsets. Figure 3, however, focuses exclusively on a subset of 1000 problems where the normal method was already outperformed by Wu-Lin-Weng's method.

- They use different evaluation methodologies to account for parameter selection. Table 2 compares fixed, parameterless methods (e_3, s_1, s_2, s_3), which require no tuning, allowing for a direct accuracy comparison. Figure 3 evaluates a parametric family of methods. For a given dataset, usually there exists a parameter combination (a_1, a_2, a_3) that matches or outperforms Wu-Lin-Weng's method. However, identifying these optimal parameters requires access to the ground truth. To avoid this bias (data leakage), we used repeated five-fold cross-validation to select the parameters.

Therefore, Figure 3 does not show that the parametric family cannot represent a better solution, but rather that we are unable to reliably select those optimal parameters without seeing the ground truth. Because cross-validation effectively penalizes overfitting and therefore yields a more conservative estimate of performance, the overall accuracies cannot be directly compared one-to-one with Table 2.

While these methodological choices are mentioned in the text of Section 4.5, we agree that the shift in methodology between Table 2 and Figure 3 could be made more prominent. We have updated the last paragraph of section 4.5 to clarify this distinction by adding this sentence: "In contrast to the direct evaluations in Table 2, this performance drop reflects the practical difficulty of parameter selection and the generalization penalty inherent to cross-validation." We have also updated Figure 3 caption.

To improve clarity, we have also included a new figure showing the same experiment without cross-validation. It is placed before the original Figure 3; consequently, the new figure is labeled as Figure 3, and the original Figure 3 has been renumbered as Figure 4. The text in the affected section was updated accordingly. 

*6. Figure 4, it appears to me that WLW is consistently lower than radial, but this arrangement of bars makes it difficult to compare directly. In the text, you say that there is no uniform advantage, I disagree for this example.*

Figure 4 compares the difference between $\boldsymbol{c}_1$ and $\boldsymbol{c}_2$ for each coupling method, rather than directly comparing WLW and radial methods. We agree that the current visualization may give the impression of a consistent difference between methods, and we have clarified this point in the text. In particular, we included the following statement: "There appears to be some advantage to method $\boldsymbol{c}_2$, although it is not uniform."

Minor Comments:
---------------

*1. When we get to the first proof, we have not yet explained the relationship (in 2.3) between the K-1 sequential odds and the matrix of pairwise classifiers. At least mention that it is coming.*

We have added a paragraph right after the proof.

*2. Page 8, there’s a phrase that says “one could advantageously use other Bregman divergence derived from a proper score” something is wrong with his sentence, but I don’t know what it’s supposed to say so I’m not sure what the problem is. The way it is written, it sounds like there’s a type of divergent called “other Bregman divergence”*

We have revised the sentence to clarify that we are referring to the broader class of divergences. The revised text now reads: "...when we are primarily concerned with accuracy, one could advantageously employ alternative members of the Bregman divergence family derived from proper scoring rules"

*3. Page 10, I think the sentence should be: "they are not canonical, but with the health of ensemble, one can construct canonical coupling methods from the arboreal coupling methods"*

We have changed the sentence to: "They are not canonical, but canonical coupling methods can be constructed from arboreal methods using an ensemble construction."

*4. Table 2, You do not define e_3 as the normal method. Add this to the caption.*

We have clarified the caption of Table 2 and modified the text preceding (24) by replacing K>3 with K≥3, so that the normal coupling method also covers the case K=3 (i.e., e_3).

*5. Page 17, “only in minority of cases” should be “only in a minority of cases”*

We agree and we have updated the text.

Additional Comment:
-------------------

*The language is generally clear, but at the end of the article there are quite a few missing English articles that make the sentences difficult to read.*

We have identified several missing articles and added them.
