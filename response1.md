
*1. P. 3, line 5: ‘density’ with respect to which measure?*

Thank you for pointing out this issue. The choice of measure does not matter, 
but it is proper to emphasize it by explicitly mentioning a measure $\mu$. We changed
notation $d_i$ to $f_i \, d\mu$ to accomodate this change.

*2. P. 4, line 1: You call $\boldsymbol{p}$ ‘multiclass distribution’. As you have mentioned prior probabilities before,
perhaps it would be better to talk here about ‘posterior probabilities’ or ‘covariate-conditional
probabilities’?*

The point of the sentence is to introduce a new notation, which is used in the proofs
of both Proposition 1 and Proposition 4. The notation is not limited to posterior predictions. 
We slightly modified wording to make it clearer.

*3. New coupling method suggested by Eq. (5)*.

This is indeed the arboreal coupling described by a star graph. We added the explicit formula you proposed 
to Section 3.2.

*4. Eq. (9): You could refer to (9) as the ‘posterior correction’ formula of Saerens et al. (2002, Eq. (2.4)).*

Thank you, we have added the reference.

*5. P. 6, line 10: Please define clearly what you mean by ‘pairwise data R’.*

We defined $\boldsymbol{R}$ in text.

*6. P. 7, line 7: “this is the case of support vector machines”. Please provide a reference.*

We have added a reference. 

*7. P. 8, line 4: $div_{\textrm{KL}}(d, d′)$. Please specify what $d$ and $d′$ stand for (e.g. probability distributions on
finite sets).*

We added description of $d$ and $d'$.

*8. P. 9, line 12: Define ‘spanning tree’.*

We added a definition of a spanning tree and a reference to a book on graph theory.

*9. Proof of Proposition 2: Define the degree of a vertex.*

We included a definition of degree of a vertex.

*10. P. 10, line 5: “The predictions will generally span the vector space $R^K$.” Can you give an argument
why this statement is true?*

We were unable to provide a rigorous argument for this intuitively clear statement, so we have modified the definition of radial method. 
We express admiration for reviewer's acuity.

*11. Eq. (14): Misleading notation: $\hat \boldsymbol{p}$ is a function of $x$, $\hat p_i$ is a constant.*

We agree, and updated the notation.

*12. P. 11, line 7: “linear combination”. What about ‘logarithmic linear combination’ for clarity? *

We opted to alter the term to "weighted combination".

*13. P. 13, line 22: How do you understand “linear separability” in this context?*

We added the following wording: "As a preliminary step we conducted analysis of linear separability of pairs of classes in the datasets." We hope this makes the goal of our analysis unambiguous.

*14. P. 17, line 36: “Again we can see, that there is no uniform advantage to either method.” This
statement is somewhat surprising since Figure 4 appears to show a clear advantage of using “the
classification model where we apply (9) to binary classifiers” even if it is not uniform.*

Perhaps we were being too cautious. Now we updated our discussion to include the following statement: "There appears to be some advantage to method $\boldsymbol{c}_2$, although it is not uniform."

*15. P. 19, line 8: “One can investigate applicability of a Bayes covariant method.” This phrase looks
somewhat strange. Perhaps you want to say something like ‘About the deployment of a Bayes
covariant method should be decided on a case by case basis.’*

We updated the recommendation section.

*16. P. 21, line 23: Definition of Π′. I’m somewhat confused because Π′ appears to be a function of the
covariates (i.e. a function of x). If this observation is correct you should make it clear that this is
not a ‘common’ change of the priors.*

We agree and we have updated the text.

*17. P. 22, line 33: Is “rmj” the name of a journal?*

Thank you for the catch. We copied Google Scholar text without verification. We have updated the conference proceedings description.
