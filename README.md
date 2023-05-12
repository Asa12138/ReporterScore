## ReporterScore

ReporterScore Functional Enrichment Method for Microbiome

The original reporter-score method from Patil, K. R. et al. PNAS 2005. In this package, we provide two methods (mixed/directed) to get reporter-score for KEGG pathway or modules:

### Install

`install.packages("devtools")`\
`devtools::install_github('Asa12138/ReporterScore',dependencies=T)`\
`library(ReporterScore)`

### mixed

"mixed" mode is the original reporter-score method from Patil, K. R. et al. PNAS 2005.

In this mode, the reporter score is **non-directional**, and the larger the reporter score, the more significant the enrichment, but it cannot indicate the up-and-down regulation information of the pathway！(Liu, L. et al. iMeta 2023.)

steps: 1. Use the Wilcoxon rank sum test to obtain the P value of the significance of each KO difference between the two groups (ie $P_{koi}$, i represents a certain KO);

2.  Using an inverse normal distribution, convert the P value of each KO into a Z value ($Z_{koi}$), the formula:

$Z_{koi}=\theta ^{-1}(1-P_{koi})$

3.  "Upgrade" KO to pathway: $Z_{koi}$, calculate the Z value of the pathway, the formula:

$Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}$

where k means A total of k KOs were annotated to the corresponding pathway;

4.  Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of $Z_{pathway}$, the formula:

$Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k$

$μ_k$ is The mean of the random distribution, $σ_k$ is the standard deviation of the random distribution.

### directed

Instead, "directed" mode is a derived version of "mixed", referenced from https://github.com/wangpeng407/ReporterScore. This approach is based on the same assumption of many differential analysis methods: the expression of most genes has no significant change.

1.  Use the Wilcoxon rank sum test to obtain the P value of the significance of each KO difference between the two groups (ie $P_{koi}$, i represents a certain KO), and then divide the P value by 2, that is, the range of (0,1] becomes (0,0.5], $P_{koi}=P_{koi}/2$;

2.  Using an inverse normal distribution, convert the P value of each KO into a Z value ($Z_{koi}$), the formula:

$Z_{koi}=\theta ^{-1}(1-P_{koi})$

since the above P value is less than 0.5, all Z values will be greater than 0;

3.  Considering whether each KO is up-regulated or down-regulated, calculate $diff\_KO$,

$Z_{koi}=-Z_{koi}\ \ \ \ (diff\_KO<0)$

so $Z_{koi}$ is greater than 0 Up-regulation, $Z_{koi}$ less than 0 is down-regulation.

4.  "Upgrade" KO to pathway: $Z_{koi}$, calculate the Z value of the pathway,

$Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}$

where k means A total of k KOs were annotated to the corresponding pathway;

5.  Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of $Z_{pathway}$, the formula:

$Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k$

$μ_k$ is The mean of the random distribution, $σ_k$ is the standard deviation of the random distribution.

The finally obtained $Z_{adjustedpathway}$ is the Reporter score value enriched for each metabolic pathway.

In this mode, the Reporter score is non-directional, and a larger positive value represents a significant up-regulation enrichment, and a smaller Negative values represent significant down-regulation enrichment.

However, the disadvantage of this mode is that when a pathway contains about the same number of significantly up-regulates KOs and significantly down-regulates KOs, the final absolute value of Reporter score may approach 0, becoming a pathway that has not been significantly enriched.

## Reference
1. Patil, K. R. & Nielsen, J. Uncovering transcriptional regulation of metabolism by using metabolic network topology. Proc Natl Acad Sci U S A 102, 2685–2689 (2005). 

2. Liu, L., Zhu, R. & Wu, D. Misuse of reporter score in microbial enrichment analysis. iMeta n/a, e95.

3. https://github.com/wangpeng407/ReporterScore
