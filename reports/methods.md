Methods
Study Design and Data Source
We conducted a retrospective analysis of antimicrobial resistance surveillance data reported through the World Health Organization Global Antimicrobial Resistance and Use Surveillance System (WHO GLASS). The analysis focused on bloodstream-infection surveillance observations reported between 2020 and 2023.
The study was designed to characterize temporal and pathogen-antibiotic-specific patterns in reported antimicrobial resistance among countries and territories contributing relevant data to GLASS. The objective was not to estimate population-representative global antimicrobial resistance prevalence, because surveillance participation, testing intensity, and reporting practices differed between countries and over time.
The final analysis included five pathogen groups:
- Acinetobacter spp.
- Escherichia coli
- Klebsiella pneumoniae
- Salmonella spp.
- Streptococcus pneumoniae
  Pseudomonas aeruginosa was considered during project development but was excluded from the final analytical scope.
Data Import and Processing
WHO GLASS source files were imported and processed in R.
Because some country and territory names contained embedded commas, including names such as "occupied Palestinian territory, including east Jerusalem," standard comma-delimited parsing did not reliably preserve row structure. A custom parsing procedure was therefore used to identify the relevant table section and reconstruct records while preserving country names containing commas.
The cleaning workflow validated study year, pathogen, antibiotic, country identifier, and resistance-count fields. Resistant counts were checked to ensure that they were non-negative and did not exceed the corresponding number of interpretable antimicrobial susceptibility tests (AST).
The reported WHO resistance percentage was independently recalculated as:
Resistance percentage = Resistant / InterpretableAST × 100
and compared with the percentage supplied in the source data as a quality-control check.
The final cleaned dataset contained:
- 9,098 country-level pathogen-antibiotic surveillance observations;
- 103 countries and territories;
- four study years (2020-2023);
- five included pathogen groups.
Unit of Analysis and Outcome Definition
The analytical unit was a country-year-pathogen-antibiotic surveillance observation.
For each observation, the dataset included the number of interpretable AST results (InterpretableAST) and the number classified as resistant (Resistant). The number classified as non-resistant was derived as:
NonResistant = InterpretableAST − Resistant
Grouped-binomial models were fitted using resistant and non-resistant counts rather than modelling reported resistance percentages as continuous outcomes.
InterpretableAST was treated as a count of reported susceptibility-test results and was not assumed to represent unique bacterial isolates or unique patients. An individual isolate may have contributed susceptibility results for more than one antibiotic. Consequently, AST counts summed across antibiotics were not interpreted as numbers of unique infections or isolates.
Descriptive Analyses
Surveillance coverage was characterized before temporal resistance patterns were interpreted. Descriptive analyses included:
- number of contributing countries and territories per year;
- total interpretable AST volume;
- participation consistency;
- pathogen-specific surveillance coverage;
- WHO-region-specific surveillance coverage;
- pooled AST-level resistance by year;
- pooled resistance patterns by pathogen;
- pathogen-antibiotic resistance patterns;
- descriptive 2020-to-2023 changes for pathogen-antibiotic combinations.
  Pooled AST-level resistance was calculated as:
  Pooled resistance percentage = ΣResistant / ΣInterpretableAST × 100
  These pooled estimates are AST-volume weighted. Countries contributing larger numbers of AST results therefore exert greater influence on pooled estimates.
  Descriptive estimates were interpreted as patterns within reported GLASS surveillance data and not as population-representative estimates of worldwide antimicrobial resistance prevalence.
  A 2023 pathogen-antibiotic heatmap was also produced to illustrate heterogeneity across specific resistance phenotypes and to demonstrate the limitations of pooling results across antibiotics within a pathogen.
Eligibility for Pathogen-Antibiotic Inferential Analyses
Coverage of all pathogen-antibiotic combinations was evaluated before fitting the primary models. Forty-eight distinct combinations were identified in the cleaned dataset.
For the primary inferential analysis, a pathogen-antibiotic combination was retained if it:
- was represented in all four study years;
- contained observations in both 2020 and 2023;
- included at least 20 contributing countries or territories in every study year.
  The minimum-country criterion was selected during analytical development after examining the distribution of country coverage. Its purpose was to avoid drawing clustered inferential conclusions from combinations supported by only a small number of country clusters in one or more years. It was not treated as a universal statistical threshold or WHO-recommended cutoff.
  Thirty-seven pathogen-antibiotic combinations met the final primary eligibility criteria.
Primary Statistical Analysis
Year was modelled as a continuous variable and centred at 2020, such that:
- 2020 = 0
- 2021 = 1
- 2022 = 2
- 2023 = 3
  The temporal coefficient therefore represented the estimated change in log odds of reported resistance per additional year.
  The primary pathogen-antibiotic-specific models used grouped-binomial logistic regression with country fixed effects:
  logit(p_it) = b0 + b1(Year_t) + a_i
  where p_it represents the probability that a reported AST result was classified as resistant for a given country and year, b1 represents the annual temporal association, and a_i represents country-specific fixed effects.
  Standard errors were estimated using country-clustered robust covariance estimators to account for repeated observations and within-country dependence. Small-sample inference used the number of country clusters to determine the reference degrees of freedom.
  For each model, the annual odds ratio was calculated as:
  OR_year = exp(b1)
  The implied 2023-versus-2020 odds ratio under the continuous-year specification was calculated as:
  OR_2023:2020 = exp(3 × b1) = OR_year^3
  Pearson dispersion statistics were calculated as diagnostics for extra-binomial variation.
Overall and Pathogen-Specific Models
A broader overall model was fitted across the 37 eligible pathogen-antibiotic combinations. This model included:
- continuous year;
- pathogen-antibiotic combination fixed effects;
- country fixed effects;
- country-clustered robust standard errors.
  The model estimated a common temporal association while reducing confounding by large baseline differences between pathogen-antibiotic combinations. Because it pooled multiple pathogens and antibiotics, remained influenced by surveillance volume, and imposed a common temporal effect across heterogeneous resistance phenotypes, it was treated as a broad summary rather than as an estimate of global antimicrobial resistance change.
  Separate pathogen-specific models were also fitted. Within each pathogen, models included:
- continuous year;
- antibiotic identity;
- country fixed effects;
- country-clustered robust standard errors.
  Adjustment for antibiotic identity reduced the possibility that temporal changes in the composition of antibiotics tested would be mistaken for pathogen-level resistance changes.
Multiple Testing
The primary inferential family consisted of the 37 eligible pathogen-antibiotic-specific models.
To account for multiple comparisons, p-values were adjusted using the Benjamini-Hochberg false discovery rate procedure.
Interpretation did not rely on adjusted p-values alone. Effect magnitude, 95% confidence intervals, absolute resistance patterns, surveillance coverage, model diagnostics, and consistency across sensitivity analyses were considered together.
Binomial Mixed-Effects Robustness Analysis
To examine sensitivity to alternative treatment of country-level heterogeneity, grouped-binomial generalized linear mixed-effects models were fitted with country-specific random intercepts.
For individual pathogen-antibiotic combinations, the model took the form:
logit(p_it) = b0 + b1(Year_t) + u_i
where u_i represents a country-level random intercept.
Analogous adjusted mixed-effects models were fitted at the overall and pathogen-specific levels.
Model convergence, singularity, and Pearson dispersion were assessed. The mixed-effects models were treated as a robustness and diagnostic analysis rather than as replacements for the primary models because substantial residual extra-binomial variation remained under the ordinary binomial distribution.
Beta-Binomial Sensitivity Analysis
Because substantial overdispersion remained in the ordinary binomial mixed-effects models, pathogen-antibiotic temporal associations were additionally estimated using beta-binomial mixed models with country-level random intercepts.
The beta-binomial distribution allows the underlying probability of resistance to vary more than expected under a simple binomial sampling model and therefore provides an explicit model for extra-binomial heterogeneity.
Beta-binomial models were evaluated using optimizer convergence and Hessian diagnostics. Benjamini-Hochberg adjustment was again applied across the 37 primary combinations within this sensitivity analysis.
Stable-Country Sensitivity Analysis
To evaluate whether changing surveillance participation influenced temporal results, a stable-country sensitivity analysis was performed separately for each pathogen-antibiotic combination.
Countries were retained only if they contributed that specific pathogen-antibiotic combination in all four study years. The primary fixed-effect grouped-binomial model with country-clustered robust standard errors was then refitted within the restricted dataset.
Because the stable set was defined separately for each pathogen-antibiotic combination, the number of countries included differed between models.
Categorical-Year Sensitivity Analysis
The primary continuous-year model assumes a constant linear change in log odds between 2020 and 2023.
To relax this assumption, year was additionally modelled as a categorical variable with 2020 as the reference year. The 2023 coefficient therefore provided a direct comparison of 2023 with 2020 without constraining the intermediate 2021 and 2022 values to follow a linear log-odds trajectory.
Country fixed effects and country-clustered robust standard errors were retained.
Equal-Country-Weighted Descriptive Sensitivity Analysis
Primary pooled summaries weight countries according to their reported AST volume.
To evaluate the influence of this weighting, country-level resistance percentages were first calculated separately for each country, pathogen-antibiotic combination, and year. These country-specific percentages were then averaged with equal weight assigned to each country.
A second equal-weight analysis was restricted to countries reporting the relevant pathogen-antibiotic combination in all four study years.
These summaries were treated as descriptive sensitivity analyses because they target a different estimand from AST-volume-weighted pooled resistance. Specifically, they describe the average resistance percentage among contributing countries rather than the proportion resistant among all reported AST results.
Assessment of Robustness
Sensitivity analyses were interpreted jointly rather than as competing approaches from which the most statistically significant model could be selected.
Greater confidence was placed in a temporal pattern when its direction and approximate magnitude remained similar across:
- the primary country-fixed-effect model with country-clustered robust inference;
- the beta-binomial mixed model;
- stable-country restriction;
- categorical-year modelling;
- equal-country-weighted descriptive analyses.
  Differences between methods were treated as evidence about sensitivity to surveillance composition, variance assumptions, or model specification rather than as reasons to preferentially select one result based on statistical significance.
Software and Reproducibility
All data preparation, descriptive analysis, statistical modelling, sensitivity analyses, and figure generation were conducted in R.
The analytical workflow was separated into reproducible scripts covering:
- package setup;
- data inspection;
- data cleaning;
- exploratory analysis;
- primary statistical analysis;
- mixed-effects modelling;
- sensitivity analyses;
- final figure generation.
  Key packages included:
- tidyverse and here for data processing and project paths;
- sandwich and lmtest for country-clustered robust inference;
- lme4 for binomial mixed-effects models;
- glmmTMB for beta-binomial mixed modelling.
  Intermediate analytical tables were exported to the project outputs/ directory, while final figures were generated programmatically from the corresponding processed data and model outputs.
Interpretation Framework
All estimates were interpreted as patterns observed in reported WHO GLASS surveillance data.
The analysis did not treat pooled resistance percentages as population-representative estimates of global antimicrobial resistance prevalence and did not interpret observed temporal associations as causal effects.
The final interpretation considered:
- effect size;
- uncertainty;
- multiple-testing correction;
- surveillance coverage;
- model diagnostics;
- sensitivity to analytical assumptions;
- consistency between relative and absolute patterns;
- limitations of the surveillance data.