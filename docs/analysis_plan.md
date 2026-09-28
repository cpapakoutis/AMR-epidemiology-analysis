WHO GLASS Antimicrobial Resistance Analysis
1. Project Aim
This project investigates temporal, geographic, pathogen-specific, and antibiotic-specific patterns in antimicrobial resistance reported through the WHO Global Antimicrobial Resistance and Use Surveillance System (GLASS).
The analysis focuses on bloodstream-infection surveillance data reported between 2020 and 2023.
The aim is not to estimate a single population-representative measure of "global antimicrobial resistance." Instead, the project examines patterns observed among countries and territories contributing relevant data to WHO GLASS and evaluates how robust observed temporal patterns are to different analytical assumptions.
2. Primary Research Question
How did reported antimicrobial resistance patterns change between 2020 and 2023 among countries and territories contributing bloodstream-infection surveillance data to WHO GLASS?
3. Secondary Research Questions
Temporal patterns
How did observed antimicrobial resistance change between 2020 and 2023?
Pathogen-specific patterns
Did temporal resistance patterns differ between bacterial pathogens?
Antibiotic-specific patterns
Which sufficiently represented pathogen-antibiotic combinations showed evidence of increasing or decreasing resistance over time?
Geographic and surveillance patterns
How did country participation and surveillance coverage vary across years and WHO regions?
Robustness of temporal patterns
Did conclusions about temporal change remain similar when alternative analytical assumptions were used, including:
- country fixed effects with country-clustered robust inference;
- country-level random intercepts;
- beta-binomial mixed models;
- categorical rather than continuous year effects;
- restriction to countries contributing a given pathogen-antibiotic combination in all four years;
- descriptive analyses giving countries equal rather than AST-volume-based weight?
4. Study Scope
Data source
WHO Global Antimicrobial Resistance and Use Surveillance System (GLASS).
Years
2020-2023.
Infection type
Bloodstream infections.
Geographic scope
Countries and territories contributing relevant observations to WHO GLASS during the study period.
Unit of analysis
The cleaned analytical dataset contains country-level pathogen-antibiotic surveillance observations.
Each observation includes information such as:
- country or territory;
- WHO region;
- year;
- pathogen;
- antibiotic;
- number of interpretable antimicrobial susceptibility tests (InterpretableAST);
- number of reported resistant AST results (Resistant);
- resistance percentage.
InterpretableAST is treated as a count of reported susceptibility-test results. It should not automatically be interpreted as a count of unique bacterial isolates or unique patients, because the available data do not establish one-to-one correspondence at those levels.
The final cleaned dataset contains:
- 9,098 analytical observations;
- 103 countries and territories;
- 4 years;
- 5 included pathogens.
5. Selected Pathogens
The final analysis includes:
- Acinetobacter spp.
- Escherichia coli
- Klebsiella pneumoniae
- Salmonella spp.
- Streptococcus pneumoniae
Pseudomonas aeruginosa was considered during project development but was excluded from the final project scope.
6. Main Outcome
The main outcome is reported antimicrobial resistance.
For each country-pathogen-antibiotic observation, resistance is represented using:
- the number of resistant AST results;
- the number of interpretable AST results;
- the corresponding resistance percentage.
For grouped-binomial statistical models:
NonResistant = InterpretableAST − Resistant
Models therefore use resistant and non-resistant counts rather than treating reported resistance percentages as independent continuous observations.
Because the same bacterial isolate may be tested against more than one antibiotic, summed AST counts across antibiotics should not be interpreted as numbers of unique infections or unique isolates.
7. Descriptive Analysis
The descriptive analysis examined:
- country participation over time;
- total reported interpretable AST volume over time;
- pooled AST-level resistance by year;
- pooled resistance patterns by pathogen;
- pathogen-antibiotic resistance patterns;
- surveillance coverage by pathogen and year;
- surveillance coverage by WHO region and year;
- country participation consistency;
- descriptive 2020-to-2023 changes for pathogen-antibiotic combinations.
For pooled summaries:
Pooled resistance (%) = total resistant AST results / total interpretable AST results × 100
These estimates are AST-volume weighted. Countries contributing larger numbers of AST results therefore have greater influence.
Descriptive estimates are interpreted as patterns within reported GLASS surveillance data and not as population-representative estimates of worldwide AMR prevalence.
8. Primary Inferential Analysis
The primary inferential analysis evaluated temporal changes while accounting for repeated observations contributed by the same countries.
8.1 Continuous year
Year was centred at 2020:
- 2020 = 0
- 2021 = 1
- 2022 = 2
- 2023 = 3
The main continuous-year coefficient therefore represents the estimated change in log odds of reported resistance per additional year.
8.2 Country adjustment
The primary models used:
- country fixed effects, to account for time-invariant differences in baseline resistance between countries;
- country-clustered robust standard errors, to account for within-country dependence and variance misspecification.
8.3 Overall adjusted model
The broad overall model included:
- continuous year;
- pathogen-antibiotic combination fixed effects;
- country fixed effects;
- country-clustered robust standard errors.
This model estimated a common temporal association across adequately represented pathogen-antibiotic combinations.
Because it pools multiple pathogens and antibiotics and remains influenced by surveillance volume, it is treated as a broad summary rather than as a single estimate of global AMR change.
8.4 Pathogen-specific models
Within each pathogen, models included:
- continuous year;
- antibiotic identity;
- country fixed effects;
- country-clustered robust standard errors.
Adjustment for antibiotic identity reduces the risk that changes in which antibiotics were tested are interpreted as changes in resistance.
8.5 Pathogen-antibiotic-specific models
Separate models were fitted for each eligible pathogen-antibiotic combination:
resistance ~ year + country fixed effects
with standard errors clustered by country.
These models form the primary inferential analysis family.
9. Eligibility for Pathogen-Antibiotic Models
A pathogen-antibiotic combination was eligible for the primary inferential analysis if it:
- was represented in all four study years;
- had observations in both 2020 and 2023;
- had at least 20 contributing countries in every year.
The minimum-country criterion was chosen pragmatically to avoid drawing clustered inferential conclusions from combinations with very sparse country representation.
Of 48 pathogen-antibiotic combinations in the cleaned data:
- 37 met the primary eligibility criteria.
A less restrictive threshold of at least 10 countries per year was considered during sensitivity planning but was not used to define the primary analysis family.
10. Multiple Testing
The primary pathogen-antibiotic analysis involved 37 statistical tests.
The Benjamini-Hochberg false discovery rate (BH-FDR) procedure was applied across these 37 predefined models.
Interpretation considered:
- odds ratios;
- 95% confidence intervals;
- absolute resistance patterns;
- BH-adjusted p-values;
- surveillance coverage;
- model diagnostics;
- consistency across sensitivity analyses;
- epidemiological plausibility.
Crossing an FDR threshold was not treated as sufficient evidence on its own.
11. Mixed-Effects Robustness Analysis
A separate binomial generalized linear mixed-effects analysis was conducted using country-level random intercepts.
Overall model
resistance ~ year + pathogen-antibiotic combination + (1 | country)
Pathogen-specific models
resistance ~ year + antibiotic + (1 | country)
Pathogen-antibiotic models
resistance ~ year + (1 | country)
These models were used as a robustness and diagnostic analysis rather than as the primary inferential approach.
The random-intercept GLMMs produced temporal point estimates very similar in direction to the primary country-fixed-effect models, but substantially narrower standard errors and much more widespread statistical significance.
Substantial residual overdispersion remained in many models, indicating that the ordinary binomial variance assumption did not adequately represent the heterogeneity in the surveillance data.
The nominal p-values from these ordinary binomial GLMMs were therefore not used as the principal inferential evidence.
12. Sensitivity Analyses
Four major sensitivity analyses were conducted.
12.1 Beta-binomial mixed models
Because the ordinary binomial GLMMs showed substantial extra-binomial variation, pathogen-antibiotic models were re-estimated using beta-binomial mixed models with country random intercepts.
The beta-binomial distribution allows greater between-observation variability than an ordinary binomial model.
All 37 primary combinations were successfully modelled, with acceptable Hessian diagnostics.
12.2 Stable-country analysis
For each pathogen-antibiotic combination, countries were restricted to those reporting that combination in all four study years.
The same general country-fixed-effect and country-clustered robust modelling approach was then repeated.
This analysis directly addresses the possibility that temporal patterns are driven by countries entering or leaving the surveillance dataset.
12.3 Categorical-year analysis
Year was modelled as a categorical variable with 2020 as the reference year.
The 2023 coefficient therefore estimates the direct 2023-versus-2020 comparison without requiring resistance to follow a constant linear trend on the log-odds scale through 2021 and 2022.
12.4 Equal-country-weighted descriptive analysis
For each pathogen-antibiotic combination and year, resistance was first calculated separately within each contributing country.
Country-level percentages were then averaged without weighting by AST volume.
A second version restricted the calculation to countries contributing that combination in all four years.
These analyses answer a different descriptive question from the pooled AST-weighted summaries:
What happened to the average reported country-level resistance percentage rather than to the pooled set of reported AST results?

13. Interpretation of Sensitivity Analyses
Sensitivity analyses were used to assess whether substantive conclusions were robust to reasonable analytical choices.
A finding was considered more convincing when:
- the direction remained consistent;
- effect estimates were of similar magnitude;
- confidence intervals remained compatible with the same substantive conclusion;
- results were not dependent solely on changing country participation;
- the pattern remained visible under alternative weighting or variance assumptions.
Sensitivity analyses were not used to search for whichever model produced the smallest p-value.
14. Primary Temporal Signals
Two pathogen-antibiotic combinations met the predefined BH-FDR < 0.05 criterion in the primary analysis:
- Escherichia coli - meropenem
- Klebsiella pneumoniae - levofloxacin
Escherichia coli - meropenem
The primary model estimated an odds ratio of approximately 1.13 per year.
The signal remained directionally consistent under:
- beta-binomial modelling;
- stable-country restriction;
- categorical-year analysis;
- equal-country weighting;
- equal-country weighting among stable reporters.
This combination showed the most consistently supported temporal increase across the primary and sensitivity analyses.
Klebsiella pneumoniae - levofloxacin
The primary model estimated an odds ratio of approximately 1.09 per year.
The increasing direction and approximate magnitude were also consistent across the sensitivity analyses.
Multiplicity-adjusted evidence was somewhat weaker in the stable-country analysis, although the estimated association remained positive.
These findings are interpreted as temporal associations in the reported surveillance data, not as causal effects or population-representative estimates of global resistance change.
15. Epidemiological Interpretation
Potential explanations for observed resistance patterns may include:
- antimicrobial selection pressure;
- antibiotic consumption;
- healthcare-associated infections;
- healthcare access;
- surveillance-system differences;
- diagnostic practices;
- laboratory capacity;
- changing country participation;
- COVID-19-related disruption;
- broader One Health factors.
These factors are not directly tested by the current WHO GLASS dataset.
They may therefore be discussed as possible contextual explanations, but they should not be presented as causal conclusions from this analysis unless additional external data are incorporated.
16. Key Limitations
The analysis must be interpreted in light of:
- changing country participation between years;
- non-random surveillance coverage;
- unequal numbers of AST results contributed by countries;
- variation in laboratory and surveillance capacity;
- differences in pathogen and antibiotic composition;
- repeated observations within countries;
- extra-binomial heterogeneity;
- inability to identify unique patients or isolates from the available aggregate AST counts;
- possible repeated contribution of the same isolate across multiple antibiotics;
- incomplete comparability of surveillance systems between countries;
- the short four-year study period;
- multiple statistical testing;
- the observational nature of surveillance data;
- lack of direct measurement of many potential explanatory factors.
Country adjustment and sensitivity analyses reduce some concerns but cannot eliminate these limitations.
17. Descriptive Versus Inferential Evidence
The project explicitly distinguishes between:
Descriptive evidence
What patterns are visible in the reported surveillance records?
Examples include:
- changing country participation;
- pooled AST resistance;
- pathogen-level resistance patterns;
- annual pathogen-antibiotic percentages;
- equal-country-weighted resistance summaries.
Inferential evidence
How strongly do statistical models support temporal associations after accounting for relevant structure and uncertainty?
Examples include:
- country-fixed-effect models;
- country-clustered robust standard errors;
- BH-FDR correction;
- beta-binomial models;
- stable-country analyses;
- categorical-year analyses.
Observed associations are not automatically interpreted as causal effects.
18. Guiding Analytical Principle
The objective of the project is not to obtain statistically significant results.
The objective is to determine:
- what conclusions are reasonably supported by the surveillance data;
- how strongly those conclusions depend on analytical assumptions;
- whether important findings persist under reasonable sensitivity analyses;
- what uncertainty and limitations must accompany their interpretation.
Model complexity, p-values, and statistical significance are treated as tools for inference rather than as goals in themselves.