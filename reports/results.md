Results
Surveillance Coverage and Analytical Dataset
The final cleaned dataset contained 9,098 country-level pathogen-antibiotic surveillance observations from 103 countries and territories across 2020-2023.
Country participation increased over the study period:
- 2020: 81 countries/territories
- 2021: 83
- 2022: 88
- 2023: 95
  The total volume of reported interpretable antimicrobial susceptibility testing (AST) did not increase in parallel with country participation. Approximate annual totals were:
- 2020: 3.03 million interpretable AST results
- 2021: 3.56 million
- 2022: 3.59 million
- 2023: 3.48 million
  Thus, the number of participating countries increased across the study period, while reported testing volume increased initially and then remained broadly stable or declined slightly. This indicates that the composition of the surveillance dataset changed over time and supports the need for country-aware and stable-participation sensitivity analyses.
Descriptive Temporal Patterns
When all included pathogens, antibiotics, countries, and surveillance volumes were pooled, the proportion of reported AST results classified as resistant was:
- 18.4% in 2020
- 21.0% in 2021
- 20.5% in 2022
- 23.2% in 2023
  The pooled proportion was therefore higher in 2023 than in 2020, but the pattern was not monotonic because resistance decreased slightly between 2021 and 2022 before increasing again in 2023.
  These values are descriptive AST-level surveillance summaries. They should not be interpreted as estimates of global antimicrobial resistance prevalence because they combine different pathogens, antibiotics, countries, and surveillance volumes, and because country participation changed over time.
Pathogen-Specific Descriptive Patterns
Broad descriptive resistance patterns differed substantially between pathogens.
Acinetobacter spp. showed the highest pooled resistant-AST proportion across the study period, remaining at approximately one-half of reported AST results.
Klebsiella pneumoniae showed a sustained descriptive increase across the four years.
Escherichia coli remained at a lower pooled resistance level than K. pneumoniae and Acinetobacter spp., with a modest overall increase by 2023.
Salmonella spp. increased initially and then declined by 2023, while Streptococcus pneumoniae remained comparatively low and relatively stable.
Because these pathogen-level summaries combine multiple antibiotics, they are composition-dependent and were not treated as direct measures of a single biological resistance phenotype.
Pathogen-Antibiotic Heterogeneity
The 2023 pathogen-antibiotic analysis showed substantial heterogeneity both between pathogens and between antibiotics within the same pathogen.
This confirmed that broad pathogen-level averages can conceal important antibiotic-specific patterns and motivated pathogen-antibiotic-specific inferential modelling.
Forty-eight pathogen-antibiotic combinations were represented in the cleaned data. Of these, 37 met the predefined primary eligibility criteria of:
- representation in all four study years;
- data in both 2020 and 2023;
- at least 20 contributing countries in every study year.
  These 37 combinations formed the primary inferential analysis family.
Overall Adjusted Temporal Association
In the primary country-fixed-effect model adjusted for pathogen-antibiotic combination and using country-clustered robust standard errors, the estimated overall temporal association was close to the null:
- annual odds ratio (OR): approximately 1.01
- 95% CI: approximately 0.99-1.02
- p = 0.431
  Thus, the increase in the crude pooled AST-level resistance summary from 18.4% in 2020 to 23.2% in 2023 was not reproduced as a clear common temporal effect after accounting for country and pathogen-antibiotic composition.
  This contrast suggests that the crude pooled increase was influenced substantially by surveillance composition and weighting rather than representing a uniform increase across all adequately represented pathogen-antibiotic combinations.
Pathogen-Specific Adjusted Temporal Associations
After adjustment for antibiotic identity and country, pathogen-specific models showed the following annual associations:
- Acinetobacter spp.: OR 0.951, 95% CI 0.906-0.999, p = 0.0436, BH-FDR = 0.109
- Escherichia coli: OR 1.02, 95% CI 0.997-1.03, p = 0.105, BH-FDR = 0.175
- Klebsiella pneumoniae: OR 1.03, 95% CI approximately 1.00-1.06, p = 0.0282, BH-FDR = 0.109
- Salmonella spp.: OR 1.06, 95% CI 0.923-1.21, p = 0.408, BH-FDR = 0.408
- Streptococcus pneumoniae: OR 0.875, 95% CI 0.683-1.12, p = 0.284, BH-FDR = 0.355
  None of the five pathogen-specific associations remained below BH-FDR 0.05 after multiplicity adjustment.
  These results support the interpretation that temporal changes were concentrated in specific pathogen-antibiotic combinations rather than representing uniform pathogen-wide trends.
Primary Pathogen-Antibiotic Temporal Results
Across the 37 predefined pathogen-antibiotic models, two combinations met the BH-FDR < 0.05 criterion in the primary analysis:
1. Escherichia coli - meropenem
2. Klebsiella pneumoniae - levofloxacin
Escherichia coli - Meropenem
The primary model estimated:
- OR per year: 1.13
- 95% CI: 1.07-1.19
- p = 0.0000599
- BH-FDR = 0.00222
- implied OR for 2023 versus 2020: 1.43
- 95% CI for 2023 versus 2020: 1.21-1.70
  The pooled AST-weighted resistance percentage increased from approximately 0.47% in 2020 to 1.77% in 2023, an absolute increase of approximately 1.30 percentage points.
  The relative increase was therefore substantial despite the low absolute resistance percentage.
Klebsiella pneumoniae - Levofloxacin
The primary model estimated:
- OR per year: 1.09
- 95% CI: 1.03-1.15
- p = 0.00213
- BH-FDR = 0.0395
- implied OR for 2023 versus 2020: 1.28
- 95% CI for 2023 versus 2020: 1.10-1.50
  The pooled AST-weighted resistance percentage increased from approximately 12.3% in 2020 to 17.5% in 2023, an absolute increase of approximately 5.23 percentage points.
  This combination therefore showed both a consistent relative increase and a larger absolute change than E. coli-meropenem.
Ordinary Binomial Mixed-Effects Analysis
The ordinary binomial mixed-effects analysis produced point estimates that were highly consistent in direction with the primary models.
Across all 37 pathogen-antibiotic combinations:
- direction agreed between the primary fixed-effect models and the binomial GLMMs for 37/37 combinations;
- 2/37 combinations had BH-FDR < 0.05 in the primary analysis;
- 30/37 combinations had BH-FDR < 0.05 in the ordinary binomial GLMM analysis;
- the same two primary signals were below BH-FDR 0.05 in both approaches.
  However, the ordinary binomial GLMMs produced substantially narrower confidence intervals and widespread statistical significance despite persistent residual overdispersion.
  The overall adjusted GLMM estimated an annual log-odds coefficient of approximately 0.00602, essentially identical to the primary overall model, but with a much smaller standard error and p < 10^-14. Pearson dispersion for this model was approximately 67.7.
  Among the 37 pathogen-antibiotic GLMMs, median Pearson dispersion was approximately 3.83, with a maximum above 10.
  These results indicate that the main discrepancy between the two approaches was uncertainty rather than direction or point-estimate magnitude. Because substantial extra-binomial variation remained, nominal p-values from the ordinary binomial GLMMs were not treated as the primary inferential evidence.
Beta-Binomial Sensitivity Analysis
All 37 beta-binomial mixed models fitted successfully with positive-definite Hessians.
Eleven of the 37 combinations had BH-FDR < 0.05 under the beta-binomial specification.
For the two primary signals:
E. coli - meropenem
- OR per year: 1.12
- 95% CI: 1.04-1.21
- BH-FDR = 0.0150
Klebsiella pneumoniae - Levofloxacin
- OR per year: 1.11
- 95% CI: 1.04-1.19
- BH-FDR = 0.00772
  Thus, both primary signals remained positive and of similar magnitude after explicitly modelling extra-binomial variation.
Stable-Country Sensitivity Analysis
For each pathogen-antibiotic combination, countries were restricted to those reporting that specific combination in all four years.
The number of stable countries ranged from 10 to 64 across the 37 combinations, with a median of 37. Thirty-two combinations retained at least 20 stable countries.
Only one of the 37 stable-country models had BH-FDR < 0.05.
For the two primary signals:
E. coli - meropenem
- 54 stable countries
- OR per year: 1.13
- 95% CI: 1.07-1.19
- BH-FDR = 0.000631
Klebsiella pneumoniae - Levofloxacin
- 30 stable countries
- OR per year: 1.08
- 95% CI: 1.03-1.14
- BH-FDR = 0.0705
  The E. coli-meropenem signal therefore remained strongly supported after restricting the analysis to consistently reporting countries.
  For K. pneumoniae-levofloxacin, the estimated association remained positive and of similar magnitude, but the BH-adjusted p-value moved above 0.05.
Categorical-Year Sensitivity Analysis
Treating year as a categorical variable removed the assumption of a constant linear annual change in log odds.
Three of the 37 combinations had BH-FDR < 0.05 for the direct 2023-versus-2020 comparison.
For the primary signals:
E. coli - meropenem
- OR for 2023 versus 2020: 1.55
- 95% CI: 1.24-1.94
- BH-FDR = 0.00725
Klebsiella pneumoniae - Levofloxacin
- OR for 2023 versus 2020: 1.30
- 95% CI: 1.12-1.51
- BH-FDR = 0.0119
  Both combinations therefore remained positively associated with calendar time when the linear-year assumption was relaxed.
Equal-Country-Weighted Descriptive Sensitivity
Equal-country weighting was used to examine whether the direction of change depended on high-volume surveillance systems dominating pooled AST-level estimates.
For E. coli-meropenem:
- equal-country-weighted 2020-to-2023 change: +2.10 percentage points
- stable-country equal-weight change: +1.31 percentage points
  For K. pneumoniae-levofloxacin:
- equal-country-weighted 2020-to-2023 change: +5.93 percentage points
- stable-country equal-weight change: +7.39 percentage points
  Both primary signals therefore retained an increasing direction when each country was given equal descriptive weight.
Cross-Sensitivity Consistency
Across all 37 primary pathogen-antibiotic combinations, 21 showed the same direction of change across all major sensitivity approaches.
The number of combinations meeting BH-FDR < 0.05 varied substantially by analytical specification:
- primary country-fixed-effect robust analysis: 2
- beta-binomial mixed models: 11
- stable-country analysis: 1
- categorical-year analysis: 3
  This variability indicates that many apparent temporal associations were sensitive to model specification, surveillance composition, or weighting.
  In contrast, the two primary signals retained the same increasing direction across all major sensitivity analyses.
Overall Interpretation of the Primary Signals
Escherichia coli-meropenem showed the most consistently supported temporal increase in the project.
Its direction and approximate magnitude remained similar under:
- the primary country-fixed-effect model;
- beta-binomial modelling;
- restriction to countries reporting in all four years;
- categorical-year modelling;
- equal-country weighting;
- stable-country equal weighting.
  Klebsiella pneumoniae-levofloxacin also showed a consistent increasing direction and similar effect magnitude across analyses. However, multiplicity-adjusted evidence was weaker in the stable-country analysis, despite the stable-country odds ratio remaining above 1 with a 95% confidence interval excluding 1.
  The results therefore support a more cautious distinction between the two signals: E. coli-meropenem was the most consistently supported increase, while K. pneumoniae-levofloxacin showed a robust directional increase with somewhat greater sensitivity of the multiplicity-adjusted evidence to analytical specification.
  Neither result should be interpreted as a causal effect or as a population-representative global trend. Both describe temporal associations within the countries and surveillance systems contributing relevant WHO GLASS bloodstream-infection data during 2020-2023.