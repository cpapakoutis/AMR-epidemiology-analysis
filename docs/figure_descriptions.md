
Figure Descriptions - WHO GLASS Antimicrobial Resistance Analysis
This document is the central figure catalogue for the project. It records what each figure asks, how it was calculated, what it shows, how it should be interpreted, and the main limitations that should accompany it.

Figure 1 - Country participation in the analysed WHO GLASS data
File: figures/01_country_participation_over_time.png

Question Addressed
How did the number of countries and territories contributing relevant bloodstream-infection surveillance data change between 2020 and 2023?

Data and Calculation
For each year, the cleaned WHO GLASS dataset was grouped by Year, and the number of unique country/territory ISO3 codes was counted using n_distinct(Iso3).
Observed participation:

2020: 81 countries/territories
2021: 83
2022: 88
2023: 95
What the Figure Shows
The number of participating countries/territories increased across the study period, with the largest increase occurring between 2022 and 2023.

Interpretation
The surveillance population is not identical across years. An apparent temporal change in pooled resistance could therefore reflect both:

genuine changes in resistance among reporting surveillance systems, and
changes in which countries/territories contributed data.
This is one reason later analyses include country-aware modelling and a sensitivity analysis restricted to countries contributing data across all four years.
Main Caveats
A higher number of participating countries does not necessarily mean more complete surveillance within each country.
Countries differ substantially in surveillance intensity, laboratory capacity, population size, healthcare systems, and numbers of AST results contributed.
Country participation alone does not describe surveillance volume.
Portfolio / Manuscript Caption
Figure 1. Country participation in the analysed WHO GLASS bloodstream-infection surveillance data, 2020-2023. The number of countries and territories contributing relevant observations increased from 81 in 2020 to 95 in 2023. Because the set of contributing surveillance systems changed over time, temporal comparisons may partly reflect changes in surveillance composition.

Figure 2 - Volume of reported antimicrobial susceptibility testing
File: figures/02_reported_AST_volume_over_time.png

Question Addressed
How did the total volume of interpretable antimicrobial susceptibility testing represented in the analytical dataset change between 2020 and 2023?

Data and Calculation
For each year, InterpretableAST was summed across all included country-pathogen-antibiotic observations.
Approximate totals:

2020: 3.03 million interpretable AST results
2021: 3.56 million
2022: 3.59 million
2023: 3.48 million
What the Figure Shows
Reported testing volume increased substantially from 2020 to 2021, remained similar in 2022, and decreased slightly in 2023.
This differs from the country-participation pattern, which increased throughout the entire period.

Interpretation
The dataset changes over time in more than one way. The number of participating countries increased, but the total volume of AST results did not increase monotonically.
This demonstrates why surveillance composition must be considered when interpreting pooled resistance trends.

Main Caveats
InterpretableAST represents susceptibility-test results, not necessarily unique bacterial isolates or unique patients.
A single isolate may contribute results for multiple antibiotics.
Countries with high testing volumes contribute much more heavily to pooled AST-level summaries.
Portfolio / Manuscript Caption
Figure 2. Volume of interpretable antimicrobial susceptibility testing represented in the analysed WHO GLASS data, 2020-2023. Total reported interpretable AST results increased from approximately 3.03 million in 2020 to 3.56 million in 2021 and 3.59 million in 2022, before decreasing slightly to 3.48 million in 2023. Changes in surveillance volume did not mirror changes in country participation exactly.

Figure 3 - Pooled proportion of reported AST results classified as resistant
File: figures/03_pooled_resistant_AST_over_time.png

Question Addressed
What proportion of all reported interpretable AST results in the included dataset were classified as resistant in each year?

Data and Calculation
For each year:
pooled resistance (%) = total resistant AST results / total interpretable AST results × 100
Observed descriptive values:

2020: 18.4%
2021: 21.0%
2022: 20.5%
2023: 23.2%
What the Figure Shows
The pooled proportion of AST results classified as resistant was higher in 2023 than in 2020, but the pattern was not monotonic: it increased in 2021, decreased slightly in 2022, and increased again in 2023.

Interpretation
This is a descriptive AST-level surveillance summary. It is useful for showing the overall pattern in the collected surveillance records and for motivating the later inferential analysis.
It is not an estimate of the percentage of infections globally that are antimicrobial resistant.

Main Caveats
The pooled estimate combines:

different countries and territories,
different pathogens,
different antibiotics,
different surveillance volumes,
and changing surveillance participation.
Because AST counts are summed, observations from high-volume surveillance systems have greater influence.
The measure should therefore be described as the pooled proportion of reported AST results classified as resistant, not as global AMR prevalence.
Portfolio / Manuscript Caption
Figure 3. Pooled proportion of reported antimicrobial susceptibility test results classified as resistant, 2020-2023. Across the included WHO GLASS bloodstream-infection surveillance observations, the pooled proportion of AST results classified as resistant increased from 18.4% in 2020 to 23.2% in 2023, with a small decline between 2021 and 2022. The measure pools different pathogens, antibiotics, countries and surveillance volumes and should not be interpreted as a population-representative estimate of global AMR prevalence.

Figure 4 - Pooled resistant AST results by pathogen
File: figures/04_pooled_AST_resistance_by_pathogen.png

Question Addressed
How did pooled AST-level resistance summaries differ between the five included pathogens across 2020-2023?

Data and Calculation
Within each pathogen and year:
pooled resistance (%) = total resistant AST results / total interpretable AST results × 100
All antibiotics reported for that pathogen were combined.

What the Figure Shows
The five pathogens display substantially different descriptive resistance patterns.
Notable visual patterns include:

Acinetobacter spp. has the highest pooled resistant-AST proportion throughout the study period and remains broadly around one-half of reported AST results.
Klebsiella pneumoniae shows a sustained increase across the four years.
Escherichia coli remains at a lower level than K. pneumoniae and Acinetobacter spp., with a modest overall increase by 2023.
Salmonella spp. rises initially and then declines by 2023.
Streptococcus pneumoniae remains comparatively low and relatively stable.
Interpretation
The figure is useful for identifying broad pathogen-specific patterns and motivating pathogen-specific statistical analyses.
However, each pathogen-level line combines results from several antibiotics. A change in the mix of antibiotics tested can therefore change the pooled value even if resistance to individual antibiotics does not change.

Main Caveats
Antibiotic composition differs between pathogens.
Antibiotic-testing composition may also change over time within a pathogen.
The figure is weighted by AST volume.
These lines should not be interpreted as a single biological measure of "resistance" for each pathogen.
Portfolio / Manuscript Caption
Figure 4. Pooled proportion of reported AST results classified as resistant by pathogen, 2020-2023. Descriptive resistance patterns differed markedly between pathogens. Acinetobacter spp. showed the highest pooled values throughout the period, while Klebsiella pneumoniae displayed a sustained increase. Because results are pooled across antibiotics within each pathogen, temporal differences may partly reflect changes in antibiotic-testing composition.

Figure 5 - Reported resistance by pathogen and antibiotic, 2023
File: figures/05_pathogen_antibiotic_heatmap_2023.png

Question Addressed
How did reported resistance differ across specific pathogen-antibiotic combinations in 2023?

Data and Calculation
The dataset was restricted to 2023.
For each pathogen-antibiotic combination:
resistance (%) = total resistant AST results / total interpretable AST results × 100
Values were pooled across contributing countries and territories.
Grey cells represent pathogen-antibiotic combinations without observations in the analysed 2023 data.

What the Figure Shows
Resistance varies substantially both between pathogens and between antibiotics within the same pathogen.
The heatmap also demonstrates that the set of antibiotics represented differs across pathogens. This is important because it explains why a single pooled pathogen-level resistance percentage can conceal substantial antibiotic-specific heterogeneity.

Interpretation
This figure provides a more scientifically interpretable view of resistance than a single all-antibiotic pooled estimate because each coloured cell refers to a defined pathogen-antibiotic combination.
The figure should primarily be used to identify patterns for subsequent targeted analysis rather than to rank pathogens or antibiotics without considering surveillance coverage.

Main Caveats
Cells are weighted by the volume of reported AST results.
Countries contributing to one pathogen-antibiotic combination may differ from those contributing to another.
The number of countries and AST results underlying each cell can vary considerably.
The figure does not display uncertainty.
A visually high resistance percentage may be based on much less surveillance information than another cell.
Missing combinations should not be interpreted as zero resistance.
Portfolio / Manuscript Caption
Figure 5. Reported antimicrobial resistance by pathogen-antibiotic combination in 2023. Cell values represent the pooled proportion of interpretable AST results classified as resistant across contributing WHO GLASS surveillance systems. Grey cells indicate combinations without observations in the analysed 2023 data. Resistance patterns varied markedly between specific pathogen-antibiotic combinations, highlighting the limitations of interpreting resistance only at an all-antibiotic pathogen level.

Figure 6 - Temporal patterns in the two primary pathogen-antibiotic signals
File: figures/06_selected_temporal_trajectories.png

Question Addressed
How did reported resistance change from 2020 to 2023 for the pathogen-antibiotic combinations identified by the predefined primary inferential analysis, and how sensitive were the descriptive trajectories to surveillance weighting and country participation?

Data and Calculation
The figure includes the two pathogen-antibiotic combinations with BH-FDR < 0.05 in the predefined Script 04 primary analysis:

Escherichia coli - meropenem
Klebsiella pneumoniae - levofloxacin
For each combination and year, two descriptive summaries were calculated.
Pooled AST-weighted estimate
resistance (%) = total resistant AST results / total interpretable AST results × 100
Countries contributing more AST results therefore have greater influence.

Stable-country equal-weight estimate
Countries were first restricted to those reporting that specific pathogen-antibiotic combination in all four study years. Resistance percentages were then calculated separately for each country and year, followed by an unweighted mean across countries.
Each stable country therefore contributes equally regardless of AST volume.
The figure uses independent y-axis scales across the two panels because the absolute resistance levels differ substantially between the two pathogen-antibiotic combinations.

What the Figure Shows
For E. coli-meropenem, both the pooled AST-weighted trajectory and the stable-country equal-weight trajectory increase from 2020 to 2023.
The pooled resistance percentage increased from approximately 0.47% in 2020 to 1.77% in 2023. In the equal-weight analysis restricted to stable countries, the 2020-to-2023 increase was approximately 1.31 percentage points.
For K. pneumoniae-levofloxacin, the pooled resistance percentage increased from approximately 12.3% in 2020 to 17.5% in 2023. The equal-weight stable-country analysis also increased, with an estimated 2020-to-2023 change of approximately 7.39 percentage points.
The difference between pooled AST-weighted and equal-country estimates also demonstrates that absolute resistance percentages depend strongly on how surveillance systems are weighted.

Interpretation
The figure supports the directional robustness of the two primary signals.
For both pathogen-antibiotic combinations, the direction of change remains upward when:

all reported AST results are pooled according to surveillance volume, and
the analysis is restricted to consistently reporting countries and gives each country equal weight.
This reduces concern that the observed direction is solely a consequence of changing country participation or a small number of very high-volume surveillance systems.
The two lines should not be expected to overlap because they estimate different descriptive quantities.
Main Caveats
The pooled AST-weighted estimate gives more influence to countries contributing more AST results.
The stable-country equal-weight estimate answers a different question and is not a replacement for the primary model.
Country-level surveillance remains heterogeneous even after equal weighting.
Stable-country restriction improves temporal comparability but reduces the set of contributing surveillance systems.
AST results are not necessarily unique isolates or unique patients.
Independent y-axis scales improve within-combination readability but mean that vertical distances should not be compared directly between panels.
Selection of the two displayed combinations was based on the predefined Script 04 BH-FDR criterion, not on visual appearance or the later sensitivity analyses.
Portfolio / Manuscript Caption
Figure 6. Temporal patterns in the two primary pathogen-antibiotic signals, 2020-2023. Pooled AST-weighted resistance percentages are compared with equal-weight estimates restricted to countries reporting each pathogen-antibiotic combination in all four study years. Both Escherichia coli-meropenem and Klebsiella pneumoniae-levofloxacin showed increasing trajectories under both descriptive weighting approaches. Panels use independent y-axis scales. Percentages describe reported surveillance data and should not be interpreted as population-representative global prevalence.

Figure 7 - Sensitivity analyses of the two primary temporal signals
Current file: figures/07_primary_signal_sensitivity_comparison.png
Final report label: Figure 7

Question Addressed
Do the two primary temporal signals retain a similar direction and approximate magnitude under alternative reasonable modelling assumptions?

Data and Calculation
The figure compares four modelling approaches for the two combinations identified by the primary Script 04 analysis.

Primary robust model
Country fixed effects with country-clustered robust standard errors.

Beta-binomial model
Country random intercept with a beta-binomial distribution to allow extra-binomial variation.

Stable-country model
Country fixed-effects model with country-clustered robust standard errors, restricted to countries reporting that specific pathogen-antibiotic combination in all four study years.

Categorical-year model
Country fixed-effects model with country-clustered robust standard errors treating year as a categorical variable, allowing direct comparison of 2023 with 2020 without assuming a linear annual trend.
The primary robust, beta-binomial and stable-country models estimate an odds ratio per year. For the figure, these were converted to a common 2023-versus-2020 scale:
OR 2023 vs 2020 = (OR per year)^3
The corresponding confidence-interval limits were transformed in the same way.
The categorical-year model directly estimates 2023 versus 2020 and therefore required no transformation.

What the Figure Shows
All four modelling approaches produced odds ratios above 1 for both primary signals.
For E. coli-meropenem:

primary robust OR/year: approximately 1.13
beta-binomial OR/year: approximately 1.12
stable-country OR/year: approximately 1.13
categorical 2023-versus-2020 OR: approximately 1.55
For K. pneumoniae-levofloxacin:
primary robust OR/year: approximately 1.09
beta-binomial OR/year: approximately 1.11
stable-country OR/year: approximately 1.08
categorical 2023-versus-2020 OR: approximately 1.30
The estimates therefore remain directionally consistent across substantially different modelling assumptions.
Interpretation
This is the principal robustness figure for the inferential findings.
For E. coli-meropenem, the increasing association remained supported across the primary robust analysis, beta-binomial model, stable-country restriction and categorical-year analysis.
For K. pneumoniae-levofloxacin, the direction and approximate effect magnitude were also consistent across approaches. The stable-country analysis retained an odds ratio above 1 with a confidence interval above 1, although its BH-adjusted p-value was approximately 0.071.
The purpose of the figure is not to identify the method producing the smallest p-value. Instead, it asks whether the substantive conclusion changes when important modelling assumptions are altered.

Main Caveats
The models are not identical and do not estimate perfectly interchangeable quantities.
Logistic-model odds ratios are conditional on the model specification and should not be interpreted as risk ratios.
The primary, beta-binomial and stable-country effects shown in the figure are three-year transformations of a continuous-year model and therefore inherit the assumption of a constant log-odds slope across 2020-2023.
The categorical-year model avoids that linearity assumption but uses a different parameterisation.
Stable-country restriction reduces changing-participation bias but changes the analysed surveillance population.
Beta-binomial modelling addresses extra-binomial variation but does not remove all possible surveillance heterogeneity or selection bias.
Consistency across sensitivity analyses strengthens robustness but does not establish causality or population representativeness.
Portfolio / Manuscript Caption
Figure 7. Sensitivity analyses of the two primary temporal signals. Temporal associations for Escherichia coli-meropenem and Klebsiella pneumoniae-levofloxacin were compared across the primary country-fixed-effect model with country-clustered robust standard errors, beta-binomial mixed models, stable-country analyses and categorical-year models. Continuous-year estimates were transformed to a common 2023-versus-2020 odds-ratio scale. Both pathogen-antibiotic combinations retained an increasing direction across the alternative modelling approaches, although the strength of statistical evidence varied between specifications.

Supplementary Figure S1 - Primary temporal estimates across all pathogen-antibiotic combinations
Current file: figures/S1_primary_temporal_forest_plot.png
Final report label: Supplementary Figure S1

Question Addressed
What were the estimated annual temporal associations and 95% confidence intervals for all 37 pathogen-antibiotic combinations included in the predefined primary inferential analysis?

Data and Calculation
The figure displays the primary Script 04 estimate for each of the 37 eligible pathogen-antibiotic combinations.
For each combination, a grouped-binomial logistic regression was fitted with:

continuous year centred at 2020,
country fixed effects,
country-clustered robust standard errors.
Eligibility required:
observations in all four study years,
representation in both 2020 and 2023,
at least 20 contributing countries in every year.
The x-axis shows the odds ratio for reported resistance per additional year on a logarithmic scale.
The vertical dashed line represents OR = 1.
BH false-discovery-rate correction was applied across the 37 predefined pathogen-antibiotic models.
Filled points indicate combinations with BH-FDR < 0.05.
What the Figure Shows
Most estimates are relatively close to OR = 1 and/or have confidence intervals that include 1 after accounting for country clustering.
Two combinations met the predefined BH-FDR < 0.05 criterion:

Escherichia coli - meropenem
Klebsiella pneumoniae - levofloxacin
Several additional combinations had point estimates above or below 1, but uncertainty was greater once country-clustered robust inference and multiplicity correction were applied.
The Salmonella spp. estimates generally have wider confidence intervals than the better-supported E. coli and K. pneumoniae combinations, consistent with lower surveillance coverage for several Salmonella combinations.
Interpretation
This figure is important because it displays the entire predefined primary analysis family rather than only the two combinations that crossed the FDR threshold.
It therefore provides context for the selected primary signals and reduces the risk of presenting only favourable results.
The figure also shows why interpretation should consider effect magnitude, confidence intervals, surveillance coverage and sensitivity analyses rather than p-values alone.

Main Caveats
Odds ratios are model-based measures of association, not absolute resistance percentages.
Continuous year assumes a constant change in log odds per year over the four-year period.
The models remain based on reported surveillance data and are not population-representative estimates of global resistance.
Country fixed effects control for time-invariant differences between countries but do not address all time-varying surveillance differences.
Country-clustered robust standard errors address within-country dependence and variance misspecification but cannot model unobserved isolate-level dependence because isolate identifiers are unavailable.
The same isolate may contribute AST results for more than one antibiotic.
Statistical significance after FDR correction does not by itself establish clinical importance or causality.
Portfolio / Manuscript Caption
Supplementary Figure S1. Primary annual temporal associations across 37 pathogen-antibiotic combinations reported to WHO GLASS, 2020-2023. Points show odds ratios for reported resistance per additional year from country-fixed-effect grouped-binomial models with country-clustered robust standard errors; horizontal lines show 95% confidence intervals. The dashed vertical line indicates OR = 1. Filled points identify combinations with BH-FDR < 0.05 across the 37 predefined primary analyses. The figure presents the complete primary analysis family rather than only statistically selected results.

Final Figure Set
The current main-report figure structure is:

Figure 1: Country participation over time
Figure 2: Reported AST volume over time
Figure 3: Pooled proportion of reported AST results classified as resistant
Figure 4: Pooled resistant AST results by pathogen
Figure 5: Reported resistance by pathogen and antibiotic in 2023
Figure 6: Temporal patterns in the two primary pathogen-antibiotic signals
Figure 7: Sensitivity analyses of the two primary temporal signals
Supplementary:
Supplementary Figure S1: Primary temporal estimates across all 37 pathogen-antibiotic combinations
This ordering creates a progression from surveillance coverage, to broad descriptive patterns, to pathogen-antibiotic specificity, to primary inferential findings, and finally to robustness.
Figure Interpretation Rules
Throughout the project, figure titles and captions should describe the quantity actually calculated.
In particular:

pooled AST-level percentages should not be labelled as global AMR prevalence;
InterpretableAST should be described as a count of reported susceptibility-test results, not assumed to represent unique isolates or patients;
cross-antibiotic summaries should be identified as composition-dependent;
pathogen-antibiotic-specific estimates should be preferred when discussing specific resistance phenotypes;
surveillance coverage and changing country participation should accompany temporal interpretation;
descriptive figures should not be presented as causal evidence;
odds ratios should be distinguished from absolute percentage changes;
sensitivity analyses should be used to assess robustness to reasonable modelling choices rather than to select the smallest p-value;
BH-FDR thresholds should not be treated as the sole criterion for scientific importance;
changing surveillance composition, unequal AST volume and residual heterogeneity should remain explicit limitations of the project.