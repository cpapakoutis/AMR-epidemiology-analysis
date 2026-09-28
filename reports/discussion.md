Discussion
Overview
This analysis examined temporal and pathogen-antibiotic-specific antimicrobial resistance patterns reported through WHO GLASS bloodstream-infection surveillance data between 2020 and 2023.
The central finding was that broad pooled resistance summaries and more carefully adjusted temporal models did not tell the same story.
At the descriptive level, the pooled proportion of reported AST results classified as resistant increased from 18.4% in 2020 to 23.2% in 2023. However, after accounting for country and pathogen-antibiotic composition, the overall adjusted temporal association was close to the null.
This contrast highlights a central challenge in interpreting surveillance data: an apparent change in a pooled summary can reflect not only changes in resistance, but also changes in which countries report, how much testing they contribute, and which pathogen-antibiotic combinations are represented.
The analysis therefore supports a cautious interpretation of broad AMR summaries and places greater emphasis on specific pathogen-antibiotic combinations and on robustness across alternative analytical approaches.
Interpreting the Pooled Increase
The crude pooled AST-level resistance percentage increased over the study period, but several features of the dataset make this quantity strongly composition-dependent.
First, the number of participating countries increased from 81 in 2020 to 95 in 2023. The surveillance population was therefore not identical across years.
Second, AST volume did not change in parallel with country participation. Total interpretable AST results increased substantially between 2020 and 2021, remained similar in 2022, and declined slightly in 2023. This means that both the number of participating countries and the relative contribution of individual surveillance systems changed over time.
Third, pooled summaries combine different pathogens and antibiotics. Because resistance prevalence differs substantially between pathogen-antibiotic combinations, a change in testing composition can alter the pooled percentage even if the underlying resistance level within specific combinations does not change.
These features explain why the crude increase from 18.4% to 23.2% should not be interpreted as evidence that "global AMR increased by 4.8 percentage points." The adjusted overall model estimated only a very small common temporal association, with an annual odds ratio close to 1 and no strong statistical evidence for a uniform increase across adequately represented pathogen-antibiotic combinations.
The difference between crude and adjusted results is therefore itself an important finding. It demonstrates that broad surveillance summaries can be highly sensitive to the composition and weighting of the contributing data.
Pathogen-Level Patterns
The descriptive pathogen-level analysis showed substantial heterogeneity.
Acinetobacter spp. had the highest pooled resistant-AST proportion across the study period, while Klebsiella pneumoniae showed a sustained descriptive increase. Escherichia coli remained lower overall, Salmonella spp. showed a non-monotonic pattern, and Streptococcus pneumoniae remained comparatively low.
However, these pathogen-level summaries pooled across multiple antibiotics and were therefore also composition-dependent.
The pathogen-specific adjusted models weakened considerably after controlling for antibiotic identity, country, and multiple testing. None of the five pathogen-level temporal associations remained below BH-FDR 0.05.
This supports the interpretation that temporal resistance patterns were not uniform across entire pathogens and were better understood at the pathogen-antibiotic level.
Primary Pathogen-Antibiotic Findings
Two pathogen-antibiotic combinations met the predefined BH-FDR < 0.05 criterion in the primary analysis:
- Escherichia coli - meropenem
- Klebsiella pneumoniae - levofloxacin
  These findings were not interpreted solely on the basis of statistical significance. Their effect sizes, absolute resistance patterns, coverage, model diagnostics, and behaviour across sensitivity analyses were considered together.
Escherichia coli - Meropenem
E. coli-meropenem showed the most consistently supported temporal increase in the project.
The primary model estimated an annual odds ratio of approximately 1.13, corresponding to an implied 2023-versus-2020 odds ratio of approximately 1.43 under the continuous-year model.
The pooled AST-weighted resistance percentage increased from approximately 0.47% in 2020 to 1.77% in 2023.
This finding is important to interpret in both relative and absolute terms. The relative increase was substantial, but the absolute resistance percentage remained comparatively low.
The increasing direction persisted across:
- the primary country-fixed-effect model;
- beta-binomial modelling;
- stable-country restriction;
- categorical-year modelling;
- equal-country weighting;
- stable-country equal weighting.
  This consistency suggests that the observed increase was not solely an artefact of changing country participation, AST-volume weighting, or the ordinary binomial variance assumption.
Klebsiella pneumoniae - Levofloxacin
K. pneumoniae-levofloxacin also showed a consistent increasing pattern.
The primary model estimated an annual odds ratio of approximately 1.09, while pooled AST-weighted resistance increased from approximately 12.3% in 2020 to 17.5% in 2023.
Compared with E. coli-meropenem, the absolute increase was larger.
The direction and approximate magnitude remained similar across beta-binomial, categorical-year, stable-country, and equal-weight analyses.
However, the stable-country BH-adjusted p-value moved above 0.05, despite the estimated odds ratio remaining above 1 and its 95% confidence interval remaining above 1.
This illustrates why significance thresholds should not be treated as the sole criterion for interpretation. The signal remained directionally stable, but the multiplicity-adjusted evidence was more sensitive to the reduced stable-country dataset.
What the Mixed-Effects Models Revealed
The ordinary binomial random-intercept GLMMs produced one of the most important methodological lessons in the project.
The GLMM and primary fixed-effect models produced very similar point estimates and the same direction for all 37 pathogen-antibiotic combinations.
However, the ordinary binomial GLMMs produced much smaller standard errors and far more statistically significant results.
Thirty of the 37 GLMM models had BH-FDR < 0.05, compared with only two in the primary robust analysis.
This difference occurred despite substantial residual overdispersion.
The overall adjusted GLMM had a Pearson dispersion of approximately 67.7, and many individual pathogen-antibiotic models also showed dispersion well above 1.
The key implication is that adding a country random intercept did not fully account for the extra-binomial heterogeneity in the surveillance data.
The GLMMs therefore demonstrated that model complexity alone does not guarantee appropriate uncertainty estimation. In this setting, the variance assumptions mattered at least as much as the point-estimate specification.
This is why the ordinary binomial GLMM p-values were treated as diagnostic and robustness information rather than as the primary inferential evidence.
Beta-Binomial Sensitivity Analysis
The beta-binomial analysis directly addressed the extra-binomial variation observed in the ordinary GLMMs.
All 37 beta-binomial models converged successfully with positive-definite Hessians.
The number of BH-FDR-positive combinations decreased from 30 under the ordinary binomial GLMM to 11 under the beta-binomial specification.
This supports the interpretation that the ordinary binomial models had been too confident in the presence of excess heterogeneity.
Importantly, both primary signals remained positive and of similar magnitude under the beta-binomial specification.
This increased confidence that the two primary findings were not merely products of under-modelled binomial variance.
Stable-Country Sensitivity Analysis
Changing country participation was one of the most important potential sources of bias in the project.
The stable-country analysis restricted each pathogen-antibiotic model to countries reporting that specific combination in all four years.
This substantially reduced the possibility that temporal differences were driven simply by countries entering or leaving the dataset.
Only one of the 37 stable-country models remained below BH-FDR 0.05.
This indicates that many apparent associations were sensitive to the composition of the reporting countries.
However, E. coli-meropenem remained strongly supported, and K. pneumoniae-levofloxacin retained a similar positive effect estimate.
The stable-country results therefore strengthened confidence in the direction of the two primary signals while also demonstrating that many other apparent temporal associations were less robust.
Categorical-Year Sensitivity Analysis
The primary continuous-year models assumed a constant linear change in log odds between 2020 and 2023.
Because the study covered only four years, this assumption could be restrictive.
The categorical-year analysis removed that assumption and directly compared 2023 with 2020.
Both primary signals remained positive and statistically supported under this alternative parameterisation.
This suggests that the main findings were not dependent on imposing a straight-line temporal trend through 2021 and 2022.
Equal-Country Weighting
The pooled AST-level summaries give more influence to countries contributing larger numbers of tests.
This may be appropriate if the estimand is the proportion resistant across all reported AST results, but it can also allow a small number of high-volume surveillance systems to dominate descriptive trends.
The equal-country-weighted analysis addressed a different question: what happened to the average country-level resistance percentage among contributors?
Both primary signals remained positive under equal-country weighting and under equal weighting restricted to stable reporters.
This further supports the conclusion that the direction of these two signals was not driven solely by surveillance volume.
Why Many Associations Were Model-Sensitive
Across all 37 pathogen-antibiotic combinations, only 21 maintained the same direction across all major sensitivity approaches.
The number of BH-FDR-positive results varied substantially by model:
- primary robust analysis: 2
- ordinary binomial GLMM: 30
- beta-binomial analysis: 11
- stable-country analysis: 1
- categorical-year analysis: 3
  This variability should not be treated as a failure of the analysis.
  Instead, it demonstrates that AMR surveillance data are structurally complex and that temporal conclusions can depend heavily on assumptions about:
- country participation;
- surveillance volume;
- within-country dependence;
- overdispersion;
- year parameterisation;
- and weighting.
  The sensitivity analysis therefore provides an important safeguard against over-interpreting any single model.
Strengths
Several features strengthen the analysis.
First, the workflow was explicitly designed to distinguish descriptive surveillance summaries from inferential conclusions.
Second, the project evaluated surveillance participation and AST volume before interpreting temporal change.
Third, pathogen-antibiotic eligibility criteria were defined before interpreting the final model results, reducing the risk of highlighting extremely sparse combinations.
Fourth, the primary analysis accounted for repeated observations within countries using country fixed effects and country-clustered robust standard errors.
Fifth, multiple testing was addressed using BH-FDR correction across the predefined 37-model family.
Sixth, sensitivity analyses challenged the major assumptions identified during the analysis, including:
- changing country participation;
- linearity of the year effect;
- AST-volume weighting;
- country heterogeneity;
- extra-binomial variation.
  Finally, the two main findings were interpreted using both relative measures and absolute resistance changes rather than p-values alone.
Limitations
The project also has important limitations.
Surveillance Representativeness
WHO GLASS data are surveillance data rather than a probability sample of the global population.
Countries differ in:
- surveillance coverage;
- laboratory capacity;
- healthcare access;
- diagnostic practices;
- testing intensity;
- reporting completeness;
- case mix.
  The analysis therefore cannot estimate population-representative global AMR prevalence.
Changing Participation
The set of reporting countries changed between 2020 and 2023.
Although country adjustment and stable-country sensitivity analyses reduce this concern, they cannot remove all time-varying differences in surveillance systems.
Aggregate AST Counts
The available data contain aggregate AST counts rather than patient- or isolate-level records.
InterpretableAST should therefore not be treated as a count of unique patients or isolates.
The same isolate may contribute results for multiple antibiotics, creating dependence across antibiotic-specific observations that cannot be modelled directly without isolate identifiers.
Residual Heterogeneity
Substantial overdispersion was present in several model specifications.
The beta-binomial sensitivity analysis addressed part of this heterogeneity, but no model can account for all unmeasured surveillance differences.
Short Time Series
The analysis covers only four calendar years.
This limits the ability to distinguish sustained long-term trends from short-term fluctuations.
Multiple Comparisons
Thirty-seven pathogen-antibiotic combinations were tested in the primary family.
BH-FDR correction reduces the expected proportion of false discoveries among selected results, but multiplicity-adjusted thresholds do not transform findings into definitive biological conclusions.
Observational Interpretation
The analysis is observational.
Potential drivers such as antibiotic consumption, infection-control practices, healthcare disruption, COVID-19-related changes, laboratory capacity, and One Health factors were not directly linked to the analytical dataset.
They may provide context but cannot be inferred as causes of the observed patterns.
Implications
The project demonstrates that interpretation of AMR surveillance data benefits from moving beyond a single pooled resistance percentage.
A crude pooled increase can coexist with a near-null adjusted common temporal effect when the reporting population, testing volume, and pathogen-antibiotic composition change over time.
For surveillance analysis, this supports several practical principles:
- report surveillance coverage alongside resistance estimates;
- distinguish AST-weighted summaries from country-weighted summaries;
- avoid treating pooled cross-antibiotic resistance as a single biological quantity;
- examine pathogen-antibiotic combinations directly;
- account for repeated country observations;
- assess overdispersion;
- correct for multiple testing;
- use sensitivity analyses to evaluate whether conclusions depend on modelling choices.
Conclusion
Reported antimicrobial resistance patterns in WHO GLASS bloodstream-infection surveillance data changed between 2020 and 2023, but the changes were heterogeneous and strongly dependent on pathogen, antibiotic, surveillance composition, and analytical assumptions.
The crude pooled proportion of resistant AST results increased over the study period, but this did not translate into a clear common adjusted temporal increase across adequately represented pathogen-antibiotic combinations.
Among 37 predefined pathogen-antibiotic analyses, Escherichia coli-meropenem and Klebsiella pneumoniae-levofloxacin were the two combinations meeting the primary BH-FDR criterion.
E. coli-meropenem showed the most consistently supported increase across all major sensitivity analyses.
K. pneumoniae-levofloxacin also showed a stable increasing direction and similar effect magnitude across analyses, although the strength of multiplicity-adjusted evidence was somewhat more sensitive to the stable-country restriction.
Overall, the analysis illustrates both the value and the limitations of international AMR surveillance data. Robust interpretation requires attention not only to resistance percentages, but also to surveillance participation, weighting, dependence, heterogeneity, uncertainty, and sensitivity to analytical choices.