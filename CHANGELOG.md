# nf-core/differentialabundance: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.5.0

### `Added`

- [[#266](https://github.com/nf-core/differentialabundance/pull/266)] - Fix logging by specifying assays to log ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#259](https://github.com/nf-core/differentialabundance/pull/259)] - Bump gtf2featureannotation to fix GTF handling error ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#257](https://github.com/nf-core/differentialabundance/pull/257)] - Added maxquant profile to nextflow.config to make it available ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#254](https://github.com/nf-core/differentialabundance/pull/254)] - Some parameter changes, added qbic credits ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#250](https://github.com/nf-core/differentialabundance/pull/250)] - Template update for nf-core/tools v2.13.1 ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#244](https://github.com/nf-core/differentialabundance/pull/244)] - Add pipeline params for matrixfilter NA options ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#241](https://github.com/nf-core/differentialabundance/pull/241)] - Template update for nf-core/tools v2.13 ([@WackerO](https://github.com/WackerO), review by [@nvnieuwk](https://github.com/nvnieuwk))
- [[#228](https://github.com/nf-core/differentialabundance/pull/228)] - Template update for nf-core/tools v2.12 ([@nf-core-bot](https://github.com/nf-core-bot), review by [@pinin4fjords](https://github.com/pinin4fjords), [@WackerO](https://github.com/WackerO))
- [[#222](https://github.com/nf-core/differentialabundance/pull/222)] - Add rounding to all numeric report tables ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#219](https://github.com/nf-core/differentialabundance/pull/219)] - Template update for nf-core/tools v2.11.1 ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#197](https://github.com/nf-core/differentialabundance/pull/197)] - Add contributor info to report ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [[#267](https://github.com/nf-core/differentialabundance/pull/267)] - Whitespace fix, remove TODO, also update changelog for release release 1.5.0 ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#265](https://github.com/nf-core/differentialabundance/pull/265)] - GSEA- pngs and htmls in same place ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#257](https://github.com/nf-core/differentialabundance/pull/257)] - Fixed FILTER_DIFFTABLE module, updated PROTEUS module to better handle whitespace in prefix param, made docs clearer ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#254](https://github.com/nf-core/differentialabundance/pull/254)] - Made differential_file_suffix optional ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#240](https://github.com/nf-core/differentialabundance/pull/240)] - Publish GSEA reports ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#231](https://github.com/nf-core/differentialabundance/pull/231)] - Update GSEA module to fix butterfly plot bug ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#226](https://github.com/nf-core/differentialabundance/pull/226)] - Fix DESEQ2_NORM in modules.config ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#221](https://github.com/nf-core/differentialabundance/pull/221)] - Update shinyngs modules to address density plots issue ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [[#223](https://github.com/nf-core/differentialabundance/pull/223)] - tabulartogseacls fixes ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [[#213](https://github.com/nf-core/differentialabundance/pull/213)] - Fix volcano plot legend ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#210](https://github.com/nf-core/differentialabundance/pull/210)] - Include Affy timeout fix ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#208](https://github.com/nf-core/differentialabundance/pull/208)] - Fix resource issues and bump versions ([@pinin4fjords](https://github.com/pinin4fjords), review by [@sguizard](https://github.com/sguizard))

### `Changed`

- [[#256](https://github.com/nf-core/differentialabundance/pull/256)] - Release 1.5.0 ([@WackerO](https://github.com/WackerO), review by [@maxulysse](https://github.com/maxulysse), [@pinin4fjords](https://github.com/pinin4fjords))
- [[#264](https://github.com/nf-core/differentialabundance/pull/264)] - Change FILTER_DIFFTABLE to python because AWK does not correctly filter reliably ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#232](https://github.com/nf-core/differentialabundance/pull/232)] - Mention missing dots in volcano plot, change rounding, turn off rounding by default ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))

## v1.4.0 - 2023-11-27

### `Added`

- [[#203](https://github.com/nf-core/differentialabundance/pull/203)] - Transcript lengths for DESeq2 ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse))
- [[#199](https://github.com/nf-core/differentialabundance/pull/199)] - Add gprofiler2 module and local differential table filtering module ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#193](https://github.com/nf-core/differentialabundance/pull/193)] - Add DESeq2 text to report ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#192](https://github.com/nf-core/differentialabundance/pull/192)] - Add scree plot in report ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#189](https://github.com/nf-core/differentialabundance/pull/189)] - Add DE models to report ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#188](https://github.com/nf-core/differentialabundance/pull/188)] - Add option to cluster all features ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#198](https://github.com/nf-core/differentialabundance/pull/200)] - Document correct RNAseq matrix usage ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

### `Fixed`

- [[#201](https://github.com/nf-core/differentialabundance/issues/201)] - DESeq2/Limma update, fix incorrect column names issue ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#191](https://github.com/nf-core/differentialabundance/issues/191)] - Fix sample metadata table in the html report not paginating ([@davidecarlson](https://github.com/davidecarlson), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#190](https://github.com/nf-core/differentialabundance/pull/190)] - Fix GSEA indent in report ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))

### `Changed`

- [[#194](https://github.com/nf-core/differentialabundance/pull/194)] - Change report volcano colors ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#188](https://github.com/nf-core/differentialabundance/pull/188)] - Update min nextflow version ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))

## v1.3.1 - 2023-10-26

### `Fixed`

- [[#183](https://github.com/nf-core/differentialabundance/pull/183)] - Fix logging for dendrograms ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

### `Changed`

- [[#182](https://github.com/nf-core/differentialabundance/pull/182)] - Fixed Jon's employer ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

## v1.3.0 - 2023-10-24

### `Added`

- [[#124](https://github.com/nf-core/differentialabundance/pull/124)] - Template update for nf-core/tools v2.8 ([@pinin4fjords](https://github.com/pinin4fjords), review by [@jasmezz](https://github.com/jasmezz))
- [[#129](https://github.com/nf-core/differentialabundance/pull/129)] - Module updates to fit with recent registry changes ([@pinin4fjords](https://github.com/pinin4fjords), review by [@maxulysse](https://github.com/maxulysse), [@adamrtalbot](https://github.com/adamrtalbot))
- [[#130](https://github.com/nf-core/differentialabundance/pull/130)] - Document reasons for lack of differential expression ([@pinin4fjords](https://github.com/pinin4fjords), review by [@jfy133](https://github.com/jfy133))
- [[#131](https://github.com/nf-core/differentialabundance/pull/131)] - Improve gtf to table configurability ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#136](https://github.com/nf-core/differentialabundance/pull/136)] - Added support for non-Affymetrix arrays via automatic download of SOFT matrices in GEO ([@azedinez](https://github.com/azedinez), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#137](https://github.com/nf-core/differentialabundance/pull/137)] - Add `--sizefactors_from_controls` and `--gene_id_col` for DESeq2 module to modules.config ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#145](https://github.com/nf-core/differentialabundance/pull/145)] - Template update for nf-core/tools v2.9 ([@nf-core-bot](https://github.com/nf-core-bot), review by [@pinin4fjords](https://github.com/pinin4fjords), [@WackerO](https://github.com/WackerO))
- [[#147](https://github.com/nf-core/differentialabundance/pull/147)] - Add Maxquant analysis module ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#166](https://github.com/nf-core/differentialabundance/issues/166)] - Output a parameter-resolved R Markdown document, as well as rendered HTML ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#119](https://github.com/nf-core/differentialabundance/issues/119)] - Document sample sheet for Affy arrays ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#165](https://github.com/nf-core/differentialabundance/issues/165)] - Update subway map ([@pinin4fjords](https://github.com/pinin4fjords), review by [@FriederikeHanssen](https://github.com/FriederikeHanssen))
- [[#135](https://github.com/nf-core/differentialabundance/issues/135)] - workaround OPENBLAS using all cores problem ([@pinin4fjords](https://github.com/pinin4fjords), review by [@sateeshperi](https://github.com/sateeshperi))
- [[#176](https://github.com/nf-core/differentialabundance/pull/176)] - bump shinyngs ([@pinin4fjords](https://github.com/pinin4fjords), review by )

### `Fixed`

- [[#116](https://github.com/nf-core/differentialabundance/issues/116)] - Skip outlier detection with low replication ([@pinin4fjords](https://github.com/pinin4fjords), review by [@nvnieuwk](https://github.com/nvnieuwk))
- [[#122](https://github.com/nf-core/differentialabundance/pull/126)] - Add spaces to satisfy nf-core download for singularity ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#127](https://github.com/nf-core/differentialabundance/issues/127)] - [Bug] Can't pass samplesheet with -c file.config , or -params-file params.yml or directly with --input samplesheet.csv ([@ctuni](https://github.com/ctuni), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#138](https://github.com/nf-core/differentialabundance/issues/138)]- Fix bugs with --control_features and --sizefactors_from_controls ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#133](https://github.com/nf-core/differentialabundance/issues/133)] - Sample exclusion options fail due to contrast-wise normalisation ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#160](https://github.com/nf-core/differentialabundance/issues/160)]- Fix merge conflicts for Template update 2.10 by nf-core-bot ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#164](https://github.com/nf-core/differentialabundance/pull/164)] - Rlog + other small fixes ([@pinin4fjords](https://github.com/pinin4fjords), review by [@drpatelh](https://github.com/drpatelh))
- [[#174](https://github.com/nf-core/differentialabundance/pull/174)] - Fix metro map ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

### `Changed`

- [[#179](https://github.com/nf-core/differentialabundance/issues/179)] - Removed shiny app error message for proteus runs ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#159](https://github.com/nf-core/differentialabundance/issues/159)] - CUSTOM/MATRIXFILTER module update ([@WackerO](https://github.com/WackerO), review by [@suzannejin](https://github.com/suzannejin))
- [[#154](https://github.com/nf-core/differentialabundance/issues/154)] - RMARKDOWNNOTEBOOK env update ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#151](https://github.com/nf-core/differentialabundance/issues/151)] - Module update ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#147](https://github.com/nf-core/differentialabundance/pull/147)] - RMARKDOWNNOTEBOOK env update, SHINYNGS and CUSTOM update ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))

## v1.2.0 - 2023-04-19

### `Added`

- [[#97](https://github.com/nf-core/differentialabundance/issues/97)] - Allow for subsetting of samples for specific contrasts ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@danhalligan-hx](https://github.com/danhalligan-hx), review by [@WackerO](https://github.com/WackerO))
- [[#105](https://github.com/nf-core/differentialabundance/pull/105)] - Enabled multiple GMT/GMX files for GSEA ([@WackerO](https://github.com/WackerO), reported by [@grst](https://github.com/grst), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#108](https://github.com/nf-core/differentialabundance/issues/108)] - Add shiny app generation (starting feature set) ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#110](https://github.com/nf-core/differentialabundance/pull/110)] - Add shiny app outputs to tower.yml ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO), [@maxulysse](https://github.com/maxulysse))
- [[#149]](https://github.com/nf-core/differentialabundance/pull/149) - Update README.md - add ref to nf-core/rnaseq and Affymetrix ([@smoe](https://github.com/smoe), review by [@pinin4fjords](https://github.com/pinin4fjords))

### `Fixed`

- [[#95](https://github.com/nf-core/differentialabundance/issues/95)] - Pipeline doesn't check for gene sets file specification when GSEA is activated ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@danhalligan-hx](https://github.com/danhalligan-hx), review by [@FriederikeHanssen](https://github.com/FriederikeHanssen))
- [[#93](https://github.com/nf-core/differentialabundance/issues/93)] - Shouldn't be re-using the single exploratory palette across multiple informative variables ([@pinin4fjords](https://github.com/pinin4fjords), review by [@matthdsm](https://github.com/matthdsm))

## v1.1.1 - 2023-03-02

### `Fixed`

- [[#89](https://github.com/nf-core/differentialabundance/pull/89)] - Sanitise for differential ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

## v1.1.0 - 2023-02-22

### `Added`

- [[#63](https://github.com/nf-core/differentialabundance/issues/63)] - Add git CI matrix for different test profiles to run ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#60](https://github.com/nf-core/differentialabundance/pull/60)] - Add Affymetrix analysis modules, observation_name_col, make GTF optional, closing [[#47](https://github.com/nf-core/differentialabundance/issues/47)], [[#46](https://github.com/nf-core/differentialabundance/issues/46)], [[#66](https://github.com/nf-core/differentialabundance/issues/66)] ([@pinin4fjords](https://github.com/pinin4fjords), review by [@apeltzer](https://github.com/apeltzer))
- [[#70](https://github.com/nf-core/differentialabundance/pull/70)] - Array integration final bits: density plot improvement, output paths ([@pinin4fjords](https://github.com/pinin4fjords), review by [@matthdsm](https://github.com/matthdsm), [@WackerO](https://github.com/WackerO))

### `Fixed`

- [[#55](https://github.com/nf-core/differentialabundance/pull/55)] - Don't drop dimensions in report with single informative variable ([@pinin4fjords](https://github.com/pinin4fjords), review by [@mashehu](https://github.com/mashehu))
- [[#57](https://github.com/nf-core/differentialabundance/issues/57)] - Update module shinyngs/validatefomcomponents and bump other shinyngs modules to resolve [[#56](https://github.com/nf-core/differentialabundance/issues/56)] ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#65](https://github.com/nf-core/differentialabundance/pull/65)] - Update shinyngs modules to latest to resolve palette issues reported by [@Shellfishgene](https://github.com/Shellfishgene) ([@pinin4fjords](https://github.com/pinin4fjords), review by [@mashehu](https://github.com/mashehu))
- [[#67](https://github.com/nf-core/differentialabundance/issues/67)] - Error when contrast is blocking by multiple variables ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@clstacy](https://github.com/clstacy), review by [@jfy133](https://github.com/jfy133))
- [[#71](https://github.com/nf-core/differentialabundance/pull/71)] - Add palette options, bring full check.names fixes from upstream, closing [[#68](https://github.com/nf-core/differentialabundance/issues/68)] reported by [@Shellfishgene](https://github.com/Shellfishgene), [[#69](https://github.com/nf-core/differentialabundance/issues/69)] reported by [@clstacy](https://github.com/clstacy). ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@clstacy](https://github.com/clstacy), review by [@matthdsm](https://github.com/matthdsm), [@WackerO](https://github.com/WackerO))
- [[#73](https://github.com/nf-core/differentialabundance/pull/73)] - Make `--matrix`, `--input`, `--features` accept either CSV,TXT or TSV files as input ([@apeltzer](https://github.com/apeltzer), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#75](https://github.com/nf-core/differentialabundance/pull/75)] - Bump deseq2/differential for check.names fix ([@pinin4fjords](https://github.com/pinin4fjords), review by [@apeltzer](https://github.com/apeltzer))
- [[#76](https://github.com/nf-core/differentialabundance/pull/76)] - Fix up non-gtf operation using matrix as annotation source ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@apeltzer](https://github.com/apeltzer), review by [@apeltzer](https://github.com/apeltzer))
- [[#82](https://github.com/nf-core/differentialabundance/pull/82)] - Address v1.1.0 pre-release feedback (subway map fixes) ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@mashehu](https://github.com/mashehu))

## v1.0.1 - 2023-01-25

- [[#49](https://github.com/nf-core/differentialabundance/pull/49)] - Add citation fixes, missing logos, output detail, and trigger Zenodo ([@pinin4fjords](https://github.com/pinin4fjords), review by [@apeltzer](https://github.com/apeltzer), [@jfy133](https://github.com/jfy133))

## v1.0.0 - 2023-01-23

Initial release of nf-core/differentialabundance, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [[#1](https://github.com/nf-core/differentialabundance/issues/3)] - Set up initial modules, define and validate workflow inputs ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#4](https://github.com/nf-core/differentialabundance/issues/4)] - Add input checking ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#5](https://github.com/nf-core/differentialabundance/issues/5)] - Add differential analysis ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#6](https://github.com/nf-core/differentialabundance/issues/6)] - Add exploratory plotting ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#7](https://github.com/nf-core/differentialabundance/issues/7)] - Add differential plotting ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#8](https://github.com/nf-core/differentialabundance/issues/8)] - Establish outputs/ reports ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#12](https://github.com/nf-core/differentialabundance/issues/12)] - Handle spike sequences correctly ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#15](https://github.com/nf-core/differentialabundance/issues/15)] - Add filtering module ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#21](https://github.com/nf-core/differentialabundance/issues/21)] - Gene set analysis ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#22](https://github.com/nf-core/differentialabundance/issues/22)] - Complete docs- README etc ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#23](https://github.com/nf-core/differentialabundance/issues/23)] - Complete minimal report content ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#38](https://github.com/nf-core/differentialabundance/issues/38)] - Complete full-size test data setup ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))
- [[#41](https://github.com/nf-core/differentialabundance/issues/41)] - Update citations ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

### `Fixed`

### `Dependencies`

### `Deprecated`
