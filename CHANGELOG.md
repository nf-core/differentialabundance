# nf-core/differentialabundance: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.1 - 2023-03-02

### `Fixed`

- [[#89](https://github.com/nf-core/differentialabundance/pull/89)] - Sanitise for differential ([@pinin4fjords](https://github.com/pinin4fjords), review by [@WackerO](https://github.com/WackerO))

## v1.1.0 - 2023-02-22

### `Added`

- [[#63](https://github.com/nf-core/differentialabundance/issues/63)] - Add git CI matrix for different test profiles to run ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#60](https://github.com/nf-core/differentialabundance/pull/60)] - Add Affymetrix analysis modules, observation_name_col, make GTF optional, closing [[#47](https://github.com/nf-core/differentialabundance/issues/47)], [[#46](https://github.com/nf-core/differentialabundance/issues/46)], [[#66](https://github.com/nf-core/differentialabundance/issues/66)] ([@pinin4fjords](https://github.com/pinin4fjords), review by [@apeltzer](https://github.com/apeltzer))
- [[#70](https://github.com/nf-core/differentialabundance/pull/70)] - Array integration final bits: density plot improvement, output paths ([@pinin4fjords](https://github.com/pinin4fjords), review by [@matthdsm](https://github.com/matthdsm), [@WackerO](https://github.com/WackerO))

### `Fixed`

- [[#55](https://github.com/nf-core/differentialabundance/pull/55)] - Don't drop dimensions in report with single informative variable ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@mashehu](https://github.com/mashehu))
- [[#57](https://github.com/nf-core/differentialabundance/issues/57)] - Update module shinyngs/validatefomcomponents and bump other shinyngs modules to resolve [[#56](https://github.com/nf-core/differentialabundance/issues/56)] ([@WackerO](https://github.com/WackerO), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#65](https://github.com/nf-core/differentialabundance/pull/65)] - Update shinyngs modules to latest to resolve palette issues reported by [@Shellfishgene](https://github.com/Shellfishgene) ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@mashehu](https://github.com/mashehu))
- [[#67](https://github.com/nf-core/differentialabundance/issues/67)] - Error when contrast is blocking by multiple variables ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@clstacy](https://github.com/clstacy), review by [@jfy133](https://github.com/jfy133))
- [[#71](https://github.com/nf-core/differentialabundance/pull/71)] - Add palette options, bring full check.names fixes from upstream, closing [[#68](https://github.com/nf-core/differentialabundance/issues/68)] reported by [@Shellfishgene](https://github.com/Shellfishgene), [[#69](https://github.com/nf-core/differentialabundance/issues/69)] reported by [@clstacy](https://github.com/clstacy). ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@clstacy](https://github.com/clstacy), review by [@matthdsm](https://github.com/matthdsm), [@WackerO](https://github.com/WackerO))
- [[#73](https://github.com/nf-core/differentialabundance/pull/73)] - Make `--matrix`, `--input`, `--features` accept either CSV,TXT or TSV files as input ([@apeltzer](https://github.com/apeltzer), review by [@pinin4fjords](https://github.com/pinin4fjords))
- [[#75](https://github.com/nf-core/differentialabundance/pull/75)] - Bump deseq2/differential for check.names fix ([@pinin4fjords](https://github.com/pinin4fjords), reviewed by [@apeltzer](https://github.com/apeltzer))
- [[#76](https://github.com/nf-core/differentialabundance/pull/76)] - Fix up non-gtf operation using matrix as annotation source ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@apeltzer](https://github.com/apeltzer), review by [@apeltzer](https://github.com/apeltzer))
- [[#82](https://github.com/nf-core/differentialabundance/pull/82)] - Address v1.1.0 pre-release feedback (subway map fixes) ([@pinin4fjords](https://github.com/pinin4fjords), reported by [@mashehu](https://github.com/mashehu))

## v1.0.1 - 2023-01-25

- [[#49](https://github.com/nf-core/differentialabundance/pull/49) - Add citation fixes, missing logos, output detail, and trigger Zenodo ([@pinin4fjords](https://github.com/pinin4fjords), review by [@apeltzer](https://github.com/apeltzer), [@jfy133](https://github.com/jfy133))

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
