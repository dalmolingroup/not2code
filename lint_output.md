## `nf-core pipelines lint` overall result: Failed :x:

Posted for pipeline commit 70ad359

```diff
+| ✅ 148 tests passed       |+
#| ❔   3 tests were ignored |#
!| ❗  19 tests had warnings |!
-| ❌  25 tests failed       |-
```

<details>

### :x: Test failures:

* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `.github/workflows/nf-test.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `.github/actions/get-shards/action.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `.github/actions/nf-test/action.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `nf-test.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `tests/default.nf.test`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_cpus`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_memory`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_time`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Lines for loading custom profiles not found. File should contain:
```groovy
// Load nf-core custom profiles from different institutions
includeConfig params.custom_config_base && (!System.getenv('NXF_OFFLINE') || !params.custom_config_base.startsWith('http')) ? "${params.custom_config_base}/nfcore_custom.config" : "/dev/null"
```
* [nf_test_content](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nf_test_content) - 'tests/nextflow.config' does not exist
* [nf_test_content](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nf_test_content) - 'nf-test.config' does not exist
* [actions_awsfulltest](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_awsfulltest) - `.github/workflows/awsfulltest.yml` is not triggered correctly
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `tpm_threshold` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `pfam_db` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `multiqc_replace_names` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `multiqc_sample_names` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `mstrg_prep_script` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `plek` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `lncselect` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `filter_by_tpm` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `reference_genome` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `reference_gtf` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `validationSkipDuplicateCheck` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `validationS3PathCheck` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Param `monochromeLogs` from `nextflow config` not found in nextflow_schema.json

### :heavy_exclamation_mark: Test warnings:

* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `conf/igenomes_ignored.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found: `ro-crate-metadata.json`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - nf-validation has been detected in the pipeline. Please migrate to nf-schema: https://nextflow-io.github.io/nf-schema/latest/migration_guide/
* [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README did not have a Nextflow minimum version badge.
* [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README did not have an nf-core template version badge.
* [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README contains the placeholder `zenodo.XXXXXXX`. This should be replaced with the zenodo doi (after the first release).
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Remove this line if you don't need a FASTA file_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `awsfulltest.yml`: _You can customise AWS full pipeline tests as required_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `ci.yml`: _You can customise CI pipeline run tests as required_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `methods_description_template.yml`: _#Update the HTML below to your preferred methods description, e.g. add publication citation for this pipeline_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `base.config`: _Check the defaults for all processes_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `base.config`: _Customise requirements for specific processes._
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `test_full.config`: _Specify the paths to your full test data ( on nf-core/test-datasets or directly in repositories, e.g. SRA)_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `test_full.config`: _Give any required params for the test so that command line flags are not needed_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `output.md`: _Write this documentation describing your workflow's output_
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `usage.md`: _Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website._
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Optionally add in-text citation tools to this list._
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Optionally add bibliographic entries to this list._
* [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!_

### :grey_question: Tests ignored:

* [files_unchanged](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_unchanged) - Required pipeline config not found - {'manifest.contributors'}
* [actions_nf_test](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_nf_test) - '.github/workflows/nf-test.yml' not found
* [rocrate_readme_sync](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/rocrate_readme_sync) - `ro-crate-metadata.json` not found

### :white_check_mark: Tests passed:

* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.gitattributes`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.gitignore`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.nf-core.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.prettierignore`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.prettierrc.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `CHANGELOG.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `CITATIONS.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `CODE_OF_CONDUCT.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `LICENSE` or `LICENSE.md` or `LICENCE` or `LICENCE.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `nextflow_schema.json`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `nextflow.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `README.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/.dockstore.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/CONTRIBUTING.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/bug_report.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/config.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/feature_request.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/PULL_REQUEST_TEMPLATE.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/branch.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting_comment.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.html`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.txt`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/sendmail_template.txt`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/nf-core-nottocode_logo_light.png`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/modules.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/test.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/test_full.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/images/nf-core-nottocode_logo_light.png`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/images/nf-core-nottocode_logo_dark.png`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/output.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/usage.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `main.nf`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/multiqc_config.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/base.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/igenomes.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/awstest.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/awsfulltest.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `modules.json`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/bug_report.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/feature_request.md`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/workflows/push_dockerhub.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.markdownlint.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.nf-core.yaml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.yamllint.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `bin/markdown_to_html.r`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `conf/aws.config`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `docs/images/nf-core-nottocode_logo.png`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Checks.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Completion.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/NfcoreTemplate.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Utils.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Workflow.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/WorkflowMain.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/WorkflowNottocode.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `parameters.settings.json`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `pipeline_template.yml`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `Singularity`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/nfcore_external_java_deps.jar`
* [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.travis.yml`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Found nf-validation plugin
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.name`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.nextflowVersion`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.description`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.version`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.homePage`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.cpus`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.memory`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.time`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.outdir`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.input`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.mainScript`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.nf_required_version`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.container`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.singleEnd`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.igenomesIgnore`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.name`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.enable_conda`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``timeline.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``report.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``trace.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``dag.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``manifest.name`` began with ``nf-core/``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ``manifest.homePage`` began with https://github.com/nf-core/
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``dag.file`` ended with ``.html``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ``manifest.nextflowVersion`` started with >= or !>=
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config ``manifest.version`` ends in ``dev``: ``1.0dev``
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `params.custom_config_version` is set to `master`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `params.custom_config_base` is set to `https://raw.githubusercontent.com/nf-core/configs/master`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - nextflow.config contains configuration profile `test`
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_version= master
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_base= https://raw.githubusercontent.com/nf-core/configs/master
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_cpus= 16
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_memory= 128.GB
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_time= 240.h
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.publish_dir_mode= copy
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_multiqc_email_size= 25.MB
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.validate_params= true
* [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.pipelines_testdata_base_path= https://raw.githubusercontent.com/nf-core/test-datasets/
* [actions_awstest](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_awstest) - '.github/workflows/awstest.yml' is triggered correctly
* [actions_awsfulltest](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_awsfulltest) - `.github/workflows/awsfulltest.yml` does not use `-profile test`
* [pipeline_if_empty_null](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_if_empty_null) - No `ifEmpty(null)` strings found
* [plugin_includes](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/plugin_includes) - No wrong validation plugin imports have been found
* [pipeline_name_conventions](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_name_conventions) - Name adheres to nf-core convention
* [template_strings](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/template_strings) - Did not find any Jinja template strings (0 files)
* [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Schema lint passed
* [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Schema title + description lint passed
* [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Input mimetype lint passed: 'text/csv'
* [system_exit](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/system_exit) - No `System.exit` calls found
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: awsfulltest.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: awstest.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: branch.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: ci.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: clean-up.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: download_pipeline.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: fix-linting.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting_comment.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: release-announcements.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: release-announcments.yml
* [merge_markers](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/merge_markers) - No merge markers found in pipeline files
* [modules_json](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_json) - Only installed modules found in `modules.json`
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` found and not ignored.
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_section_order`
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `export_plots`
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_comment`
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` follows the ordering scheme of the minimally required plugins.
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains a matching 'report_comment'.
* [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains 'export_plots: true'.
* [modules_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_structure) - modules directory structure is correct 'modules/nf-core/TOOL/SUBTOOL'
* [local_component_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/local_component_structure) - local subworkflows directory structure is correct 'subworkflows/local/TOOL/SUBTOOL'
* [base_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/base_config) - `conf/base.config` found and not ignored.
* [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `conf/modules.config` found and not ignored.
* [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `FASTQC` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `MULTIQC` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `HMMER_HMMPRESS` found in `conf/modules.config` and Nextflow scripts.
* [nfcore_yml](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nfcore_yml) - Repository type in `.nf-core.yml` is valid: `pipeline`
* [nfcore_yml](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nfcore_yml) - nf-core version in `.nf-core.yml` is set to the latest version: `2.11.dev0`

### Run details

* nf-core/tools version 3.5.2
* Run at `2026-04-16 00:10:59`

</details>
