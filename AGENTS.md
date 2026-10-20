# AGENTS.md

Metagenomics-Toolkit (MGTK): a Nextflow DSL2 workflow for Illumina and ONT metagenomes
(QC, assembly, binning, annotation, plasmids, dereplication, co-occurrence, metabolic modelling, EMGB export).
The full docs are in `docs/` (published with mkdocs-material). Read `docs/developer_guidelines.md` before changing code.

## Layout

- `main.nf`: all entry points. `wFullPipeline` (per-sample + aggregation), `wAggregatePipeline`, `wExportPipeline`,
  and one entry per module (`wShortReadAssembly`, `wMagAttributes`, `wCooccurrence`, ...). `_wConfigurePipeline` forces internal settings.
- `modules/<module>/`: workflows and processes per module. Module names and versions live in the `modules` section of `nextflow.config`.
- `templates/*.sh`: shell bodies of processes, loaded with `template('x.sh')`.
  Two escaping styles exist: `shell:` processes use `!{var}` for Nextflow and `$var` for bash;
  `script:` processes use `${var}` for Nextflow and `\$var` for bash. Check which one a process uses before editing its template.
- `lib/Utils.groovy`, `lib/Output.groovy`: shared helpers (resources, time limits, output paths). Put reusable code here.
- `default/`: production configs. `example_params/`: demo/test configs only. `test_data/`: sample sheets.
- `scripts/test_*.sh`, `scripts/modules/*/test_*.sh`: tests; they are what CI runs (`.github/workflows/workflow_modules.yml`).
- `docs/`: user documentation. Code and YAML are pulled in with `---8<--- "path"` includes, so the referenced files must exist.

## Conventions

- Commits: `type(scope): message`, e.g. `fix(assembly): ...` (git-chglog builds the changelog). Branch off `dev`; `master` is for releases.
- Workflow names start with `w` (internal: `_w`), processes with `p`, camelCase everywhere, no magic numbers.
- Every process needs a resource label and a time limit (`Utils.setTimeLimit`), reads tool options from
  `params.steps.<module>.<process>.additionalParams`, and publishes `.command.*` logs but never `.command.run`.
- Outputs go to `SAMPLE/RUN_ID/MODULE/VERSION/TOOL/` (`AGGREGATED/...` for aggregation). Config and output layout are the public API:
  bump the module version in `nextflow.config` (semver) when they change. Version bumps are usually left to the release step; ask first.
- File and contig naming follow `docs/module_specification.md`.
- Do not put secrets in config YAMLs; the config is part of the pipeline output.
- New tools that need a database follow the download strategy in `docs/database.md` and `docs/developer_guidelines.md`.
- When adding or changing a tool, update its page in `docs/modules/` and add or adjust a test script and its CI step.

## Running tests

Nextflow is not installed system-wide. Install it with `make nextflow` (downloads the version from `VERSIONS.txt`
into the repo root; git ignores it). Java must be installed. Docker is required for the `standard` profile.

```bash
PATH=$PWD:$PATH bash scripts/modules/<module>/test_<name>.sh
```

- Tests download their inputs from GitHub/openstack, so they need network access.
- A run can take more than 10 minutes. Start it as a background Bash command; the foreground limit is 10 minutes.
- Scripts use `-resume` and finish with `make check`, which fails on any `IGNORE` in the trace.
  After changing a template, only that process reruns.
- Test runs leave `work/`, `output/`, `logs_<module>/`, `.nextflow/` and `.nextflow.log*` in the repo root.
  Don't commit them. Do not delete them on your own: ask the user first, and only clean up if they accept.
- Verify results, not just the exit code: open the files under `output/` and check them against what the change should do.
- Tests that need large databases (e.g. GTDB-Tk) are usually skipped locally and covered by CI. Say so when you skip one.

## Docs

```bash
make wiki_venv    # creates ./wiki_venv (git-ignored) and installs wiki_scripts/requirements.txt
make dev_wiki     # serves the docs locally with mkdocs
wiki_venv/bin/mkdocs build --config-file wiki_scripts/mkdocs.yml -d <site dir>   # static build, no Makefile target
```

Use the Makefile targets rather than creating a venv by hand. For a static build, write the site outside the repo (`-d`).
A clean build with no warnings is the check; also grep the built pages for leftover `---8<---` markers.
`make wiki_publish` pushes to GitHub Pages: never run it without being asked.

## Gotchas

- Templates are often copied from another tool. Check a copied template for leftover tool names, labels and paths.
- Use `find … -exec grep … {} +` rather than a bare glob with `grep`: an unmatched glob makes grep fail or read stdin.
- A tool's results only reach downstream steps if its output channel is mixed into the channel those steps read.
  When adding a tool, check every place the sibling tools' channels are used.
