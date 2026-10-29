<#
  GCZ-RL - reproduce the paper with one command (Windows PowerShell).

  powershell -ExecutionPolicy Bypass -File .\RUN_ALL.ps1                 # Part 1 + Part 2 + checks + paper numbers
  powershell -ExecutionPolicy Bypass -File .\RUN_ALL.ps1 -Part 2         # primary LOSO only
  powershell -ExecutionPolicy Bypass -File .\RUN_ALL.ps1 -Fresh          # also recompute ICP geometry and SLD
  powershell -ExecutionPolicy Bypass -File .\RUN_ALL.ps1 -Slow           # also Table 3 (about 25 min), the original notebook, and the ICP seed analysis (about 100 min)
  powershell -ExecutionPolicy Bypass -File .\RUN_ALL.ps1 -Pdf            # also compile the manuscript and the supplement (needs latexmk, e.g. MiKTeX)

  Afterwards, compare every output with the committed version:
  .venv\Scripts\python.exe scripts\compare_outputs.py

  Raw data: fetched automatically into dataset\raw unless GCZRL_RAW_ROOT points to an existing copy.
  Python: .venv in this folder (created if missing) unless PY points to another interpreter.
#>
param([ValidateSet('1', '2', 'all')][string]$Part = 'all', [switch]$Fresh, [switch]$Slow, [switch]$Pdf)

$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Root
$Py = if ($env:PY) { $env:PY } else { Join-Path $Root '.venv\Scripts\python.exe' }
$env:PYTHONIOENCODING = 'utf-8'; $env:MPLBACKEND = 'Agg'; $env:PYTHONHASHSEED = '0'
New-Item -ItemType Directory -Force (Join-Path $Root 'logs') | Out-Null
$Log = Join-Path $Root ("logs\run_{0}.log" -f (Get-Date -Format 'yyyyMMdd_HHmmss'))

function Say([string]$t) { $l = "[{0}] {1}" -f (Get-Date -Format 'HH:mm:ss'), $t; Write-Host $l; Add-Content $Log $l }

function Run([string]$Name, [string[]]$PyArgs) {
    Say "START $Name"
    $prev = $ErrorActionPreference; $ErrorActionPreference = 'Continue'
    & $Py @PyArgs 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $Log -Append | Where-Object { $_ -match '^\[|PASS|FAIL|MISMATCH|verified|reproduced|Error|Traceback' }
    $code = $LASTEXITCODE; $ErrorActionPreference = $prev
    if ($code -ne 0) { Say "FAILED $Name (exit $code) - see $Log"; exit $code }
    Say "OK $Name"
}

if (-not (Test-Path $Py)) {
    Say 'Creating Python 3.11 environment (.venv) ...'
    if (Get-Command uv -ErrorAction SilentlyContinue) {
        & uv venv --python 3.11 .venv; & uv pip install --python $Py -r requirements.txt
    } else {
        & py -3.11 -m venv .venv; & $Py -m pip install -r requirements.txt
    }
}

$RawRoot = if ($env:GCZRL_RAW_ROOT) { $env:GCZRL_RAW_ROOT } else { Join-Path $Root 'dataset\raw' }
if (-not (Test-Path (Join-Path $RawRoot 'VerSe_Ligaments_GT'))) {
    Run 'Fetch and verify raw data' @('scripts\fetch_data.py', '--out', $RawRoot)
}
if ($Fresh -and -not (Test-Path (Join-Path $Root 'code\third_party\SpineLib'))) {
    Run 'Fetch SpineLib (SLD)' @('scripts\fetch_spinelib.py')
}

$f = @(); if ($Fresh) { $f = @('--fresh') }
if ($Part -in @('1', 'all')) {
    Run 'Part 1 (development, 57 vertebrae)' (@('part1_dissertation_57\run_part1.py') + $f)
    if ($Slow) {
        Run 'Original notebook (matching diagnostic, alignment)' @('part1_dissertation_57\original_notebook\run_original_notebook.py')
        Run 'Original notebook tables' @('part1_dissertation_57\original_notebook\extract_tables.py')
        Run 'Table 3 regressor comparison' @('paper\table3_regressor_comparison.py')
        Run 'ICP seed sensitivity (10 runs, about 100 min)' @('paper\icp_seed_sensitivity.py', 'run', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        Run 'ICP seed sensitivity summary' @('paper\icp_seed_sensitivity.py', 'summary')
    }
}
if ($Part -in @('2', 'all')) {
    Run 'Part 2 (primary LOSO, 14 vertebrae)' (@('part2_patient_loso_14\run_part2.py') + $f)
    Run 'Part 2 leakage test' @('part2_patient_loso_14\run_leakage_test.py')
}
if ($Part -eq 'all') {
    Run 'Verification' @('code\verify_all.py')
    Run 'Paper numbers' @('paper\make_paper_numbers.py')
    Run 'Supplementary Fig. S7 (previous qualitative figure)' @('paper\fig6_qualitative_loso.py')
    Run 'Fig. 7 and S7: drawn points = saved predictions' @('paper\check_fig6.py')
    Run 'Graphical abstract' @('paper\graphical_abstract.py')
    Run 'Revision analyses (component evaluation)' @('paper\revision_analyses.py')
    Run 'Figures' @('paper\figures_main.py')
    Run 'Region verification (Supplementary S14)' @('paper\region_verification.py')
    Run 'Margin verification (Supplementary S16, Fig. S6)' @('paper\margin_verification.py')
    Run 'Anatomical figures (Figs. 7 and 8)' @('paper\figures_anatomy.py')
    Run 'Development figure (Fig. 11)' @('paper\fig11_development.py')
    Run 'ITL step-by-step trace (Supplementary Table S18)' @('paper\itl_trace.py')
    Run 'ITL independent checks and S03 leakage test (Supplementary S19)' @('paper\itl_double_check.py')
    Run 'ICP seed table (from the saved runs)' @('paper\icp_seed_sensitivity.py', 'summary')
    Run 'Revision tables' @('paper\revision_tables.py')
    Run 'Reproducibility check' @('paper\reproducibility_check.py')
    Run 'Numbers of the text against the outputs' @('paper\check_text_numbers.py')
    Run 'Every cell of the typed tables against the outputs' @('paper\check_tables.py')
    if ($Pdf) { Run 'Manuscript figures and tables, PDFs' @('paper\build_manuscript_assets.py', '--pdf') }
    else { Run 'Manuscript figures and tables' @('paper\build_manuscript_assets.py') }
    Run 'Unit tests' @('-m', 'unittest', 'discover', '-s', 'tests')
}
Say 'DONE'
Say "Compare the outputs with the committed versions: $Py scripts\compare_outputs.py"
