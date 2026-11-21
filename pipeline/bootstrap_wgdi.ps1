[CmdletBinding()]
param(
    [string]$PythonExe = "",
    [string]$EnvDir = "",
    [switch]$Force,
    [switch]$SkipExternalCheck,
    [switch]$SkipPal2nal,
    [switch]$AllowUnpinnedFallback
)

$ErrorActionPreference = "Stop"
$env:PYTHONUTF8 = "1"
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$python = if ([string]::IsNullOrWhiteSpace($PythonExe)) {
    (Get-Command python -ErrorAction Stop).Source
}
else {
    (Resolve-Path -LiteralPath $PythonExe).Path
}
$arguments = @(
    (Join-Path $PSScriptRoot "code/tools/bootstrap_wgdi.py"),
    "--project-root", $PSScriptRoot
)
if (-not [string]::IsNullOrWhiteSpace($EnvDir)) { $arguments += @("--env-dir", $EnvDir) }
if (-not [string]::IsNullOrWhiteSpace($PythonExe)) { $arguments += @("--python-exe", $PythonExe) }
if ($Force) { $arguments += "--force" }
if ($SkipExternalCheck) { $arguments += "--skip-external-check" }
if ($SkipPal2nal) { $arguments += "--skip-pal2nal" }
if ($AllowUnpinnedFallback) { $arguments += "--allow-unpinned-fallback" }
& $python @arguments
exit $LASTEXITCODE
