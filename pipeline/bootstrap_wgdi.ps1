<#
.SYNOPSIS
    Create a local WGDI runtime for this workflow repository.

.DESCRIPTION
    Creates .wgdi-env next to this script and installs the Python packages WGDI
    needs.  The launcher bin/wgdi detects .wgdi-env automatically, so the
    repository can be copied to another machine and executed without editing
    code or hard-coded interpreter paths.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1
    powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1 -EnvDir D:\wgdi-venv
#>
[CmdletBinding()]
param(
    [string]$PythonExe = "",
    [string]$EnvDir = "",
    [switch]$Force,
    [switch]$SkipExternalCheck
)

$ErrorActionPreference = "Stop"
$projectRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($EnvDir)) {
    $EnvDir = Join-Path $projectRoot ".wgdi-env"
}
$requirements = Join-Path $projectRoot "requirements-wgdi-validated.txt"

function Resolve-Python {
    param([string]$Requested)
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        return (Resolve-Path -LiteralPath $Requested).Path
    }
    foreach ($name in @("python", "python3", "py")) {
        $command = Get-Command $name -ErrorAction SilentlyContinue
        if ($null -ne $command) {
            return $command.Source
        }
    }
    throw "找不到 Python。请安装 Python 3.11+ 或用 -PythonExe 指定解释器。"
}

$python = Resolve-Python $PythonExe
Write-Host "Python: $python"

if ((Test-Path -LiteralPath $EnvDir) -and $Force) {
    Write-Host "删除已有环境: $EnvDir"
    Remove-Item -LiteralPath $EnvDir -Recurse -Force
}
if (-not (Test-Path -LiteralPath $EnvDir)) {
    & $python -m venv $EnvDir
}
$venvPython = Join-Path $EnvDir "Scripts\python.exe"
if (-not (Test-Path -LiteralPath $venvPython)) {
    $venvPython = Join-Path $EnvDir "bin/python"
}
if (-not (Test-Path -LiteralPath $venvPython)) {
    throw "虚拟环境创建失败：$EnvDir"
}

& $venvPython -m pip install --upgrade pip

$installed = $false
if (Test-Path -LiteralPath $requirements) {
    Write-Host "安装固定版本依赖: $requirements"
    & $venvPython -m pip install -r $requirements
    $installed = ($LASTEXITCODE -eq 0)
}
if (-not $installed) {
    Write-Host "固定版本安装失败，改用 WGDI 官方依赖（不锁定版本）。" -ForegroundColor Yellow
    & $venvPython -m pip install wgdi matplotlib numpy scipy pandas biopython tabulate
    if ($LASTEXITCODE -ne 0) {
        throw "WGDI Python 依赖安装失败。"
    }
}

& $venvPython -c "import wgdi, wgdi.run; print('wgdi import OK')"
if ($LASTEXITCODE -ne 0) {
    throw "WGDI 安装校验失败。"
}

Write-Host ""
Write-Host "Python 运行时已就绪: $venvPython" -ForegroundColor Green
Write-Host "启动方式: bin/wgdi <wgdi 参数>（自动使用 $EnvDir）"

if (-not $SkipExternalCheck) {
    Write-Host ""
    Write-Host "外部程序（DIAMOND/BLAST/MUSCLE/yn00）需另行安装并在 PATH 中可用："
    foreach ($tool in @("diamond", "blastp", "makeblastdb", "muscle", "yn00", "perl")) {
        $found = Get-Command $tool -ErrorAction SilentlyContinue
        if ($null -ne $found) {
            Write-Host ("  [OK]   {0} -> {1}" -f $tool, $found.Source)
        } else {
            Write-Host ("  [MISS] {0}" -f $tool) -ForegroundColor Yellow
        }
    }
}

Write-Host ""
Write-Host "自检: powershell -File .\run_all_pipeline.ps1 -CheckOnly -SkipPreflight"
