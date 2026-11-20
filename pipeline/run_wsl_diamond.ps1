[CmdletBinding()]
param(
    [string]$Distro = "Ubuntu",
    [ValidateRange(1, 128)]
    [int]$Threads = 8,
    [switch]$TestOnly,
    [string]$ProjectRoot = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $PSScriptRoot
}
else {
    (Resolve-Path -LiteralPath $ProjectRoot).Path
}
$windowsProjectPath = $projectRoot.Replace("\", "/")
$wslProjectRoot = (
    & wsl -d $Distro -- wslpath -a -u $windowsProjectPath |
        Select-Object -First 1
).Trim()

if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($wslProjectRoot)) {
    throw "无法将项目路径转换为 WSL 路径：$projectRoot"
}

& wsl -d $Distro -- /usr/local/bin/diamond version
if ($LASTEXITCODE -ne 0) {
    throw "WSL 中的官方 DIAMOND 不可用：/usr/local/bin/diamond"
}

if ($TestOnly) {
    & wsl -d $Distro -- bash "$wslProjectRoot/scripts/test_wsl_diamond.sh"
    exit $LASTEXITCODE
}

& wsl -d $Distro -- env `
    "DIAMOND=/usr/local/bin/diamond" `
    "DIAMOND_THREADS=$Threads" `
    bash "$wslProjectRoot/1run_diamond.sh"
exit $LASTEXITCODE
