$ErrorActionPreference = "Stop"

$amrex = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$source = Get-Content -LiteralPath (Join-Path $amrex "Src/EB/AMReX_EB_StateRedistItracker.cpp") -Raw

function Assert-Matches([string]$Text, [string]$Pattern, [string]$Message) {
    if ($Text -notmatch $Pattern) {
        throw $Message
    }
}

Assert-Matches $source 'StateRedistCandidateIsValid' `
    "MakeITracker must validate every candidate against the physical domain and vfrac."
Assert-Matches $source 'StateRedistFindValidCandidate' `
    "MakeITracker must replace an invalid preferred direction with a legal candidate."
Assert-Matches $source 'vfrac\(ii,jj,kk\)\s*>\s*Real\(0\.0\)' `
    "A StateRedist candidate must have positive fluid volume."

# Reproduce the ylo failure pattern.  The first x+ neighbor leaves the merged
# volume slightly below 0.5.  The normal-directed second candidate is y-, but
# j=-1 is outside this non-periodic domain.  A legal y+ neighbor must be used.
$target = 0.5
$sumVol = 0.499177
$domainJLo = 0
$cellJ = 0
$preferredSecondJOffset = -1
$fallbackSecondJOffset = 1

if (($cellJ + $preferredSecondJOffset) -ge $domainJLo) {
    throw "The regression setup must make the preferred y- candidate illegal."
}
if (($cellJ + $fallbackSecondJOffset) -lt $domainJLo) {
    throw "The fallback y+ candidate must remain inside the domain."
}

$fallbackVfrac = 1.0
if (($sumVol + $fallbackVfrac) -lt $target) {
    throw "The legal y+ fallback must satisfy the target volume threshold."
}

Write-Host "StateRedist ylo candidate-filter regression checks passed."
