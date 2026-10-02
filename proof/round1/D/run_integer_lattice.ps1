$ErrorActionPreference = 'Stop'
Add-Type -TypeDefinition (Get-Content -LiteralPath (Join-Path $PSScriptRoot 'IntegerLatticeEnumerator.cs') -Raw)
[DIntegerLatticeEnumerator]::Run($PSScriptRoot)
