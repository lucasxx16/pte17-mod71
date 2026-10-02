$ErrorActionPreference='Stop'
Add-Type -TypeDefinition (Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'audit_isqrt.cs')) -CompilerOptions '/checked+'
[EAuditIsqrt]::Run((Join-Path $PSScriptRoot 'isqrt_audit_result.json'))
