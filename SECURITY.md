# Security policy

## Supported versions

Only the latest release is supported.

## What the tool does with data

Genotype files are read locally, or inside the browser tab in the Shinylive site. The only network use is the opt-in BrAPI loader on the command line, which sends the token only in the `Authorization` header (docs/adr/0024).

## Reporting a vulnerability

Use GitHub private vulnerability reporting on this repository (the Security tab, "Report a vulnerability"). Please do not open a public issue for a vulnerability. Expect a response within 14 days. There is no bug bounty.
