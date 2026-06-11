# Contributing

Contributions are welcome when they improve reproducibility, portability, or
clarity of the analysis workflow.

## Good Contributions

- fix a script so it runs on a clean machine
- document missing package versions or command-line steps
- add a small synthetic test input for a workflow
- improve README, citation, or data-availability notes
- report a reproducibility issue with the exact command and error message

## Please Do Not Add

- private sample metadata
- raw sequencing files or large generated outputs
- unpublished manuscript drafts
- institution-only paths or credentials
- personally identifying information

## Pull Request Checklist

1. Keep changes focused on one reproducibility problem.
2. Use synthetic or public example inputs when possible.
3. Confirm that no private data or large binary files are staged.
4. Mention which scripts or notes you tested.
5. Link any related issue.

## Local Review

```powershell
git status --short
rg -n -i "password|token|secret|credential|private|sample_id" .
```
