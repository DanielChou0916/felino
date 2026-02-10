# Contributing to Felino

Thanks for your interest in contributing to Felino! We welcome contributions including bug reports, documentation improvements, new benchmark/tutorial cases, and code changes.

## Before you start
- Please check the documentation: https://danielchou0916.github.io/felino.github.io/
- For questions, open a GitHub Issue with the label **question**.

## Reporting bugs
When reporting a bug, please include:
1. Your OS (Ubuntu/WSL), compiler, and MPI (mpich/openmpi) info
2. MOOSE version (or moose-dev package version) and Felino commit hash
3. The input file(s) needed to reproduce the issue (or a minimal reproduction case)
4. The full error output / stack trace
5. Expected behavior vs. observed behavior

## Suggesting enhancements
- Open an Issue with the label **enhancement**
- Describe the use case, expected input-file interface, and why it is needed
- If possible, propose a minimal example input

## Development setup
Felino is a MOOSE-based application.
1. Follow the installation guide in README (Ubuntu/WSL recommended)
2. Build Felino:
   - `make -j 4`

## Tests and validation
- Run the included installation test:
  - `./felino_test.sh`
- For changes that affect physics or outputs, please include:
  - a new/updated tutorial case, or
  - a small regression-style example and expected output (describe what to check)

## Pull request workflow
1. Fork the repo and create a feature branch:
   - `git checkout -b feature/short-description`
2. Make changes with clear commit messages
3. Update docs/tutorials when user-facing behavior changes
4. Open a Pull Request and describe:
   - what changed
   - why it changed
   - how to test it (commands + expected output)

## License
By contributing, you agree that your contributions will be licensed under the project license (see LICENSE).
