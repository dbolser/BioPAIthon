# CLAUDE.md

**Read [AGENTS.md](AGENTS.md).** It is the single source of truth for working
on BioPAIthon, and it applies to you exactly as it applies to every other
contributor.

This file exists only so that Claude Code finds its way there. Everything that
was once here — project layout, setup, how to run the tests, style rules, the
standard a contribution has to meet — now lives in `AGENTS.md`, so that humans,
computational intelligences, AIs and chimps all read the same instructions
rather than drifting apart into per-tool dialects.

Do not add project guidance to this file. Put it in `AGENTS.md`.

Two things from there worth repeating, because they are the ones most often
skipped:

- The importable package is still `Bio`. The fork renamed the project, not the
  module. `import Bio` must keep working.
- Check before you assume. This codebase is old and much of what looks like a
  mistake is deliberate. A fluent, confident, wrong patch is the most expensive
  thing you can submit here.
