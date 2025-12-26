# Agent Protocol (Required)

This protocol standardizes context sharing, logging, snapshots, and daily notes for all AI agents working in `/Users/byron/projects`.

## 1) Startup (every session)
1. Read shared context:
   - `/Users/byron/projects/claude.md`
   - `/Users/byron/projects/.immunos/recovery/CONTEXT_RECOVERY.md`
2. Read model-specific context (if applicable):
   - `/Users/byron/projects/.immunos/model-contexts/<model>-context.md`
3. Check today’s journal(s):
   - `/Users/byron/projects/.immunos/journal/YYYY-MM-DD.md`
   - `/Users/byron/projects/daily/YYYY-MM-DD.md`

## 2) During Work (ongoing)
- Log key decisions and changes:
  - `python3 /Users/byron/projects/scripts/immunos_memory.py store --content "<summary>" --priority high --type conversation`
- Capture snapshots at major milestones:
  - `python3 /Users/byron/projects/scripts/immunos_snapshot.py create --trigger manual --summary "<what changed>"`

## 3) Session Logging (required)
- Log session summary (always):
  - `python3 /Users/byron/projects/scripts/immunos_log_session.py --model "<model>" --summary "<work summary>" --tags "<tags>" --files "<touched files>"`

## 4) Daily Notes (required)
- Append a bullet to BOTH:
  - `/Users/byron/projects/.immunos/journal/YYYY-MM-DD.md`
  - `/Users/byron/projects/daily/YYYY-MM-DD.md`
- Format: one concise bullet describing what this agent did.

## 5) Shutdown (end of session)
1. Final snapshot (if meaningful progress):
   - `python3 /Users/byron/projects/scripts/immunos_snapshot.py create --trigger daily --summary "<end-of-session checkpoint>"`
2. Confirm journals updated.

## 6) Privacy Rules
- Never commit secrets or PII.
- Store personal data in `/Users/byron/projects/claude.private.md` (gitignored).
- Keep `.immunos/` local-only unless explicitly approved.

## 7) Exceptions
- If automation is blocked or missing, document the omission in the daily note.

