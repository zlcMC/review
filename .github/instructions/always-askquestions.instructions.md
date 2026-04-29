---
description: "Call vscode_askQuestions after replies, and continue when the tool returns an answer."
applyTo: "**"
---

# Continuous AskQuestions Interaction Loop

For this user, assistant turns should normally call the `vscode_askQuestions` tool after substantive replies to keep interaction moving. The tool is an interactive checkpoint, not a forced turn-ending handoff; if it returns a selected option or free text, continue handling that answer in the same turn.

Rules:

- Use `vscode_askQuestions` after substantive answers, confirmations, transition messages, meta discussion, or tool-result summaries.
- Do not force the tool for brief user messages such as "继续", "好的", "收到", or "嗯" if a direct continuation is more useful.
- A turn does not always need a full long answer before `vscode_askQuestions`: if important information, confirmation, or a user choice is needed, give a concise explanation and ask with `vscode_askQuestions` immediately.
- If `vscode_askQuestions` returns a selected option or free text, treat that answer as the user's current request and continue answering or acting on it in the same turn.
- Never respond to a selected option with only an empty transition such as "下一步会做..." or "我接下来会...". Continue the requested explanation or task in that same reply unless more user input is genuinely required.
- If there is no useful project-specific follow-up, ask a lightweight continuation question such as whether to continue, pause, or change direction.
- The only exception is when the user explicitly says not to ask follow-up questions anymore.

Continuous interaction rules:

- After calling `vscode_askQuestions`, if the tool returns selected answers or free text, continue with the requested work or response in the same turn.
- After handling an answer returned by `vscode_askQuestions`, if you produce any visible response, call `vscode_askQuestions` again before ending the turn. Do not end just because there is no natural next step.
- Do not send an empty `final` after the tool. If the turn ends, provide a meaningful concise conclusion.
- Avoid infinite ask loops by asking a lightweight continuation question such as whether to continue, pause, or change direction. Do not skip `vscode_askQuestions` after a visible response.
