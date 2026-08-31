---
description: Generate a commit message for staged changes
---

Generate a concise, informative git commit message following these guidelines:

## Format

```
<type>(<scope>): <subject>

<body>
```

## Rules

1. **Type** (required): One of:
   - `feat`: New feature
   - `fix`: Bug fix
   - `refactor`: Code change that neither fixes a bug nor adds a feature
   - `docs`: Documentation only
   - `test`: Adding or updating tests
   - `chore`: Build process, dependencies, or tooling
   - `perf`: Performance improvement
   - `style`: Formatting, whitespace (no code change)

2. **Scope** (optional): Component or area affected

3. **Subject** (required):
   - Use imperative mood ("add" not "added")
   - No period at the end
   - Max 50 characters
   - Lowercase first letter

4. **Body** (optional):
   - Wrap at 300 characters
   - Explain *what* and *why*, not *how*
   - Separate from subject with blank line

## Instructions

1. Run `git diff --cached` to see staged changes
2. Analyze the changes and determine the appropriate type and scope
3. Write a clear subject line summarizing the change
4. Add body only if the change needs explanation

Output ONLY the commit message. Then, ask if you should commit the changes. If the user replies 'yes', proceed to run `git commit` with the generated message.
