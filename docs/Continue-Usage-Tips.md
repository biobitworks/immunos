---
title: Continue Usage Tips - Working with Local Models
tags: [continue, tips, workarounds]
created: 2025-11-30
---

# Continue Usage Tips

The `@Codebase` feature doesn't always work well with local models. Here are better ways to get accurate answers from Continue about your code.

## ✅ Method 1: Select Code in VS Code (BEST)

1. **Open the file** you want to ask about
2. **Select the code** (highlight it)
3. Press **Cmd+L** to open Continue
4. Ask your question

**Example**:
- Open `qml_ainet.py`
- Select the `_mutate` method (lines 250-274)
- Cmd+L
- Ask: "Explain this code"

✅ Continue will see your selection and respond accurately!

## ✅ Method 2: Paste Code Directly

Just paste the code into Continue chat:

**Example**:
```
Look at this code and explain how it works:

[paste your code here]

Why does it use uniform random instead of Gaussian mutation?
```

## ✅ Method 3: Reference Line Numbers

If you have a file open:

```
In the currently open file, explain the _mutate method starting at line 250
```

## ✅ Method 4: Ask About Open Files

Continue automatically sees files you have open in tabs. Ask:

```
Look at the open qml_ainet.py file. How does the mutation operator work?
```

## ❌ What Doesn't Work Well

These features don't work reliably with local Ollama models:

- `@Codebase` - Tries to index but often fails
- `@File` syntax - Hit or miss
- File path references - Often ignored

## 🎯 Best Practices

### For Code Explanation
1. Open file
2. Select specific function/class
3. Cmd+L
4. "Explain this"

### For Code Generation
```
Generate a function that does X

Use this style:
[paste example code from your codebase]
```

### For Debugging
1. Select the buggy code
2. Cmd+L
3. "Why isn't this working? What's wrong?"

### For Refactoring
1. Select code to refactor
2. Cmd+I (inline edit)
3. "Refactor this to use dataclasses"
4. Review changes

## 🔥 Pro Tips

### Tip 1: Give Context
Instead of:
```
How does mutation work?
```

Do this:
```
I'm looking at QML-AiNet which is for discrete constraint spaces.
Here's the mutation method:

[paste _mutate method]

Why does it use random.randint instead of Gaussian noise?
```

### Tip 2: Reference Your Documentation
```
I have this code:
[paste code]

And this documentation:
[paste from your Obsidian notes]

Are they consistent? Is the code doing what the docs say?
```

### Tip 3: Multi-step Questions
```
First, look at this code:
[paste _mutate from qml_ainet.py]

Now compare it to this:
[paste _mutate from opt_ainet.py]

What's the key difference and why?
```

### Tip 4: Use Custom Commands
Select code, then in chat:
- `/test` - Generate tests
- `/docstring` - Add docs
- `/security` - Security review
- `/optimize` - Performance tips

## 🎨 Example Conversations

### Good Conversation 1: Code Review
```
You: [Select _mutate method, Cmd+L]
     "Review this code for potential bugs"

Continue: [Analyzes the selected code]
          "This looks good but consider..."

You: "How would I add logging to track mutations?"

Continue: [Suggests logging additions]
```

### Good Conversation 2: Comparison
```
You: Here's QML-AiNet mutation:
     [paste code from qml_ainet.py:250-274]

     Here's Opt-AiNet mutation:
     [paste code from opt_ainet.py]

     Explain the key difference

Continue: [Accurate comparison because it has both codes]
```

### Good Conversation 3: Implementation Help
```
You: I want to add network suppression to my B Cell agent.

     Here's the current B Cell code:
     [paste relevant section]

     Here's how network suppression works in Opt-AiNet:
     [paste network_suppression method]

     How should I integrate this?

Continue: [Provides specific suggestions based on your actual code]
```

## 📊 Model Selection

### Use Qwen 2.5 Coder (7B) for:
- ✅ Code explanation
- ✅ Code generation
- ✅ Refactoring suggestions
- ✅ Adding documentation

### Use DeepSeek R1 (14B) for:
- ✅ Complex debugging
- ✅ Architectural decisions
- ✅ "Why" questions
- ✅ Comparing approaches

### Use Qwen Fast (1.5B) for:
- ✅ Autocomplete (automatic)
- ✅ Simple questions
- ✅ Quick clarifications

## 🚨 Common Issues

### Issue: Generic/Wrong Answers
**Cause**: Continue doesn't see your code
**Fix**: Select code first OR paste it into chat

### Issue: Hallucinated File Paths
**Cause**: LLM making up files
**Fix**: Paste actual code, don't rely on file references

### Issue: Outdated Information
**Cause**: Model trained on old data
**Fix**: Provide current code in the prompt

### Issue: Slow Responses
**Cause**: Large context or slow model
**Fix**: Switch to Qwen Fast (1.5B) or reduce context

## 🎓 Learning Curve

**Day 1**: Use select + Cmd+L exclusively
**Day 2**: Try custom commands (/test, /docstring)
**Day 3**: Experiment with pasting code for comparisons
**Week 1**: Build your personal prompting style

## 🔗 Related

- [[VS-Code-LLM-Integration]] - Full setup guide
- [[VS-Code-LLM-Quick-Reference]] - Keyboard shortcuts
- Config: `~/.continue/config.yaml`

---

**Last Updated**: 2025-11-30
**Continue Version**: 1.2.11
**Key Insight**: Selection + Cmd+L works better than @Codebase with local models
