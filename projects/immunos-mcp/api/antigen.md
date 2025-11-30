---
project: immunos-mcp
type: api-documentation
source: src/core/antigen.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# antigen

**Source**: [[../code-mirror/src/core/antigen.py|src/core/antigen.py]]

Antigen representation for IMMUNOS-MCP

An antigen is a data unit to be processed by the immune system,
representing text, code, structured data, or other inputs.

## Contents

### Classes
- [DataType](#datatype)
- [Antigen](#antigen)
- [AntigenBatch](#antigenbatch)

## Classes

### DataType

**Inherits from**: `Enum`

**Source**: `src/core/antigen.py:13`

Type of data contained in the antigen

---

### Antigen

**Source**: `src/core/antigen.py:23`

Represents a data unit to be classified or analyzed.

Attributes:
    data: The raw data (text, code, dict, etc.)
    data_type: Type of data
    features: Extracted feature vector
    metadata: Additional context (sender, timestamp, etc.)
    class_label: Known classification (for training)
    identifier: Unique ID

#### Methods

##### `from_text(`cls, text: str, class_label: Optional[str], metadata: Optional[Dict])` → `Antigen`

Create antigen from text.

Args:
    text: Text content
    class_label: Optional known classification
    metadata: Additional context

Returns:
    Antigen instance

##### `from_code(`cls, code: str, language: str, class_label: Optional[str], metadata: Optional[Dict])` → `Antigen`

Create antigen from code.

Args:
    code: Source code
    language: Programming language
    class_label: Optional known classification
    metadata: Additional context

Returns:
    Antigen instance

##### `from_dict(`cls, data: Dict[(str, Any)], class_label: Optional[str], metadata: Optional[Dict])` → `Antigen`

Create antigen from structured data.

Args:
    data: Structured data dictionary
    class_label: Optional known classification
    metadata: Additional context

Returns:
    Antigen instance

##### `set_features(`self, features: Dict[(str, Any)])`

Set extracted features

##### `add_metadata(`self, key: str, value: Any)`

Add metadata field

##### `get_text_content(`self)` → `str`

Get text representation of data

##### `__repr__(`self)` → `str`

*No documentation available.*

---

### AntigenBatch

**Source**: `src/core/antigen.py:134`

Batch of antigens for efficient processing

#### Methods

##### `__len__(`self)` → `int`

*No documentation available.*

##### `__iter__(`self)`

*No documentation available.*

##### `__getitem__(`self, idx: int)` → `Antigen`

*No documentation available.*

##### `from_texts(`cls, texts: List[str], labels: Optional[List[str]])` → `AntigenBatch`

Create batch from list of texts

---

## Links

- [[../code-mirror/src/core/antigen.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

