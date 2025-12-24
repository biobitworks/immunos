---
source: /Users/byron/projects/rockbeatspaper/src/utils/cn.ts
relative: rockbeatspaper/src/utils/cn.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}

```
