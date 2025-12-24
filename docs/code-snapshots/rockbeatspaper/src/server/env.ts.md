---
source: /Users/byron/projects/rockbeatspaper/src/server/env.ts
relative: rockbeatspaper/src/server/env.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";

const envSchema = z.object({
  NODE_ENV: z.enum(["development", "production"]),
  BASE_URL: z.string().optional(),
  BASE_URL_OTHER_PORT: z.string().optional(),
  ADMIN_PASSWORD: z.string(),
  JWT_SECRET: z.string(),
  OPENAI_API_KEY: z.string(),
});

export const env = envSchema.parse(process.env);

```
