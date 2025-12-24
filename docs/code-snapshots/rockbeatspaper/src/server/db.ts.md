---
source: /Users/byron/projects/rockbeatspaper/src/server/db.ts
relative: rockbeatspaper/src/server/db.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { PrismaClient } from "@prisma/client";

import { env } from "~/server/env";

const createPrismaClient = () =>
  new PrismaClient({
    log:
      env.NODE_ENV === "development" ? ["query", "error", "warn"] : ["error"],
  });

const globalForPrisma = globalThis as unknown as {
  prisma: ReturnType<typeof createPrismaClient> | undefined;
};

export const db = globalForPrisma.prisma ?? createPrismaClient();

if (env.NODE_ENV !== "production") globalForPrisma.prisma = db;

```
