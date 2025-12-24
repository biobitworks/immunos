---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/getDevices.ts
relative: rockbeatspaper/src/server/trpc/procedures/getDevices.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

export const getDevices = baseProcedure
  .input(z.object({ authToken: z.string() }))
  .query(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    const devices = await db.device.findMany({
      where: { userId },
      orderBy: { createdAt: "desc" },
    });

    return { devices };
  });

```
