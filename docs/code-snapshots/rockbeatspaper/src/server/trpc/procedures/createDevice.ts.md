---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/createDevice.ts
relative: rockbeatspaper/src/server/trpc/procedures/createDevice.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

export const createDevice = baseProcedure
  .input(
    z.object({
      authToken: z.string(),
      name: z.string().min(1),
      type: z.enum(["thermostat", "hvac", "heater", "cooler"]),
      location: z.string().min(1),
    })
  )
  .mutation(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    const device = await db.device.create({
      data: {
        userId,
        name: input.name,
        type: input.type,
        location: input.location,
        status: "active",
      },
    });

    return { device };
  });

```
