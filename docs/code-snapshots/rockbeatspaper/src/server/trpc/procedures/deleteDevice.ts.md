---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/deleteDevice.ts
relative: rockbeatspaper/src/server/trpc/procedures/deleteDevice.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

export const deleteDevice = baseProcedure
  .input(
    z.object({
      authToken: z.string(),
      deviceId: z.number(),
    })
  )
  .mutation(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    // Verify device belongs to user
    const device = await db.device.findUnique({
      where: { id: input.deviceId },
    });

    if (!device || device.userId !== userId) {
      throw new TRPCError({
        code: "NOT_FOUND",
        message: "Device not found",
      });
    }

    await db.device.delete({
      where: { id: input.deviceId },
    });

    return { success: true };
  });

```
