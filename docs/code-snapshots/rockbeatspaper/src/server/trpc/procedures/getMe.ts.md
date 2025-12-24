---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/getMe.ts
relative: rockbeatspaper/src/server/trpc/procedures/getMe.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

export const getMe = baseProcedure
  .input(z.object({ authToken: z.string() }))
  .query(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    const user = await db.user.findUnique({
      where: { id: userId },
    });

    if (!user) {
      throw new TRPCError({
        code: "NOT_FOUND",
        message: "User not found",
      });
    }

    return {
      id: user.id,
      email: user.email,
      name: user.name,
    };
  });

```
