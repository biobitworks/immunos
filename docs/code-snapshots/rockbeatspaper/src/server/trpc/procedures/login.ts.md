---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/login.ts
relative: rockbeatspaper/src/server/trpc/procedures/login.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import bcryptjs from "bcryptjs";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { signAuthToken } from "~/server/utils/auth";

export const login = baseProcedure
  .input(
    z.object({
      email: z.string().email(),
      password: z.string(),
    })
  )
  .mutation(async ({ input }) => {
    // Find user
    const user = await db.user.findUnique({
      where: { email: input.email },
    });

    if (!user) {
      throw new TRPCError({
        code: "UNAUTHORIZED",
        message: "Invalid email or password",
      });
    }

    // Verify password
    const isValidPassword = await bcryptjs.compare(
      input.password,
      user.hashedPassword
    );

    if (!isValidPassword) {
      throw new TRPCError({
        code: "UNAUTHORIZED",
        message: "Invalid email or password",
      });
    }

    // Generate auth token
    const authToken = signAuthToken(user.id);

    return {
      authToken,
      user: {
        id: user.id,
        email: user.email,
        name: user.name,
      },
    };
  });

```
