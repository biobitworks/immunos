---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/signup.ts
relative: rockbeatspaper/src/server/trpc/procedures/signup.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import bcryptjs from "bcryptjs";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { signAuthToken } from "~/server/utils/auth";

export const signup = baseProcedure
  .input(
    z.object({
      email: z.string().email(),
      name: z.string().min(1),
      password: z.string().min(8),
    })
  )
  .mutation(async ({ input }) => {
    // Check if user already exists
    const existingUser = await db.user.findUnique({
      where: { email: input.email },
    });

    if (existingUser) {
      throw new TRPCError({
        code: "CONFLICT",
        message: "User with this email already exists",
      });
    }

    // Hash password
    const hashedPassword = await bcryptjs.hash(input.password, 10);

    // Create user
    const user = await db.user.create({
      data: {
        email: input.email,
        name: input.name,
        hashedPassword,
      },
    });

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
