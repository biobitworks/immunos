---
source: /Users/byron/projects/rockbeatspaper/src/server/utils/auth.ts
relative: rockbeatspaper/src/server/utils/auth.ts
generated_at: 2025-12-23 10:28
---

```typescript
import jwt from "jsonwebtoken";
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import { env } from "~/server/env";

const tokenPayloadSchema = z.object({
  userId: z.number(),
});

export function verifyAuthToken(authToken: string): { userId: number } {
  try {
    const verified = jwt.verify(authToken, env.JWT_SECRET);
    const parsed = tokenPayloadSchema.parse(verified);
    return parsed;
  } catch (error) {
    throw new TRPCError({
      code: "UNAUTHORIZED",
      message: "Invalid or expired authentication token",
    });
  }
}

export function signAuthToken(userId: number): string {
  return jwt.sign({ userId }, env.JWT_SECRET, { expiresIn: "30d" });
}

```
