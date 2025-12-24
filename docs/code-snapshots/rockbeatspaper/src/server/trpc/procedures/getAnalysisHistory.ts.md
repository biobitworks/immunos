---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/getAnalysisHistory.ts
relative: rockbeatspaper/src/server/trpc/procedures/getAnalysisHistory.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

export const getAnalysisHistory = baseProcedure
  .input(z.object({ authToken: z.string() }))
  .query(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    const analyses = await db.thermalAnalysis.findMany({
      where: { userId },
      orderBy: { createdAt: "desc" },
      take: 10,
    });

    return {
      analyses: analyses.map((a) => ({
        id: a.id,
        analysisType: a.analysisType,
        recommendations: JSON.parse(a.recommendations),
        potentialSavings: a.potentialSavings,
        createdAt: a.createdAt,
      })),
    };
  });

```
