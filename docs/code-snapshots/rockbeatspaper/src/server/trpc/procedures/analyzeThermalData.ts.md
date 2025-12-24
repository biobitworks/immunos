---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/analyzeThermalData.ts
relative: rockbeatspaper/src/server/trpc/procedures/analyzeThermalData.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { generateObject } from "ai";
import { openai } from "@ai-sdk/openai";
import { db } from "~/server/db";
import { baseProcedure } from "~/server/trpc/main";
import { verifyAuthToken } from "~/server/utils/auth";

const analysisSchema = z.object({
  summary: z.string(),
  recommendations: z.array(
    z.object({
      title: z.string(),
      description: z.string(),
      priority: z.enum(["high", "medium", "low"]),
    })
  ),
  estimatedSavings: z.number(),
  insights: z.array(z.string()),
});

export const analyzeThermalData = baseProcedure
  .input(
    z.object({
      authToken: z.string(),
      currentTemperature: z.number(),
      targetTemperature: z.number(),
      humidity: z.number().optional(),
      deviceType: z.string(),
      location: z.string(),
      currentPowerUsage: z.number(),
    })
  )
  .mutation(async ({ input }) => {
    const { userId } = verifyAuthToken(input.authToken);

    const model = openai("gpt-4o");

    const prompt = `Analyze the following thermal/energy data and provide optimization recommendations:
    
Current Temperature: ${input.currentTemperature}°C
Target Temperature: ${input.targetTemperature}°C
${input.humidity ? `Humidity: ${input.humidity}%` : ""}
Device Type: ${input.deviceType}
Location: ${input.location}
Current Power Usage: ${input.currentPowerUsage}W

Provide actionable recommendations to optimize energy efficiency and thermal comfort. Include estimated energy savings as a percentage.`;

    const { object } = await generateObject({
      model,
      schema: analysisSchema,
      prompt,
    });

    // Save analysis to database
    const analysis = await db.thermalAnalysis.create({
      data: {
        userId,
        analysisType: "optimization",
        inputData: JSON.stringify(input),
        recommendations: JSON.stringify(object.recommendations),
        potentialSavings: object.estimatedSavings,
      },
    });

    return {
      analysisId: analysis.id,
      summary: object.summary,
      recommendations: object.recommendations,
      estimatedSavings: object.estimatedSavings,
      insights: object.insights,
    };
  });

```
