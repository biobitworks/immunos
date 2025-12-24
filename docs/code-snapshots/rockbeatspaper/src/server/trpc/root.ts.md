---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/root.ts
relative: rockbeatspaper/src/server/trpc/root.ts
generated_at: 2025-12-23 10:28
---

```typescript
import {
  createCallerFactory,
  createTRPCRouter,
} from "~/server/trpc/main";
import { signup } from "~/server/trpc/procedures/signup";
import { login } from "~/server/trpc/procedures/login";
import { getMe } from "~/server/trpc/procedures/getMe";
import { getDevices } from "~/server/trpc/procedures/getDevices";
import { createDevice } from "~/server/trpc/procedures/createDevice";
import { deleteDevice } from "~/server/trpc/procedures/deleteDevice";
import { analyzeThermalData } from "~/server/trpc/procedures/analyzeThermalData";
import { getAnalysisHistory } from "~/server/trpc/procedures/getAnalysisHistory";
import { runPiBenfordAnalysis } from "~/server/trpc/procedures/runPiBenfordAnalysis";

export const appRouter = createTRPCRouter({
  signup,
  login,
  getMe,
  getDevices,
  createDevice,
  deleteDevice,
  analyzeThermalData,
  getAnalysisHistory,
  runPiBenfordAnalysis,
});

export type AppRouter = typeof appRouter;

export const createCaller = createCallerFactory(appRouter);

```
