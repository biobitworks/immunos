---
source: /Users/byron/projects/rockbeatspaper/src/router.tsx
relative: rockbeatspaper/src/router.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createRouter as createTanStackRouter } from "@tanstack/react-router";
import { routeTree } from "./generated/tanstack-router/routeTree.gen";

export function createRouter() {
  const router = createTanStackRouter({
    routeTree,
    scrollRestoration: true,
    defaultPreload: "intent",
    defaultPendingComponent: () => <div>Loading...</div>,
  });

  return router;
}

declare module "@tanstack/react-router" {
  interface Register {
    router: ReturnType<typeof createRouter>;
  }
}

```
