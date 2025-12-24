---
source: /Users/byron/projects/rockbeatspaper/src/routes/__root.tsx
relative: rockbeatspaper/src/routes/__root.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import {
  Outlet,
  createRootRoute,
  useRouterState,
} from "@tanstack/react-router";
import { TRPCReactProvider } from "~/trpc/react";
import { Toaster } from "react-hot-toast";

export const Route = createRootRoute({
  component: RootComponent,
});

function RootComponent() {
  const isLoading = useRouterState({ select: (s) => s.isLoading });

  if (isLoading) {
    return <div>Loading...</div>;
  }

  return (
    <TRPCReactProvider>
      <Toaster position="top-right" />
      <Outlet />
    </TRPCReactProvider>
  );
}

```
