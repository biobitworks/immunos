---
source: /Users/byron/projects/rockbeatspaper/src/components/ProtectedLayout.tsx
relative: rockbeatspaper/src/components/ProtectedLayout.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { ReactNode, useEffect } from "react";
import { useNavigate } from "@tanstack/react-router";
import { useAuthStore } from "~/stores/auth";
import { Navbar } from "./Navbar";

interface ProtectedLayoutProps {
  children: ReactNode;
}

export function ProtectedLayout({ children }: ProtectedLayoutProps) {
  const { authToken } = useAuthStore();
  const navigate = useNavigate();

  useEffect(() => {
    if (!authToken) {
      navigate({ to: "/login" });
    }
  }, [authToken, navigate]);

  if (!authToken) {
    return null;
  }

  return (
    <div className="min-h-screen bg-gray-50 dark:bg-gray-900">
      <Navbar />
      <main>{children}</main>
    </div>
  );
}

```
