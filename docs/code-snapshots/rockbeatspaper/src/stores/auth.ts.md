---
source: /Users/byron/projects/rockbeatspaper/src/stores/auth.ts
relative: rockbeatspaper/src/stores/auth.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { create } from "zustand";
import { persist, createJSONStorage } from "zustand/middleware";

type AuthStore = {
  authToken: string | null;
  user: {
    id: number;
    email: string;
    name: string;
  } | null;
  setAuth: (authToken: string, user: { id: number; email: string; name: string }) => void;
  clearAuth: () => void;
};

export const useAuthStore = create<AuthStore>()(
  persist(
    (set) => ({
      authToken: null,
      user: null,
      setAuth: (authToken, user) => set({ authToken, user }),
      clearAuth: () => set({ authToken: null, user: null }),
    }),
    {
      name: "thermal-auth-storage",
      storage: createJSONStorage(() => localStorage),
    }
  )
);

```
