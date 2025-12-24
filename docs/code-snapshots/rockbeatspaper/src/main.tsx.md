---
source: /Users/byron/projects/rockbeatspaper/src/main.tsx
relative: rockbeatspaper/src/main.tsx
generated_at: 2025-12-23 10:28
---

```tsx
/// <reference types="vinxi/types/client" />

import React from "react";
import ReactDOM from "react-dom/client";
import { RouterProvider } from "@tanstack/react-router";

import "./styles.css";

import { createRouter } from "./router";

// Set up a Router instance
const router = createRouter();

const rootElement = document.getElementById("root")!;
if (!rootElement.innerHTML) {
  const root = ReactDOM.createRoot(rootElement);
  root.render(
    <React.StrictMode>
      <RouterProvider router={router} />
    </React.StrictMode>,
  );
}

```
