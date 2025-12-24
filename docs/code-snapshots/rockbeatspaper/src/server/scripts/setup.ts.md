---
source: /Users/byron/projects/rockbeatspaper/src/server/scripts/setup.ts
relative: rockbeatspaper/src/server/scripts/setup.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { execSync } from "child_process";
import path from "path";

async function setup() {
  console.log("Installing Python dependencies...");
  try {
    const requirementsPath = path.join(process.cwd(), "requirements.txt");
    execSync(`pip3 install -r ${requirementsPath}`, {
      stdio: "inherit",
      cwd: process.cwd(),
    });
    console.log("Python dependencies installed successfully");
  } catch (error) {
    console.error("Failed to install Python dependencies:", error);
    throw error;
  }
}

setup()
  .then(() => {
    console.log("setup.ts complete");
    process.exit(0);
  })
  .catch((error) => {
    console.error(error);
    process.exit(1);
  });

```
