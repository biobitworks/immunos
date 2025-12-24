---
source: /Users/byron/projects/rockbeatspaper/src/server/trpc/procedures/runPiBenfordAnalysis.ts
relative: rockbeatspaper/src/server/trpc/procedures/runPiBenfordAnalysis.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { z } from "zod";
import { TRPCError } from "@trpc/server";
import { spawn } from "child_process";
import { baseProcedure } from "~/server/trpc/main";
import path from "path";

export const runPiBenfordAnalysis = baseProcedure
  .input(
    z.object({
      numDigits: z.number().int().min(10).max(100000),
    })
  )
  .mutation(async ({ input }) => {
    return new Promise<{ output: string }>((resolve, reject) => {
      // Path to whatever.py (in project root)
      const scriptPath = path.join(process.cwd(), "whatever.py");
      
      // Spawn Python process
      const pythonProcess = spawn("python3", [scriptPath]);
      
      let stdout = "";
      let stderr = "";
      
      // Capture stdout
      pythonProcess.stdout.on("data", (data) => {
        stdout += data.toString();
      });
      
      // Capture stderr
      pythonProcess.stderr.on("data", (data) => {
        stderr += data.toString();
      });
      
      // Handle process completion
      pythonProcess.on("close", (code) => {
        if (code !== 0) {
          reject(
            new TRPCError({
              code: "INTERNAL_SERVER_ERROR",
              message: `Python script exited with code ${code}. Error: ${stderr}`,
            })
          );
        } else {
          resolve({ output: stdout });
        }
      });
      
      // Handle process errors
      pythonProcess.on("error", (error) => {
        reject(
          new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to start Python process: ${error.message}`,
          })
        );
      });
      
      // Send the number of digits as input (answer the script's prompt)
      pythonProcess.stdin.write(`${input.numDigits}\n`);
      pythonProcess.stdin.end();
    });
  });

```
